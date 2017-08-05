#include "OnlineEstimation.h"

template <class MatT>
Eigen::Matrix<typename MatT::Scalar, MatT::ColsAtCompileTime, MatT::RowsAtCompileTime> pinv(const MatT &mat, typename MatT::Scalar tolerance = typename MatT::Scalar{ 1e-4 }) // choose appropriately
{
	typedef typename MatT::Scalar Scalar;
	auto svd = mat.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
	const auto &singularValues = svd.singularValues();
	Eigen::Matrix<Scalar, MatT::ColsAtCompileTime, MatT::RowsAtCompileTime> singularValuesInv(mat.cols(), mat.rows());
	singularValuesInv.setZero();
	for (int i = 0; i < singularValues.size(); ++i) {
		if (singularValues(i) > tolerance)
		{
			singularValuesInv(i, i) = Scalar{ 1 } / singularValues(i);
		}
		else
		{
			singularValuesInv(i, i) = Scalar{ 0 };
		}
	}
	return svd.matrixV() * singularValuesInv * svd.matrixU().adjoint();
}



void OnlineEstimation::Init_OE(const int & fun_n_oe_start, const int & fun_c, const int & fun_n_fp, const int & fun_n_vel) {

	//cout<<"test initialization"<<endl; 
	n_oe_start = fun_n_oe_start;
	//index of the start of cyclic gait
	c = fun_c;

	//number of dataset for foot placement data store
	n_fp = fun_n_fp;
	//number of dataset for velcoity estimation data store
	n_vel = fun_n_vel;

	// Setting for foot placement online estimation
	// To organize  '1 to n ' to 'dxcom_s', '2 to n+1 ' to 'fp_wrt_com'
	walking_state_stack.resize(fun_n_fp + 1, 3);
	//for CoM offset, CoM Vel, CoM Pos, 
	dxcom_s.resize(fun_n_fp, 3);
	//for next foot placement with respect to CoM position
	fp_wrt_com.resize(fun_n_fp, 1);
	// Model coefficient of foot placement estimation in LIPM
	coeff_fp_LIP.resize(3, 1);
	OnlineEstimation_fp_weighting_matrix(fun_n_fp);

	// setting for foot placement online estimation
	com_state_stack.resize(fun_n_vel + 1, 3);
	//for CoM offset, CoM Vel, CoM Pos, 
	com_state_s.resize(fun_n_vel, 3);
	// Model coefficient of foot placement estimation in LIPM
	vel_f_predict_s.resize(fun_n_vel, 1);
	// Model coefficient of velocity estimation in LIPM
	coeff_vel_LIP.resize(3, 1);
	OnlineEstimation_fp_weighting_matrix(fun_n_vel);

	//obtain the current walking state for velocity estimation
	walking_state_current.resize(1, 3);
	//Estimation of final velocity(LIPM)  continuous phase
	walking_state_for_fp_predict.resize(1, 3);
	/*
	//threshold_CoM_pos = 0.001;
	//threshold_CoM_Vel_s = 0.0005;
	//threshold_CoM_Vel_f= 0.0045;
	//threshold_fp = 0.0005;
	//for current CoM velocity, targeted velocity, CoM offset;
	dataset_current_walking_state.resize(1, 3);
	weighting_matrix = MatrixXd::Identity(number_of_dataset, number_of_dataset);
	//weighting_matrix::setZero;
	statef.resize(2);
	//model_coeff.resize(3, 1);
	//model_coeff
	model_coeff = MatrixXd::Zero(3, 1);
	//model_coeff::setZero;
	*/
}

void OnlineEstimation::Init_LIPM(double &constant_height, double &steptime){
	g = 9.81;
	Zc	= constant_height;
	StepTime =	steptime;
	Tc = sqrt(Zc / g);
	tao = StepTime / Tc;
	// Model coefficient of foot placement estimation in LIPM
	coeff_fp_LIP(0, 0) = cosh(tao)/sinh(tao)*Tc;
	coeff_fp_LIP(1, 0) = -1/sinh(tao)*Tc;
	coeff_fp_LIP(2, 0) = 0;
	// Model coefficient of velocity estimation in LIPM
	coeff_vel_LIP(0, 0) = sinh(tao)/Tc;
	coeff_vel_LIP(1, 0) = cosh(tao);
	coeff_vel_LIP(2, 0) = 0;
}


double OnlineEstimation::OnlineEsitmation_main(const double &fun_com, const double &fun_dcom, const double &fun_com_offset, const double &fun_input_fp_wrt_com, const double &target_vel, int & stepindex) {
	//Place the estimation of swing foot location, then this become the global support foot location
	StanceFoot_global = fun_com + fun_input_fp_wrt_com;

	//----------------------------- Store Data ----------------------

	//use online estimation to obtain the model coefficients for veolcity estimation, continuous phase
	OnlineEsitmation_vel(fun_com - StanceFoot_global, fun_dcom, fun_com_offset, stepindex);

	//use online estimation to obtain the model coefficients for foot placement estimation, discrete transitions
	OnlineEsitmation_fp(fun_com - StanceFoot_global, fun_dcom, fun_com_offset, fun_input_fp_wrt_com, stepindex);

	//obtain the current walking state for velocity estimation
	walking_state_current(0, 0) = fun_com - StanceFoot_global;
	walking_state_current(0, 1) = fun_dcom;
	walking_state_current(0, 2) = 1;


	//----------------------------- Online Estimation ----------------------

	//Estimation of final velocity(LIPM)  continuous phase
	vel_f_predict = walking_state_current * return_coeff_vel;
	//print("vel_f_predict", vel_f_predict)
	walking_state_for_fp_predict(0, 0) = vel_f_predict(0,0);
	walking_state_for_fp_predict(0, 1) = target_vel;
	walking_state_for_fp_predict(0, 2) = 1;
	//next placement of swing foot(LIPM)  discrete transitions
	fun_return_fp_wrt_com = walking_state_for_fp_predict * coeff_fp_LIP;
	return_fp_wrt_com = fun_return_fp_wrt_com(0, 0);

	return return_fp_wrt_com;
}

void OnlineEstimation::OnlineEstimation_vel_weighting_matrix(const int & fun_n_vel) {
	// velocity regression term
	// P2 will set higher weight to lastest step
	P1_value.resize(fun_n_vel, 1);
	P1.resize(fun_n_vel, fun_n_vel);
	Q1_value.resize(3, 1);
	Q1.resize(3, 3);
	//gain_P is used to weight the influence of the regression term
	Gain_P1 = 0.1;
	//gain_Q is used to weight the influence of the regularisation term
	Gain_Q1 = 0.1;
	P1_value.resize(fun_n_vel, 1);
	P1.resize(fun_n_vel, fun_n_vel);
	Q1_value.resize(3, 1);
	Q1.resize(3, 3);
	for (int i = 0; i < fun_n_vel; i++) {
		P1_value(i, 0) = i + 1;
	}
	Q1_value(0, 0) = 1;
	Q1_value(1, 0) = 1;
	Q1_value(2, 0) = 0.01;
	P1 = P1_value.array().matrix().asDiagonal()*Gain_P1;
	Q1 = Q1_value.array().matrix().asDiagonal()*Gain_Q1;

}

void OnlineEstimation::OnlineEstimation_fp_weighting_matrix(const int & fun_n_fp) {
	// foot placement regression term
	// P1 will set higher weight to lastest step
	P2_value.resize(fun_n_fp, 1);
	P2.resize(fun_n_fp, fun_n_fp);
	Q2_value.resize(3, 1);
	Q2.resize(3, 3);
	//gain_P is used to weight the influence of the regression term
	Gain_P2 = 0.1;
	//gain_Q is used to weight the influence of the regularisation term
	Gain_Q2 = 0.1;
	for (int i = 0; i < n_fp; i++) {
		P2_value(i, 0) = i + 1;
	}
	Q2_value(0, 0) = 1;
	Q2_value(1, 0) = 1;
	Q2_value(2, 0) = 0.01;
	P2 = P2_value.array().matrix().asDiagonal()*Gain_P2;
	Q2 = Q2_value.array().matrix().asDiagonal()*Gain_Q2;
}


void OnlineEstimation::OnlineEsitmation_vel(const double &fun_com, const double &fun_dcom, const double &fun_com_offset, int & stepindex) {
	// data collection, based on Fisrt in First Out
	if (stepindex >= c){

		if (stepindex <= (n_vel + c)){
			//print("n_vel", n_vel)
			//print("n_vel+c", n_vel + c)

			int index = stepindex - c;
			//print("index", index)
			com_state_stack(index, 0) = fun_com;
			com_state_stack(index, 1) = fun_dcom;
			com_state_stack(index, 2) = fun_com_offset;

			if(index >= 1){
				com_state_s(index-1, 0) = com_state_stack(index - 1, 0);
				com_state_s(index-1, 1) = com_state_stack(index - 1, 1);
				com_state_s(index-1, 2) = com_state_stack(index - 1, 2);
				vel_f_predict_s(index - 1, 0) = com_state_stack(index, 1);
			}
		}
	}

	else{

		for (int i = 0; i < n_vel; i++) {
			com_state_stack(i, 0) = com_state_stack(i + 1, 0);
			com_state_stack(i, 1) = com_state_stack(i + 1, 1);
			com_state_stack(i, 2) = com_state_stack(i + 1, 2);
		}
		com_state_stack(n_vel, 0) = fun_com;
		com_state_stack(n_vel, 1) = fun_dcom;
		com_state_stack(n_vel, 2) = fun_com_offset;

		for (int i = 0; i < n_vel; i++) {
			com_state_s(i - 1, 0) = com_state_stack(i, 0);
			com_state_s(i - 1, 1) = com_state_stack(i, 1);
			com_state_s(i - 1, 2) = com_state_stack(i, 2);
			vel_f_predict_s(i, 0) = com_state_stack(i + 1, 1);
		}
	}

	if (stepindex < (n_oe_start + c)) {
		cout << "--------------" << endl;
		return_coeff_vel = coeff_vel_LIP;
	}
	else{
		cout << "-------online estimation start-------" << endl;
		//Tikhonov regularization
		//https: // en.wikipedia.org / wiki / Tikhonov_regularization  # Generalized_Tikhonov_regularization
		return_coeff_vel = coeff_vel_LIP + pinv(com_state_s.transpose()*P1*com_state_s + Q1)*(com_state_s.transpose()*P1*(vel_f_predict_s - com_state_s*coeff_vel_LIP));
	
	}

}

void OnlineEstimation::OnlineEsitmation_fp(const double &fun_com, const double &fun_dcom, const double &fun_com_offset, const double &fp, int & stepindex) {
	if (stepindex >= c) {
		if (stepindex <= (n_fp + c)) {
			int index = stepindex - c;
			walking_state_stack(index, 0) = fun_com;
			walking_state_stack(index, 1) = fun_dcom;
			walking_state_stack(index, 2) = fun_com_offset;
			if (index >= 1) {
				dxcom_s(index - 1, 0) = walking_state_stack(index - 1, 1);
				dxcom_s(index - 1, 1) = walking_state_stack(index, 1);
				dxcom_s(index - 1, 2) = fun_com_offset;

			}
		}
		else {

			for (int i = 0; i < n_fp; i++) {
				walking_state_stack(n_fp, 0) = walking_state_stack(i + 1, 0);
				walking_state_stack(n_fp, 1) = walking_state_stack(i + 1, 1);
				walking_state_stack(n_fp, 2) = walking_state_stack(i + 1, 2);
			}
			walking_state_stack(n_fp, 0) = fun_com;
			walking_state_stack(n_fp, 1) = fun_dcom;
			walking_state_stack(n_fp, 2) = fun_com_offset;

			for (int i = 0; i < n_fp; i++) {
				dxcom_s(i, 0) = walking_state_stack(i, 1);
				dxcom_s(i, 1) = walking_state_stack(i + 1, 1);
				dxcom_s(i, 2) = walking_state_stack(i + 1, 2);
			}
		}
	}
	if (stepindex >= (c)) {

		if (stepindex < (n_fp + c)) {

			cout << (n_fp + c) << endl;
			cout << (stepindex) << endl;
			cout << ("fp", fp) << endl;
			int index = stepindex - c;
			fp_wrt_com(index, 0) = fp;
		}
		else {
			cout << stepindex << endl;
			for (int i = 0; i < (n_fp - 1); i++) {
				fp_wrt_com(i, 0) = fp_wrt_com(i + 1, 0);

				fp_wrt_com(n_fp - 1, 0) = fp;
			}
		}
	}

	if (stepindex < (n_oe_start + c)) {
		cout << "--------------" << endl;
		return_coeff_vel = coeff_vel_LIP;
	}
	else {
		cout << "-------online estimation start-------" << endl;
		//Tikhonov regularization
		//https: // en.wikipedia.org / wiki / Tikhonov_regularization  # Generalized_Tikhonov_regularization
		return_coeff_vel = coeff_fp_LIP + pinv((dxcom_s.transpose()*P2*dxcom_s + Q2))*(dxcom_s.transpose()*P2*(fp_wrt_com - dxcom_s*coeff_fp_LIP));
	}

}

double  OnlineEstimation::vel_target(double time, double totaltime, double StepTime) {
	double dx;
	double step2target_vel = 2;
	double t1 = 2 * StepTime;
	double t2 = 2 * StepTime + totaltime / 3;
	double t3 = 2 * StepTime + totaltime * 2 / 3;
	if (time < t1) {// transist time to reach 0.1 m/s
		dx = 0.1 / t1*time;
	}
	else if (time< totaltime / 3) {
		dx = 0.1;
	}
	else if (time< t2) {
		dx = 0.2 / (2 * StepTime)*(time - totaltime / 3) + 0.1;
	}
	else if (time< totaltime * 2 / 3) {
		dx = 0.3;
	}
	else if (time < t3) {
		dx = 0.2 / (2 * StepTime)*(time - totaltime * 2 / 3) + 0.3;
	}
	else {
		dx = 0.5;
	}
	return dx;
}
