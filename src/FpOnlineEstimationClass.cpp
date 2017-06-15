#include "FpOnlineEstimationClass.h"
// online estimation for foot placement 
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


void FpOnlineEstimation::Init(const int & number_of_data_set) {
	// 	To organize  '1 to n ' to 'dataset_past_walking_state'
	// 	'2 to n+1 ' to 'dataset_next_footplacement' 
	number_of_dataset = number_of_data_set;
	dataset_past_walking_state_stack.resize(number_of_dataset + 1, 2);
	//for CoM offset, CoM Vel, CoM Pos, 
	dataset_past_walking_state.resize(number_of_dataset, 2);
	//for next foot placement with respect to CoM position
	dataset_next_footplacement.resize(number_of_dataset, 1);
	//for current CoM velocity, targeted velocity, CoM offset;
	dataset_current_walking_state.resize(1, 2);
	statef.resize(2);
	threshold_CoM_Pos = 7e-4;
	threshold_CoM_Vel = 3e-3;
	control_coefficient.resize(2, 1);

}

void FpOnlineEstimation::define_initial_LIPM_setting(double &constant_height, const double &local_com_pos, const double &local_com_vel, double &Predict_Time) {
	zc = constant_height;
	look_ahead_time = Predict_Time;
	com_pos_current = local_com_pos;
	com_vel_current = local_com_vel;
}

void FpOnlineEstimation::collect_past_walking_state_stack(const double &local_com_pos, const double &local_com_vel, const double &com_offset, int & StepIndex) {
	if (StepIndex>1) { // start to collct data when robot enter cyclic gait
		if ((StepIndex <= number_of_dataset + 2)) {
			dataset_past_walking_state_stack(StepIndex - 2, 0) = local_com_pos;
			dataset_past_walking_state_stack(StepIndex - 2, 1) = local_com_vel;
			//dataset_past_walking_state_stack(StepIndex-2, 2) = com_offset;
		}
		else {

			// test if the data set is fixed.
			//---------------------------------------------------		
			//------------------data updating--------------------

			//if((abs(ft_com[0]- data_set_xH(number_of_dataset, 1))>= threshold_CoM_Pos_x )&& (abs(ft_com[1] - data_set_x_stack(number_of_dataset, 2))) >= threshold_CoM_Vel_x){
			for (int index_stack = 0; index_stack < number_of_dataset; index_stack++) {
				dataset_past_walking_state_stack(index_stack, 0) = dataset_past_walking_state_stack(index_stack + 1, 0);
				dataset_past_walking_state_stack(index_stack, 1) = dataset_past_walking_state_stack(index_stack + 1, 1);
				//dataset_past_walking_state_stack(index_stack, 2) = com_offset;
			}
			dataset_past_walking_state_stack(number_of_dataset, 0) = local_com_pos;
			dataset_past_walking_state_stack(number_of_dataset, 1) = local_com_vel;
			//dataset_past_walking_state_stack(number_of_dataset, 2) = com_offset;
			//}
			//---------------------------------------------------	
		}
	}
}

void FpOnlineEstimation::collect_walking_state(const double &local_com_pos, const double &local_com_vel, const double &com_offset, int & StepIndex) {

	collect_past_walking_state_stack(local_com_pos, local_com_vel, com_offset, StepIndex);
	//cout<<"dataset_past_walking_state_stack"<<endl<<dataset_past_walking_state_stack<<endl;

	if (StepIndex >= number_of_dataset + 2) {

		for (int index_stack = 0; index_stack < number_of_dataset; index_stack++) {
			//collect 1 to n walking state as past walking state (CoM Vel, CoM next vel acutral,1)
			dataset_past_walking_state(index_stack, 0) = dataset_past_walking_state_stack(index_stack, 1);
			dataset_past_walking_state(index_stack, 1) = dataset_past_walking_state_stack(index_stack + 1, 1);
			//dataset_past_walking_state(index_stack, 2) = dataset_past_walking_state_stack(index_stack, 2);
		}

	}
	//cout<<"dataset_past_walking_state"<<endl<<dataset_past_walking_state<<endl;
	//cout<<"dataset_next_footplacement"<<endl<<dataset_next_footplacement<<endl;
}

void FpOnlineEstimation::collect_next_footplacement(const double &next_fp_wrt_CoM, int & StepIndex) {
	if (StepIndex>1) { // start to collct data when robot enter cyclic gait
		if ((StepIndex <= number_of_dataset + 1)) {
			dataset_next_footplacement(StepIndex - 2, 0) = next_fp_wrt_CoM;
		}
		else {
			// test if the data set is fixed.	

			//---------------------------------------------------		
			//------------------data updating--------------------	
			//if((abs(ft_com[0]- data_set_xH(number_of_dataset, 1))>= threshold_CoM_Pos_x )&& (abs(ft_com[1] - data_set_x_stack(number_of_dataset, 2))) >= threshold_CoM_Vel_x){
			for (int index_stack = 0; index_stack < (number_of_dataset - 1); index_stack++) {
				dataset_next_footplacement(index_stack, 0) = dataset_next_footplacement(index_stack + 1, 0);
			}
			dataset_next_footplacement(number_of_dataset - 1, 0) = next_fp_wrt_CoM;
			//}		
		}
	}

}

void FpOnlineEstimation::collect_current_walking_state(const double &local_com_vel, const double &targeted_vel) {
	dataset_current_walking_state(0, 0) = local_com_vel;
	dataset_current_walking_state(0, 1) = targeted_vel;//-copysign(dataset_past_walking_state_stack(number_of_dataset, 1),abs(dataset_past_walking_state_stack(number_of_dataset-1, 1)));//targeted_vel;
													   //dataset_current_walking_state(0, 2) = dataset_past_walking_state_stack(number_of_dataset, 2);// CoM offset--1
													   //http://www.cplusplus.com/reference/cmath/copysign/
}

double FpOnlineEstimation::estimate_walking_state_next_step(const Ref<MatrixXd> &current_walking_state, const Ref<MatrixXd> &control_coefficient) {

	walking_state_next_step = current_walking_state * control_coefficient;
	footplacement_predict = walking_state_next_step(0, 0);
	//cout<<"current_walking_state"<<endl<<current_walking_state<<endl<<endl;
	//<<"control_coefficient"<<endl<<control_coefficient<<endl<<endl;
	//cout <<"Online_estimation_fp: "<<footplacement_predict<<"; ";
	//statef[0]=0;
	return footplacement_predict;
}

MatrixXd FpOnlineEstimation::calculate_control_coefficient(const Ref<MatrixXd> &dataset_current, const Ref<MatrixXd> &dataset_next_step) {
	control_coefficient = pinv(dataset_current)*dataset_next_step;

	//cout << "control_coefficient" << endl << control_coefficient << endl;

	return control_coefficient;
}

void FpOnlineEstimation::StateEsimation(const double &local_com_vel, const double & targeted_vel, const int & StepIndex) {
	//targeted_velocity = targeted_vel;
	if (StepIndex <= (number_of_dataset + 2)) {
		static const auto runOnce = [] { cout << "-------------This is fp control: Initial LIPM to collect data----------------" << endl; return true; }();
		//LIPM_StateEsimation(zc,com_pos_current,com_vel_current,look_ahead_time);
	}
	/*
	else if(StepIndex <= (number_of_dataset+5)){
	static const auto runOnce = [] { cout << "-------------Start online estimation----------------"<< endl; return true;}();
	//LIPM_StateEsimation(zc,com_pos_current,com_vel_current,look_ahead_time);
	collect_current_walking_state(local_com_vel,targeted_vel);
	calculate_control_coefficient(dataset_past_walking_state,dataset_next_footplacement);
	//estimate_walking_state_next_step(dataset_current_walking_state,control_coefficient);
	}
	*/
	else {
		static const auto runOnce = [] { cout << "-------------Start online estimation----------------" << endl; return true; }();
		collect_current_walking_state(local_com_vel, targeted_vel);
		calculate_control_coefficient(dataset_past_walking_state, dataset_next_footplacement);
		estimate_walking_state_next_step(dataset_current_walking_state, control_coefficient);
		//cout << "---online estimation----" << endl << "control coeff:" << control_coefficient << endl << "footplacement estimation: " << footplacement_predict << endl << endl;

	}

}

/*
vector<double> FpOnlineEstimation::LIPM_StateEsimation(double &constant_height,const double &local_com_pos,const double &local_com_vel, double &Predict_Time){
if (time<0){
cout<<"time should not be negative"<<endl;
return statef;
}
// it applies to all cases, no matter diverging or converging
double Tc, xf, dxf;
double g =9.81;
Tc=sqrt(zc/g);
double x0 	= local_com_pos;
double dx0 	= local_com_vel;
control_coefficient(0, 0) = 0.02;
control_coefficient(1, 0) =  sinh(Predict_Time/ Tc)/ Tc;// Pos
control_coefficient(2, 0) =  cosh(Predict_Time/ Tc); // Velocity

xf=x0*cosh(Predict_Time/Tc)+dx0*Tc*sinh(Predict_Time/Tc);
dxf=x0/Tc*sinh(Predict_Time/Tc)+dx0*cosh(Predict_Time/Tc);
//cout<<"predic ahead time"<<endl;
//cout<<Predict_Time<<endl;
//cout<<Tc*sinh(Predict_Time/Tc)/Tc<<endl;
//cout<<cosh(Predict_Time/ Tc)<<endl;
statef[0]=xf;
statef[1]=dxf;
//cout<<"desired com velocity"<<dxf<<endl;
return statef;
}
*/