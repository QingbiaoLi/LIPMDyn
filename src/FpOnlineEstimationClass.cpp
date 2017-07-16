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


//http://web.stanford.edu/group/SOL/software/lsqr/

void FpOnlineEstimation::Init(const int & number_of_data_set) {
	// 	To organize  '1 to n ' to 'dataset_past_walking_vel'
	// 	'2 to n+1 ' to 'dataset_next_footplacement'
	//cout<<"test initialization"<<endl; 
	number_of_dataset = number_of_data_set;
	cyclic_gait_start = 4;
	threshold_CoM_pos = 0.001;
	threshold_CoM_Vel_s = 0.0005;
	//threshold_CoM_Vel_f= 0.0045;
	threshold_fp = 0.0005;

	//impact;
	dataset_past_walking_state_stack.resize(number_of_dataset + 1, 3);
	//for CoM offset, CoM Vel, CoM Pos, 
	dataset_past_walking_vel.resize(number_of_dataset, 3);
	//for next foot placement with respect to CoM position
	dataset_next_footplacement.resize(number_of_dataset, 1);
	//for current CoM velocity, targeted velocity, CoM offset;
	dataset_current_walking_state.resize(1, 3);
	weighting_matrix = MatrixXd::Identity(number_of_dataset, number_of_dataset);
	//weighting_matrix::setZero;
	statef.resize(2);
	//model_coeff.resize(3, 1);
	//model_coeff
	model_coeff = MatrixXd::Zero(3, 1);
	//model_coeff::setZero;
}

void FpOnlineEstimation::define_initial_LIPM_setting(double &constant_height, const double &local_com_pos, const double &local_com_vel, double &Predict_Time) {
	zc = constant_height;
	look_ahead_time = Predict_Time;
	com_pos_current = local_com_pos;
	com_vel_current = local_com_vel;
}

void FpOnlineEstimation::collect_past_walking_state_stack(const double &local_com_pos, const double &local_com_vel, const double &com_offset, int & StepIndex) {
	if (StepIndex>=cyclic_gait_start) { // start to collct data when robot enter cyclic gait
		if (StepIndex < (number_of_dataset + cyclic_gait_start + 1)) {
			int index_stack = StepIndex - cyclic_gait_start;
			//store into Xs matrix
			dataset_past_walking_state_stack(index_stack, 0) = local_com_pos;
			dataset_past_walking_state_stack(index_stack, 1) = local_com_vel;
			dataset_past_walking_state_stack(index_stack, 2) = com_offset;


			//store into X matrix
			if (index_stack>=1) {
				//for (int index_stack = 0; index_stack < (StepIndex-cyclic_gait_start-1); index_stack++) {
				//collect 1 to n walking state as past walking state (CoM Vel, CoM next vel acutral,1)
				dataset_past_walking_vel(index_stack - 1, 0) = dataset_past_walking_state_stack(index_stack - 1, 1);
				dataset_past_walking_vel(index_stack - 1, 1) = dataset_past_walking_state_stack(index_stack, 1);
				dataset_past_walking_vel(index_stack - 1, 2) = dataset_past_walking_state_stack(index_stack - 1, 2);
				//}
			}

		}
		else {
			// test if the data set is fixed.
			//---------------------------------------------------		
			//------------------data updating for stack--------------------

			for (int index_stack = 0; index_stack < number_of_dataset; index_stack++) {
				dataset_past_walking_state_stack(index_stack, 0) = dataset_past_walking_state_stack(index_stack + 1, 0);					
				dataset_past_walking_state_stack(index_stack, 1) = dataset_past_walking_state_stack(index_stack + 1, 1);
				dataset_past_walking_state_stack(index_stack, 2) = com_offset;
			}
			dataset_past_walking_state_stack(number_of_dataset, 0) = local_com_pos;
			dataset_past_walking_state_stack(number_of_dataset, 1) = local_com_vel;
			dataset_past_walking_state_stack(number_of_dataset, 2) = com_offset;

			//------------------data updating for walking state--------------------
			for (int index_stack = 0; index_stack < number_of_dataset; index_stack++) {
				//collect 1 to n walking state as past walking state (CoM Vel, CoM next vel acutral,1)
				dataset_past_walking_vel(index_stack, 0) = dataset_past_walking_state_stack(index_stack, 1);
				dataset_past_walking_vel(index_stack, 1) = dataset_past_walking_state_stack(index_stack + 1, 1);
				dataset_past_walking_vel(index_stack, 2) = dataset_past_walking_state_stack(index_stack, 2);
			}
			//---------------------------------------------------	
		}
	}
}

void FpOnlineEstimation::collect_walking_state(const double &local_com_pos, const double &local_com_vel, const double &com_offset, const double &next_fp_wrt_CoM, int & StepIndex) {
	if (StepIndex>=cyclic_gait_start) {// start to collct data when robot enter cyclic gait
	//cout << "----step index: " << StepIndex << "\t local xcom: " << local_com_pos << "\t dxcom: " << local_com_vel << "\t ft wrt CoM: " << next_fp_wrt_CoM << endl;

		if (StepIndex < number_of_dataset + cyclic_gait_start + 1) {
			// data collection during initial stage

			//collect Xs (CoM pos, CoM vel, 1) data from "c" to "c+n+1" 
			collect_past_walking_state_stack(local_com_pos, local_com_vel, com_offset, StepIndex);		
			//collect foot placement data from "c+1" to "c+n+1"
			collect_next_footplacement(next_fp_wrt_CoM, StepIndex);
	

		}
		else {
			//selection data of data under threshold
			int bool_CoM_vel_s = 1;
			int bool_CoM_vel_e = 1;
			int bool_fp = 1;
			
			for (int index_th = 0; index_th < number_of_dataset; index_th++) { 
				//threshold for CoM velocity at the beginning of touch down
				/*
				if (abs(local_com_vel - dataset_past_walking_vel(index_th, 0)) >= threshold_CoM_Vel_s) {
					//cout << "------out of velocity threshold, store in data --------" << endl;
					bool_CoM_vel_s *= 1;
				}
				else {
					//cout << "------within velocity threshold, ignore data --------" << endl;
					bool_CoM_vel_s *= 0;
				}
				//threshold for foot placment change
				
				if (abs(next_fp_wrt_CoM - dataset_next_footplacement(index_th, 0)) >= threshold_fp) {
					//cout << "------out of fp threshold, store in data--------" << endl;
					bool_fp *= 1;
				}
				else {
					//cout << "------within fp threshold, ignore data --------" << endl;
					bool_fp *= 0;
				}
				
				//threshold for CoM velocity at the end of touch down, desired velocity
				if(abs(local_com_vel- dataset_past_walking_state_stack(index_th, 1))>=threshold_CoM_Vel_f){
					//cout<<"------out of velocity threshold, store in data --------"<<endl;
				bool_CoM_vel_e *=1;
				}
				else{
					//cout<<"------within velocity threshold, ignore data --------"<<endl;
				bool_CoM_vel_e *=0;
				}
				*/
			}
			
			if (((bool_CoM_vel_s*bool_fp) == 1) || (impact == 1)) {
				//distribute the stack data into walking state
				collect_past_walking_state_stack(local_com_pos, local_com_vel, com_offset, StepIndex);
				collect_next_footplacement(next_fp_wrt_CoM, StepIndex);

				if (impact_forget < number_of_dataset) {
					impact = 1;
					impact_forget += 1;
				}
				else {
					impact = 0;

				}
				//cout << "add data" << endl;
			}
			else {
				impact_forget = 0;
			}
			/*
			cout<<"impact: "<<impact<<"impact_forget"<<impact_forget<<endl;
			cout<<"bool_com_vel: "<<bool_CoM_vel_s<<endl;
			cout<<"bool_fp: "<<bool_fp<<" fp: "<<next_fp_wrt_CoM<<endl;
			cout<<"threshold cond: "<<(bool_CoM_vel_s*bool_fp)<<endl<<endl;
			*/
			set_weighting_matrix(dataset_past_walking_vel);
		}
	}
	//cout<<"dataset_past_walking_state_stack"<<endl<<dataset_past_walking_state_stack<<endl;
	//cout<< "weighting_matrix"<<endl<<weighting_matrix<<endl;	
	//cout << "dataset_past_walking_vel" << endl << dataset_past_walking_vel << endl;
	//cout << "dataset_next_footplacement" << endl << dataset_next_footplacement << endl << endl;
}

void FpOnlineEstimation::collect_next_footplacement(const double &next_fp_wrt_CoM, int & StepIndex) {
	if (StepIndex>=cyclic_gait_start+1) { // start to collct data when robot enter cyclic gait
		if ((StepIndex < number_of_dataset + cyclic_gait_start+1)) {
			dataset_next_footplacement(StepIndex - (cyclic_gait_start + 1), 0) = next_fp_wrt_CoM;
		}
		else {
			// test if the data set is fixed.	
			//---------------------------------------------------		
			//------------------data updating--------------------	
			//if((abs(ft_com[0]- data_set_xH(number_of_dataset, 1))>= threshold_CoM_Pos_x )&& (abs(ft_com[1] - data_set_x_stack(number_of_dataset, 2))) >= threshold_CoM_Vel_x){
			for (int index_fp_stack = 0; index_fp_stack < (number_of_dataset - 1); index_fp_stack++) {
				dataset_next_footplacement(index_fp_stack, 0) = dataset_next_footplacement(index_fp_stack + 1, 0);
			}
			dataset_next_footplacement(number_of_dataset - 1, 0) = next_fp_wrt_CoM;
			//}		
		}
	}

}
void FpOnlineEstimation::set_weighting_matrix(const Ref<MatrixXd> &dataset_current) {
	//--------------------------identity matrix---------------------------
	
	for (int i = 0; i < number_of_dataset ; i++) {

	weighting_matrix(i,i)=1;

	}
	
	double average = dataset_current.col(0).mean();

	//--------------------------error variance---------------------------
//weight proportional to error variance, based on CoM vel at the beginning of touch-down
	/*
	for (int i = 0; i < number_of_dataset ; i++) {
	weighting_matrix(i,i)=pow((dataset_current(i,0)-average),2);
	}
	//normalized
	double sum = weighting_matrix.sum();
	for (int i = 0; i < number_of_dataset ; i++) {
	weighting_matrix(i,i)=weighting_matrix(i,i)/sum;
	}
	*/
//weight inversely proportional to error variance
	/*
	for (int i = 0; i < number_of_dataset ; i++) {
	weighting_matrix(i,i)=1/pow((dataset_current(i,0)-average),2);
	}
	//cout<<"weighting_matrix: "<<endl<<weighting_matrix<<endl;
	double sum = weighting_matrix.sum();
	//normalized
	for (int i = 0; i < number_of_dataset ; i++) {
	weighting_matrix(i,i)=weighting_matrix(i,i)/sum;
	}
	*/
	/*
	//-------------------time dependent-------------------------

	for (int i = 1; i <= number_of_dataset; i++) {
		double sum = (1 + number_of_dataset)*number_of_dataset / 2;
		weighting_matrix(i - 1, i - 1) = i / sum;
	}
	*/
	//time dependent-- exp(kt),t is consequece of step,0-older, 6 is last	
	/*
	for (int i = 0; i < number_of_dataset; i++) {
	weighting_matrix(i,i)=exp(i);
	}
	double sum = weighting_matrix.sum();
	for (int i = 0; i < number_of_dataset; i++) {
	weighting_matrix(i,i)=weighting_matrix(i,i)/sum;
	}
	
	cout << "weighting_matrix: " << endl << weighting_matrix << endl;
	*/
}
void FpOnlineEstimation::collect_current_walking_state(const double &local_com_vel, const double &targeted_vel) {
	dataset_current_walking_state(0, 0) = local_com_vel;
	dataset_current_walking_state(0, 1) = targeted_vel;//-copysign(dataset_past_walking_state_stack(number_of_dataset, 1),abs(dataset_past_walking_state_stack(number_of_dataset-1, 1)));//targeted_vel;
	dataset_current_walking_state(0, 2) = 1;//dataset_past_walking_state_stack(number_of_dataset, 2);// CoM offset--1
																								 //http://www.cplusplus.com/reference/cmath/copysign/
}

double FpOnlineEstimation::estimate_walking_state_next_step(const Ref<MatrixXd> &current_walking_state, const Ref<MatrixXd> &model_coeff) {
	cout<<"current_walking_state"<<endl<<current_walking_state<<endl;

	walking_state_next_step = current_walking_state * model_coeff;
	footplacement_predict = walking_state_next_step(0, 0);
	cout << "Online_estimation_fp: " << footplacement_predict << "; " << endl << endl;
	return footplacement_predict;
}

MatrixXd FpOnlineEstimation::calculate_model_coeff_fp(const Ref<MatrixXd> &dataset_current, const Ref<MatrixXd> &dataset_next_step) {
	model_coeff = pinv(weighting_matrix*dataset_current)*weighting_matrix*dataset_next_step;

	//cout<<"model_coeff"<<endl<<model_coeff<<endl;
	return model_coeff;
}


//----------------------------------------------
MatrixXd FpOnlineEstimation::calculate_model_coeff_dxf(const Ref<MatrixXd> &dataset_current, const Ref<MatrixXd> &dataset_next_step) {
	model_coeff = pinv(weighting_matrix*dataset_current)*weighting_matrix*dataset_next_step;

	//cout<<"model_coeff"<<endl<<model_coeff<<endl;
	return model_coeff;
}


void FpOnlineEstimation::StateEsimation(const double &local_com_vel, const double & targeted_vel, const int & StepIndex) {
	//targeted_velocity = targeted_vel;
	double zc = 0.5;
	double g = 9.81;
	double Tc = sqrt(zc / g);
	double StepTime1 = 0.5;
	double tao = StepTime1 / Tc;
	if (StepIndex < (number_of_dataset + cyclic_gait_start + 1)) {
		static const auto runOnce = [] { cout << "-------------This is fp control: Initial LIPM to collect data----------------" << endl; return true; }();
		//LIPM_StateEsimation(zc,com_pos_current,com_vel_current,look_ahead_time);
		model_coeff(0, 0) = Tc*cosh(tao)/sinh(tao);
		model_coeff(1, 0) = -Tc/ sinh(tao);
		model_coeff(2, 0) = 0;
		collect_current_walking_state(local_com_vel, targeted_vel);
		estimate_walking_state_next_step(dataset_current_walking_state, model_coeff);
	}
	/*
	else if(StepIndex <= (number_of_dataset+5)){
	static const auto runOnce = [] { cout << "-------------Start online estimation----------------"<< endl; return true;}();
	//LIPM_StateEsimation(zc,com_pos_current,com_vel_current,look_ahead_time);
	collect_current_walking_state(local_com_vel,targeted_vel);
	calculate_model_coeff_fp(dataset_past_walking_vel,dataset_next_footplacement);
	//estimate_walking_state_next_step(dataset_current_walking_state,model_coeff);
	}
	*/
	else {
		static const auto runOnce = [] { cout << "-------------Start online estimation----------------" << endl; return true; }();
		collect_current_walking_state(local_com_vel, targeted_vel);
		calculate_model_coeff_fp(dataset_past_walking_vel, dataset_next_footplacement);
		/*
		model_coeff(0, 0) = Tc*cosh(tao) / sinh(tao);
		model_coeff(1, 0) = -Tc / sinh(tao);
		model_coeff(2, 0) = 0;
		*/
		estimate_walking_state_next_step(dataset_current_walking_state, model_coeff);
		
	}
}

double  FpOnlineEstimation::vel_target(double time, double totaltime, double StepTime) {
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
		dx = 0.2 / (2 * StepTime)*(time- totaltime / 3) + 0.1;
	}
	else if (time< totaltime * 2 / 3) {
		dx = 0.3;
	}
	else if (time < t3) {
		dx =  0.2 / (2 * StepTime)*(time - totaltime*2 / 3) + 0.3;
	}
	else {
		dx = 0.5;
	}
	return dx;
}
