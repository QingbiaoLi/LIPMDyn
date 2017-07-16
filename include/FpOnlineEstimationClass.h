//=================================
// include guard
#ifndef FP_ONLINEESTIMATION_CLASS_H
#define FP_ONLINEESTIMATION_CLASS_H

//=================================
// included dependencies
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h> 
//#include "RTControlClass.h"

using namespace std;
using namespace Eigen;
//=================================
// the actual class
class FpOnlineEstimation {
private:
	int number_of_dataset;
	int cyclic_gait_start;
	double zc;
	double look_ahead_time;
	double g;
	double com_pos_current;
	double com_vel_current;
	static bool check_initial;

public:
	// define setting for data setting for online estimation 
	void Init(const int & number_of_data_set);

	// declaration statement for data set 
	MatrixXd model_coeff;
	MatrixXd dataset_past_walking_vel;
MatrixXd dataset_past_walking_state;// for projection from dx0 to dxf
	MatrixXd dataset_next_footplacement;
	MatrixXd dataset_current_walking_state;

	MatrixXd dataset_past_walking_state_stack;
	MatrixXd weighting_matrix;
	// double targeted_velocity;
	// declaration for threshold
	double threshold_CoM_pos;
	double threshold_CoM_Vel_s;
	double threshold_CoM_Vel_f;
	double threshold_fp;
	int impact;
	int impact_forget;

	MatrixXd walking_state_next_step;
	//vector<double> walking_state_next_step;
	vector<double> statef;
	double footplacement_predict;

	// include store CoM state into "dataset_past_walking_state_stack", and store them into "dataset_past_walking_state"
	void collect_walking_state(const double &local_com_pos, const double &local_com_vel, const double &com_offset, const double &next_fp_wrt_CoM, int & StepIndex);
	void collect_past_walking_state_stack(const double &local_com_pos, const double &local_com_vel, const double &com_offset, int & StepIndex);
	void collect_next_footplacement(const double &next_fp_wrt_CoM, int & StepIndex);
	void collect_current_walking_state(const double &local_com_vel, const double &targeted_vel);

	void set_weighting_matrix(const Ref<MatrixXd> &dataset_current);

MatrixXd calculate_model_coeff_dxf(const Ref<MatrixXd> &dataset_current, const Ref<MatrixXd> &dataset_next_step);
	MatrixXd calculate_model_coeff_fp(const Ref<MatrixXd> &dataset_current, const Ref<MatrixXd> &dataset_next_step);
	double estimate_walking_state_next_step(const Ref<MatrixXd> &current_walking_state, const Ref<MatrixXd> &model_coeff);

	void define_initial_LIPM_setting(double &constant_height, const double &local_com_pos, const double &local_com_vel, double &Predict_Time);
	//vector<double> LIPM_StateEsimation(double &constant_height, const double &local_com_pos, const double &local_com_vel, double &Predict_Time);
	double vel_target(double time, double totaltime, double StepTime);
	void StateEsimation(const double &local_com_vel, const double & targeted_vel, const int & StepIndex);

};

// For log out the ground turth data, not relate to online estimation 
// 0 for origin mode, LIMP + ZMP
// 1 for online estimation mode
extern const int control_mode;
//extern  int impact;
//extern  int impact_forget;
#endif