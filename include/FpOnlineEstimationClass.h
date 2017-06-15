#pragma once
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
	MatrixXd control_coefficient;
	MatrixXd dataset_past_walking_state;
	MatrixXd dataset_next_footplacement;
	MatrixXd dataset_current_walking_state;
	MatrixXd dataset_past_walking_state_stack;

	// double targeted_velocity;
	double threshold_CoM_Pos;
	double threshold_CoM_Vel;

	MatrixXd walking_state_next_step;
	//vector<double> walking_state_next_step;
	vector<double> statef;
	double footplacement_predict;

	// include store CoM state into "dataset_past_walking_state_stack", and store them into "dataset_past_walking_state"
	void collect_walking_state(const double &local_com_pos, const double &local_com_vel, const double &com_offset, int & StepIndex);
	void collect_past_walking_state_stack(const double &local_com_pos, const double &local_com_vel, const double &com_offset, int & StepIndex);
	void collect_next_footplacement(const double &next_fp_wrt_CoM, int & StepIndex);

	void collect_current_walking_state(const double &local_com_vel, const double &targeted_vel);

	MatrixXd calculate_control_coefficient(const Ref<MatrixXd> &dataset_current, const Ref<MatrixXd> &dataset_next_step);

	double estimate_walking_state_next_step(const Ref<MatrixXd> &current_walking_state, const Ref<MatrixXd> &control_coefficient);

	void define_initial_LIPM_setting(double &constant_height, const double &local_com_pos, const double &local_com_vel, double &Predict_Time);
	//vector<double> LIPM_StateEsimation(double &constant_height, const double &local_com_pos, const double &local_com_vel, double &Predict_Time);

	void StateEsimation(const double &local_com_vel, const double & targeted_vel, const int & StepIndex);

};

// For log out the ground turth data, not relate to online estimation 
// 0 for origin mode, LIMP + ZMP
// 1 for online estimation mode
extern const int control_mode;

#endif