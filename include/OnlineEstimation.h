#pragma once

//=================================
// include guard
#ifndef ONLINEESTIMATION_CLASS_H
#define ONLINEESTIMATION_CLASS_H

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
class OnlineEstimation {
private:

	//index of the start of online estimation after several steps
	int n_oe_start;
	//index of the start of cyclic gait
	int c;

	//number of dataset for foot placement data store
	int n_fp;
	//number of dataset for velcoity estimation data store
	int n_vel;

	//cyclic_gait_start;
	double Zc;
	double Tc;

	// LIPM model setting;	
	double StepTime;
	double tao;
	double g;
	double com_pos_current;
	double com_vel_current;
	static bool check_initial;

public:

	double StanceFoot_global;
	//double vel_f_predict;
	double	return_fp_wrt_com;
	MatrixXd vel_f_predict;
	MatrixXd fun_return_fp_wrt_com;
	MatrixXd walking_state_current;
	MatrixXd walking_state_for_fp_predict;

	// define setting for data setting for online estimation 
	void Init_OE(const int & fun_n_oe_start, const int & fun_c ,const int & fun_n_fp, const int & fun_n_vel);
	void Init_LIPM(double &constant_height, double &steptime);


	// setting for velocity online estimation
	MatrixXd com_state_stack;
	MatrixXd com_state_s;
	MatrixXd vel_f_predict_s;
	MatrixXd coeff_vel_LIP;
	MatrixXd return_coeff_vel;
	MatrixXd P1_value;
	MatrixXd P1;
	double	 Gain_P1;
	MatrixXd Q1_value;
	MatrixXd Q1;
	double	 Gain_Q1;

	// setting for foot placement online estimation
	// declaration statement for data set 
	MatrixXd walking_state_stack;
	MatrixXd dxcom_s;
	MatrixXd fp_wrt_com;
	MatrixXd coeff_fp_LIP;
	MatrixXd return_coeff_fp;
	MatrixXd P2_value;
	MatrixXd P2;
	double	 Gain_P2;
	MatrixXd Q2_value;
	MatrixXd Q2;
	double	 Gain_Q2;
	
	double OnlineEsitmation_main(const double &fun_com, const double &fun_dcom,const double &fun_com_offset, const double &fun_input_fp_wrt_com, const double &target_vel ,int & stepindex);
	void OnlineEsitmation_fp(const double &fun_com, const double &fun_dcom,const double &fun_com_offset, const double &fp, int & stepindex);
	void OnlineEstimation_fp_weighting_matrix(const int & fun_n_fp);
	void OnlineEsitmation_vel(const double &fun_com, const double &fun_dcom, const double &fun_com_offset, int & stepindex);
	void OnlineEstimation_vel_weighting_matrix(const int & fun_n_vel);
	double vel_target(double time, double totaltime, double StepTime);

	/*
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

	// setting for foot placement online estimation

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
	*/

};

// For log out the ground turth data, not relate to online estimation 
// 0 for origin mode, LIMP + ZMP
// 1 for online estimation mode
extern const int control_mode;
//extern  int impact;
//extern  int impact_forget;
#endif