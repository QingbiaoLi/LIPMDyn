/* Programmed by: Zhibin Li*/
/*Modified by: Qinbiao Li*/
// this program tests the LIPM dynamcis functions
// 24 April 2014
// updated: 29 May 2014
#include "windows.h"
#include <fstream>
#include <iostream>
#include <chrono>
#include <random>
#include "LIPMDynClass.h"
#include "OnlineEstimation.h"
#include "FilterClass.h"
#include "RobotParameter.h"

double Poly(double t0, double x0, double tf, double xf, double time);
double Poly_stance(double x0, double xf, double t0, double tf, double time);

void IK2D(double theta, double r, double *angle);
void traj_angle_attack_init();
void traj_angle_attack_update();
void traj_angle_hipknee_generate();
void footjudgement();

int main() {
	StartTime = 0.8; // Duration of inital step strategy. (time of initiation of sagittal gait)
	xcom = 0.0; // Initial CoM position in global frame.
	dxcom = 0.0; // Initial CoM velocity in global frame.
	ddxcom = 0.0; // Initial CoM acceleration in global frame.	
	xcop2 = 0; // Inital predicted foot placement in global frame.

	xcom_e = 0.0;
	dxcom_e = 0.0;
	ddxcom_e = 0.0;

	xcom_e_f = 0.0;
	dxcom_e_f = 0.0;
	ddxcom_e_f = 0.0;

	//setting for filter

	FilterClass_xcom.butterworth(samplingTime, 10.0, 1);
	FilterClass_dxcom.butterworth(samplingTime, 10.0, 1);

	tswitch = LIPM.LateralSwitchTime(zc, xcom, dxcom, yf, p_plus, p_minus, StartTime); // return time of switching phase
	StepEndTime = StartTime + StepTime1; // End time of 1st periodic step.
	
	//online estimation setting
	//Init_OE(const int & fun_n_oe_start, const int & fun_c, const int & fun_n_fp, const int & fun_n_vel) 
	OL_sagittal.Init_OE(6,4,6,3);
	OL_sagittal.Init_LIPM(zc, StepTime1);

	StepIndex = 0;
	time_cyclic = 0;

	initial_step = true;
	// main loop
	while (realtime < SimTime) {
		//  apply control at lower rate
		RTControl(samplingTime);	

		ddxcom = (xcom - xcop1) / zc*g ; // equation 4.4
		//cout << "outer_loop_theta_stance: " << atan((xcom - xcop1) / hip2ground) << endl;
		xcom += dxcom*dT + 0.5*ddxcom*dT*dT;
		dxcom += ddxcom*dT;
		realtime += dT;
		
	}

	savedata();
	std::cout << "End of simulation" << endl;
	system("pause");
	return 0;
}

void RTControl(double SamplingTime) {
	if ((int)floor(realtime / dT) % (int(SamplingTime / dT)) == 0) { 
		// now update very 20 ms. measure the data out of certain dT
		//---------------------------Add gaussian noise in ycom.-------------------
		std::random_device rd;
		std::mt19937 generator(rd());

		// values near the mean are the most likely
		// standard deviation affects the dispersion of generated values from the mean
		static double noise_std = 0.0002;// If std=0.0001, then 65% of your values will be in the range [-0.0001, 0.0001]. 
		static double mag = 1.0 / 5.0; // if mag =1/4 or 1/5, nearly all of the value would be in [-std,+std]
		std::normal_distribution<> distribution(0.0, (mag*noise_std));

		gaussian_e = distribution(generator);
	
		xcom_e_old = xcom_e;
		xcom_e = xcom + 0*gaussian_e + com_offset;// it is a realistic measurement of CoM with gaussian noise and offset.
		dxcom_e = (xcom_e - xcom_e_old) / SamplingTime;

		xcom_e_f = FilterClass_xcom.applyFilter(xcom_e);
		dxcom_e_f = FilterClass_dxcom.applyFilter(dxcom_e);
		
		//cout << "xcom: " << xcom << "\t dxcom: " << dxcom << endl;
		//cout << "xcom_e_f: " << xcom_e_f << "\t dxcom_e_f: " << dxcom_e_f << endl << endl;
		
		if (realtime <= StartTime) {
			Start1();	
			if (initial_step) {
				traj_angle_attack_init();
				initial_step = false;
				//cout << "xcop1 within start time " << steplength / 2 <<endl<<endl;
			}			
			//double remaintime = StepEndTime - realtime;
			//time_cyclic = StepTime2 - remaintime;
			//xcop1 = steplength/2;
			//cout << "xcop1 within start time " << xcop1<<endl<<endl;
		}
		else if (realtime>StartTime&&realtime <= StopTime) {
			xcop1 = CyclicGaitLateral();			
		}
		else if (realtime>StopTime) {
			double xs = 0;
			double dxs = 0;
			double kp = 50.0;
			double kd = 20.0;
			xcop1 = Stop(zc, kp, kd, xs, dxs, p_plus, p_minus);
		}
		traj_angle_hipknee_generate();
		logdata();
	}
}

double CyclicGaitLateral() {
	double target_vel = 0.5;//OL_sagittal.vel_target(realtime, SimTime, StepTime1); //
	if (realtime>StartTime && store_time.back() < StartTime) {
		// calculate xcapture1 after start period, only once	
		xcapture1 = 0;
		//insipred by equation 4.54, to predict the next lateral foot placement
		xcop1 = xcom + xcapture1;
		footjudgement();
		StepIndex += 1;	
		traj_angle_attack_update();
		//cout << "Stepindex" << StepIndex << " not in dataset.\n" <<"leg_state: "<<leg_state<< endl;
		//cout << "xcom_local: " << (xcom - xcop1) << "\t angle_bottom: " << atan(hip2ground / abs(xcom - xcop1)) << endl << endl;

	}
	if (realtime > StepEndTime) {
		StepIndex += 1;
	
		//-------------------Assign the global foot placement--------------
		//xcop1 = xcop2;
			
			theta_swing0 = atan((xcom-xcop1) / hip2ground);
			l_swing0 = sqrt(pow((xcom - xcop1), 2) + pow(hip2ground, 2));
		// Place the estimation of swing foot location, then this become the support foot location
			xcop1 = xcom + xcapture1;

			//theta_stance0 = atan((xcom - xcop1) / hip2ground);//theta_swingf;
			//l_stance0 = sqrt(pow((xcom - xcop1), 2) + pow(hip2ground, 2));

			// NEW ADD
			OL_sagittal.OnlineEsitmation_main(xcom_e_f, dxcom_e_f, 1, xcapture1,0.5, StepIndex);
			// NEW ADD
			//OL_sagittal.OnlineEsitmation_vel(xcom - xcop1, dxcom, 1, StepIndex);

			// Estimation of final velocity(LIPM)
			// xstate is with respect to the support foot // why time is change but xstate does not changes a lot
			xstate = LIPM.StateEvolution(zc,(xcom - xcop1), dxcom, StepTime2);
			pos_f_predict = xstate[0];
			vel_f_predict = xstate[1];

		// next placement of swing foot(LIPM)
			//swingfoot placement respect to prediction CoM position at the end
			// of current step

			// NEW ADD
			//OL_sagittal.OnlineEsitmation_fp(xcom - xcop1, dxcom, StepIndex, xcapture1,StepIndex);

			xcapture1 = LIPM.SagittallPos(zc,vel_f_predict, target_vel, StepTime2);

			//theta_stancef = atan((pos_f_predict) / hip2ground);
			//l_stancef = sqrt(pow(pos_f_predict, 2) + pow(hip2ground, 2));
			theta_swingf  = atan((-xcapture1) / hip2ground);
			l_swingf = sqrt(pow(hip2ground, 2)+pow(-xcapture1, 2) );

		//-----------------------------Update the time---------------------
		footjudgement();
		StepEndTime += StepTime2;
		StepTime1 = StepTime2; // assign the next step time to current step time


		//-----------------------------Update stance/swing foot status------------
		//traj_angle_attack_update();
		//cout << "Stepindex" << StepIndex << " in dataset.\n" << "leg_state: " << leg_state << endl ;
		//cout << "xcom_local: " << (xcom - xcop1) <<"\t angle_bottom: "<< atan(hip2ground/abs(xcom - xcop1))<<endl << endl;
	}
		 //-------------------check data-------------------
		//cout << "x accel: " << ddxcom*zc / (xcom - xcop1) << endl;
		//cout << "Stepindex" << StepIndex << "  in dataset" << endl;
		//cout << "----xcom: " << xcom << "\t dxcom: " << dxcom << endl;
		//cout << "----xcom_e_f: " << xcom_e_f << "\t dxcom_e_f: " << dxcom_e_f << endl << endl;
	
	//  below  apply FP control at lower rate
	double remaintime = StepEndTime - realtime;
	time_cyclic = StepTime2 - remaintime;

	/*
	xstate = LIPM.StateEvolution(zc, xcom - xcop1, dxcom, remaintime);  // xstate is with respect to the support foot // why time is change but xstate does not changes a lot
	xcapture1 = LIPM.SagittallPos(zc, xstate[1], target_vel, StepTime2);// ideally, vel in and out are same in sagittal, 0.5 is target velocity
	xcop2 = xcop1 + xstate[0] + xcapture1;  // xcop2 is the global position
	//cout << xcop2 << endl;
	*/
	
	/*
	//xcop_peredicted_next_step_e_f = xcop1 + xstate[0] + xcapture1;

	xstate_e_f = LIPM.StateEvolution(zc, xcom_e_f - xcop1, dxcom_e_f, time);  // ystate is with respect to the support foot																										 //ystate_noise_filter = LIPM.StateEvolution(zc, ycom_noise_filter-ycop1, dycom, time);
	xcapture1_e_f_predict = LIPM.SagittallPos(zc, xstate_e_f[1], 0.5, StepTime2);
	xcapture1_e_f_observe= xcop2- xcom_e_f;
	//xcapture2_noise_filter_observe = xcop2-xcom_noise;
	//xcapture2_noise_filter_observe = xstate[0]+xcapture2-xstate_noise_filter[0];
	xcop_peredicted_next_step_e_f = xcop1 + xstate_e_f[0] + xcapture1_e_f_observe;
	//xcop_peredicted_next_step_e_f = xcop1 + xstate_e_f[0] + xcapture2_e_f_predict;
	*/
	return xcop1;
}

void traj_angle_attack_init() {


	//theta_swing0 = 0;//atan((xcom - xcop1) / hip2ground);
	//l_swing0 = sqrt(pow((0), 2) + pow(hip2ground, 2));

	//theta_swingf = atan((-0.4) / hip2ground);
	//l_swingf = sqrt(pow(hip2ground, 2) + pow(-0.4, 2));

	
	/*define gait parameters*/
	phi = deg2rad(25);//inter leg angle
	theta_slope = deg2rad(6);
	lift = 0.15;//0.05;
	Tstep = StepTime2;

	//r_max = sqrt(upperleg*upperleg + lowerleg*lowerleg - 2 * upperleg*lowerleg*cos(pi - 5.0 / 180.0*pi));
	r_max = sqrt(upperleg*upperleg + lowerleg*lowerleg - 2 * upperleg*lowerleg*cos(pi - theta_slope));
	l_swing_min = fullleg - lift;// 5 cm less	
						   //step length, distance between swing foot and support foot, cos law
	steplength = sqrt(2 * r_max*r_max - 2 * r_max*r_max*cos(phi));

	//double angle_bottom = 0.5*(pi - phi);//等腰三角形两边底角
	angle_bottom = 0.5*(pi - phi);//等腰三角形两边底角 ,pi/2;//
	double angle_attack = angle_bottom + theta_slope;//支撑脚和水平面锐夹角
	theta_stance0 = -(0.5*pi - angle_attack);//支撑脚和重力线夹角 theta(0)<0 theta theta_stance0->theta_stancef
											 //theta_stance0 = atan2(cos(phi)-ratio,sin(phi));	

	theta_stancef = theta_stance0 + phi;//支撑脚和重力线夹角 theta(f)
	double angle_swing = 0.5*pi - 0.5*phi - theta_slope;//angle between swing leg and level ground
	theta_swing0 = theta_stancef;//swing leg angle from gravity line(0), start from '+'
	theta_swingf = theta_stance0;//swing leg angle from gravity line(0), end with '-', right hand rule
	r_strike = r_max*sin(angle_swing) / sin(angle_attack);//sine law

	// 0->00, 1->10, 2->01, 3->11; left stance=1; right stance=1
	leg_state = 1;//10

	l_stance0 = r_strike;
	l_stancef = r_max;
	l_swing0  = r_max;
	l_swingf  = r_strike;
	if (theta_stance0 >= 0)
	{
		cout << "attack angle more than 90 degree, reduce virtual slope" << endl;
		Sleep(2000);//system("PAUSE");
		exit(0);
	}

	////display

	//cout << "r_max" << r_max << endl;
	//cout << "l_swing_min" << l_swing_min << endl;
	//cout << "theta_stance0 \t" << theta_stance0*180.0 / pi << "\t degree" << endl;
	//cout << "theta_stancef \t" << theta_stancef*180.0 / pi << "\t degree" << endl;
	//cout << "Tc \t" << Tc << " s" << endl;
	//cout << "Tstep \t" << Tstep << " s" << endl;
	//cout << "r_max-l_swing_min \t" << r_max - l_swing_min << " m" << endl;
	
}

void traj_angle_hipknee_generate() {
	//----------------------------stance foot trajectory----------------------------
	double T_stance0 = 0.5*Tstep;//here to tune the time for leg extension
	double T_stancef = Tstep;//0.95*

	
	//if (time_cyclic <= T_stancef)
	if (time_cyclic <= Tstep)
	{
		theta_stance = atan((xcom - xcop1) / hip2ground);//Poly_stance(0, theta_stance0, T_stancef, theta_stancef, time_cyclic);
		r_stance = sqrt(pow(hip2ground, 2) + pow((xcom - xcop1), 2));//Poly(T_stance0, r_strike, T_stancef, r_max, time_cyclic);
	}
	//cout << "time_cyclic: " << time_cyclic << endl;
	//cout << "xcom local: " << (xcom - xcop1) <<"\t theta_stance: "<< theta_stance << atan(0.5)<<endl<<endl;

	/*
	else
	{
		theta_stance = theta_stancef;
		r_stance = theta_stancef;//hip2ground /cos(theta_stancef);
	}
	*/
	/*

	if (time_cyclic<T_stancef)
	{
		r_stance = r_strike;
	}
	else if ((time_cyclic >= T_stance0) && (time_cyclic<T_stancef))
	{
		r_stance = sqrt(pow(zc,2)+pow((xcom-xcop1),2));//Poly(T_stance0, r_strike, T_stancef, r_max, time_cyclic);
	}
	else
	{
		r_stance = r_max;
	}
	*/

	//----------------------------swing foot trajectory, polynomial method----------------------------
	// polar coordinate method
	double T_swing0 = 0.5*Tstep;//here to tune the time for leg extension
	double T_swingf = Tstep;
	

	if (time_cyclic<T_swingf)// swing angle control
	{	
		// Poly(double t0,double x0,double tf,double xf, double time)
		theta_swing = Poly(theta_swing0, theta_swingf, 0, T_swingf, time_cyclic);
	}
	else if (time_cyclic >= T_swingf)
	{
		theta_swing = theta_swingf;
	}

	if (time_cyclic <= T_swing0)// swing leg length control
	{
		r_swing = Poly( l_swing0, l_swing_min, 0.0, T_swing0, time_cyclic);
	}
	else if ((time_cyclic>T_swing0) && (time_cyclic <= T_swingf))
	{
		r_swing = Poly(l_swing_min,  l_swingf, T_swing0, T_swingf, time_cyclic);
	}
	else if (time_cyclic>T_swingf)
	{
		r_swing = l_swingf;
	}

	if (leg_state == 1)
	{
		theta_left = theta_stance;
		r_left = r_stance;
		theta_right = theta_swing;
		r_right = r_swing;
	}
	else if (leg_state == 3)
	{
		theta_left = theta_swing;
		r_left = r_swing;
		theta_right = theta_stance;
		r_right = r_stance;
	}
	/*
	cout << "time_cyclic: " << time_cyclic << endl;
	cout << "xcom local: " << (xcom - xcop1) << endl;
	cout << "theta_left: " << theta_left << "\t r_left: " << r_left << endl;
	cout << "theta_right: " << theta_right << "\t r_right: " << r_right << endl << endl;
	*/
	//double angleL[3];
	//double angleR[3];
	IK2D(theta_left, r_left, angleL);
	IK2D(theta_right, r_right, angleR);


	/*
	store_traj[0][loop] = angleL[0];
	store_traj[1][loop] = angleL[1];
	//store_traj[2][loop]=angleL[2];
	store_traj[2][loop] = angleR[0];
	store_traj[3][loop] = angleR[1];
	//store_traj[5][loop]=angleR[2];
	*/
}

void footjudgement() {
	if (leg_state == 1)//if previous is left single support
	{
		leg_state = 3;//set now as right single support
	}
	else if (leg_state == 3)//if previous is right single support
	{
		leg_state = 1;//set now as left single support
	}

}

void traj_angle_attack_update(){

}
double Poly_stance(double t0, double x0, double tf, double xf, double time) {
	//time is the current time
	double dx;
	double T;
	double t;
	double s;// output
	dx = xf - x0;
	T = tf - t0;
	t = time - t0;

	if (t >= 0 && t <= T)
	{
		s = x0 + atan((xcom-xcop1)/ hip2ground);
	}
	else if (t<0)
	{
		s = x0;
	}
	else if (t>T)
	{
		s = xf;
	}
	//cout << "theta_stance: " << s << endl;
	return s;
}

void IK2D(double theta, double l0, double *angle)
{
	// the defination of angles are with respect to the local axis based on right handed rule in the conventional x-y-z coordinate
	double lh;
	double ls;
	//double C --- l0 //the distance between hip and ankle joint
	double l0_reset;//reset C if C is out of range    
	double kneeMin;
	double kneeMax;
	double kneeExtentionMax;
	double kneeExtentionMin;
	//double alpha;// vector angle from hip to foot in radian 
	//double beta;//knee angle in radian
	double gama;//angle between upper leg and the vector from hip to foot in radian

	lh = upperleg;//upper leg
	ls = lowerleg;//lower leg	
	//C = l0;

	kneeMin = deg2rad(0);//deg2rad(1)=0.017453292519943
	kneeMax = deg2rad(180);//0.6;// 
	kneeExtentionMax = sqrt(lh*lh + ls*ls - 2 * lh*ls*cos(pi - kneeMin));
	kneeExtentionMin = sqrt(lh*lh + ls*ls - 2 * lh*ls*cos(pi - kneeMax));

	if (l0 >= kneeExtentionMax)
	{
		l0_reset = kneeExtentionMax;
		angle[1] = kneeMin;
	}
	else if (l0 <= kneeExtentionMin)
	{
		cout << "Out of KneeExtentionMin " << endl;
		l0_reset = kneeExtentionMin;
		angle[1] = kneeMax;
	}
	else
	{
		l0_reset = l0;
		angle[1] = pi - acos((lh*lh + ls*ls - l0_reset*l0_reset) / (2.0*lh*ls));
	}
	//gama = asin(ls*sin(pi-angle[1])/C_reset);//sine law
	gama = acos((lh*lh + l0_reset*l0_reset - ls*ls) / (2.0*lh*l0_reset));//cos law
	angle[0] = theta + (-gama);// hip angle in degree
	angle[2] = -(angle[0] + angle[1]);//ankle angle
}

double Poly( double x0, double xf, double t0, double tf, double time)
{
	//time is the current time
	double dx;
	double T;
	double t;
	double s;// output
	dx = xf - x0;
	T = tf - t0;
	t = time - t0;

	if (t >= 0 && t <= T)
	{
		s = x0 + dx*0.5*(1 - cos(pi / T*t));
	}
	else if (t<0)
	{
		s = x0;
	}
	else if (t>T)
	{
		s = xf;
	}
	return s;
}

void logdata() {
	//1-5
	store_time.push_back(realtime);
	store_xcom.push_back(xcom);
	store_dxcom.push_back(dxcom);
	store_ddxcom.push_back(ddxcom);
	store_xcop.push_back(xcop1);
	//6-10
	store_xcom_e_f.push_back(xcom_e_f);
	store_dxcom_e_f.push_back(dxcom_e_f);
	store_xcop_e_f.push_back(xcop_e_f);
	store_xcapture1_e_f_predict.push_back(xcapture1_e_f_predict);
	store_xcapture1_e_f_observe.push_back(xcapture1_e_f_observe);
	//10
	/*store_fp_wrt_com.push_back(OL_sagittal.return_fp_wrt_com);
	store_com_wrt_stance_foot.push_back(xcom - xcop1);
	store_coeff_v0.push_back(OL_sagittal.return_coeff_vel(0, 0));
	store_coeff_vd.push_back(OL_sagittal.return_coeff_vel(1, 0));
	store_coeff_com_offset.push_back(OL_sagittal.return_coeff_vel(2, 0));
*/
	store_fp_wrt_com.push_back(0);
	store_com_wrt_stance_foot.push_back(0);
	store_coeff_v0.push_back(0);
	store_coeff_vd.push_back(0);
	store_coeff_com_offset.push_back(0);

	// traj file
	store_traj_LHip.push_back(angleL[0]);
	store_traj_LKnee.push_back(angleL[1]);
	store_traj_RHip.push_back(angleR[0]);
	store_traj_RKnee.push_back(angleR[1]);
	
	store_theta_stance.push_back(theta_stance);
	store_l_stance.push_back(r_stance);
	store_theta_swing.push_back(theta_swing);
	store_l_swing.push_back(r_swing);

	store_theta_left.push_back(theta_left);
	store_r_left.push_back(r_left);
	store_theta_right.push_back(theta_right);
	store_r_right.push_back(r_right);

	//temporary store 
	store_xcom_e.push_back(xcom_e);
	store_dxcom_e.push_back(dxcom_e);
}

void savedata() {
	ofstream file;
	file.open("../../data/LateralGait.txt");

	for (unsigned int i = 0; i<store_time.size(); i++) {
		// Format for output data
		file << store_time[i] << "\t";
		file << store_xcom[i] << "\t";
		file << store_dxcom[i] << "\t";
		file << store_ddxcom[i] << "\t";
		file << store_xcop[i] << "\t";
		file << store_xcom_e_f[i] << "\t";
		file << store_dxcom_e_f[i] << "\t";
		file << store_xcop_e_f[i] << "\t";
		file << store_xcapture1_e_f_predict[i]<< "\t";
		file << store_xcapture1_e_f_observe[i] << "\t";
		//11
		file << store_fp_wrt_com[i] << "\t";
		file << store_com_wrt_stance_foot[i] << "\t";
		file << store_coeff_v0[i] << "\t";
		file << store_coeff_vd[i] << "\t";
		file << store_coeff_com_offset[i] << "\t";
		file << store_xcom_e[i] << "\t";
		file << store_dxcom_e[i] << "\n";
	};
	file.close();

	file.open("../../data/Traj.txt");
	for (unsigned int i = 0; i<store_time.size(); i++) {
		// Format for output data
		file << store_traj_LHip[i] << "\t";
		file << store_traj_LKnee[i] << "\t";
		file << store_traj_RHip[i] << "\t";
		file << store_traj_RKnee[i] << "\t";
		// temporary logout for check data
		file << store_theta_stance[i] << "\t";
		file << store_l_stance[i]	<< "\t";
		file << store_theta_swing[i] << "\t";
		file << store_l_swing[i]	<< "\t";

		file << store_theta_left[i] << "\t";
		file << store_r_left[i]		<< "\t";
		file << store_theta_right[i] << "\t";
		file << store_r_right[i]	<< "\n";
		//file << store_time[i] << "\n";
	};
	std::cout << "Data saved:)" << endl;
}


double Start1()
{
	
	if (realtime>tswitch && store_time.back() < tswitch)
	{
		double t2;
		//the period that the CoM takes a trip from (xcom-p_plus,dxcom) to (yf-p_plus,dxf(esitmated by  VelocityAtPositionCV))
		t2 = LIPM.RemainingTimePositionSymmetry(zc, xcom - p_plus, dxcom, yf - p_plus);
		cout << t2 << "\t" << StartTime << "\t" << tswitch + t2 << endl;
		StartTime = tswitch + t2;
		StepEndTime = StartTime + StepTime2; // end time of 1st periodic step
	}
	/*
	if (time_cyclic <= Tstep)
	{
		xcom = Poly(0, -0.2, 0, Tstep, time_cyclic);
	}
	*/
	//xcop1 = 0;
	//footjudgement();
	//ddxcom = (xcom-xcop1)/zc*g;
	return xcop1;
}

double Stop(double zc, double kp, double kd, double xs, double dxs, double p_plus, double p_minus)
{
	// note that stop should be enabled in the next coming event when COM cross zero
	// so i should insert another piece of code to judge
	double ddx, maxAcc, minAcc;
	ddx = kp*(xs - xcom) + kd*(dxs - dxcom);//desired position -- xs and velocity -- dxs, are 0

	maxAcc = g*(xcom - p_minus) / zc;//equation 4.4
	minAcc = g*(xcom - p_plus) / zc;

	if (ddx>maxAcc) ddx = maxAcc;
	else if (ddx<minAcc) ddx = minAcc;

	double xcop = xcom - ddx / g*zc;

	return xcop;
}

double Sign(double x)
{
	if (x >= 0) return 1.0;
	else return -1.0;
}

