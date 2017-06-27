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
#include "FpOnlineEstimationClass.h"
#include "FilterClass.h"


#define DEGTORAD(x)  x*M_PI/180.0
#define RADTODEG(x)  x*180.0/M_PI
void logdata();
void savedata();
void RTControl(double SamplingTime);
double Start1(); // output cop
double CyclicGaitLateral(); // output cop
double Stop(double zc, double kp, double kd, double xs, double dxs, double p_plus, double p_minus); // output cop
double Sign(double x);

LIPMDynClass LIPM;
FpOnlineEstimation OL_sagittal;

vector<double> xstate(2, 0);
vector<double> xstate_e_f(2, 0);
double StepTime1 = 0.5;
double StepTime2 = 0.5;
double StartTime;
double SimTime = 10.0;
double StopTime = SimTime - 0.0;

double dT = 0.0002;
double zc = 0.5;
double g = 9.81;
double realtime = 0.0;
double xcom, dxcom, ddxcom;
double xcop1, xcop2, xcapture1, xcapture2, ysway;//xcop1,xcop2 is the foot placement is 
double StepEndTime;
double PredictTime; // predicted landing time

double yf = 0.0;
double p_plus = 0.10;// 0.10;0.0;// probably the lateral step size  
double p_minus = -0.10; //-0.10;
double tswitch;

vector<double> store_time(1, 0);
vector<double> store_xcom(1, 0);
vector<double> store_dxcom(1, 0);
vector<double> store_ddxcom(1, 0);
vector<double> store_xcop(1, 0);
vector<double> store_xcom_e_f(1, 0);
vector<double> store_dxcom_e_f(1, 0);
vector<double> store_xcop_e_f(1, 0);
vector<double> store_xcapture1_e_f_predict(1, 0);
vector<double> store_xcapture1_e_f_observe(1, 0);

double dis_time = 0;
double dis_duration = 0.1;
double dis_A;
bool dis_enable = 0;

//additional code for saggittal, online estimation

int StepIndex;
double samplingTime = 0.004;
int data_set_number = 6;
double com_offset = 0.0;
double fp_wrt_com;

//add noise
double gaussian_e;
//static double mag;
//static double xcom_e_old, xcom_e, dxcom_e;
double xcom_e_old, xcom_e, dxcom_e, xcom_e_f_old;
double xcop_peredicted_next_step_e_f;
double ddxcom_e;
double xcapture2_e_predict;
double xcapture2_e_observe;
double xcom_e_old_filter, xcom_e_f, dxcom_e_f;
double xcop_e_pure_f;
double xcop_e_f;
double ddxcom_e_f;
double xcapture1_e_f_predict;
double xcapture1_e_f_observe;

FilterClass FilterClass_xcom, FilterClass_dxcom, FilterClass_xcop, FilterClass_xcapture2, FilterClass_xstate;

vector<double> store_fp_wrt_com(1, 0);
vector<double> store_com_wrt_stance_foot(1, 0);
vector<double> store_coeff_v0(1, 0);
vector<double> store_coeff_vd(1, 0);
vector<double> store_coeff_com_offset(1, 0);

int main() {
	StartTime = 0.8; // Duration of inital step strategy. (time of initiation of lateral gait)
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

	//FilterClass_xcop.butterworth(samplingTime, 10.0, 1);
	//FilterClass_xcapture2.butterworth(samplingTime, 5.0, 1);
	//FilterClass_xstate.butterworth(samplingTime, 10.0, 1);

	tswitch = LIPM.LateralSwitchTime(zc, xcom, dxcom, yf, p_plus, p_minus, StartTime); // return time of switching phase
	StepEndTime = StartTime + StepTime1; // End time of 1st periodic step.

	//online estimation setting
	OL_sagittal.Init(data_set_number);
	StepIndex = 0;

	while (realtime < SimTime) {
		RTControl(samplingTime);	//  apply control at lower rate

		ddxcom = (xcom - xcop1) / zc*g ; // equation 4.4
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
	if ((int)floor(realtime / dT) % (int(SamplingTime / dT)) == 0) { // now update very 20 ms. measure the data out of certain dT
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
		xcom_e = xcom + gaussian_e + com_offset;// it is a realistic measurement of CoM with gaussian noise and offset.
		dxcom_e = (xcom_e - xcom_e_old) / SamplingTime;

		//ddxcom_e = (xcom_e-xcop_peredicted_next_step_e)/zc*g;
		//xcom_e_old_filter = xcom_e_f;
		//dxcom_e_f=(xcom_e_f-ycom_e_old_filter)/SamplingTime;
		//dxcom_e_f=FilterClass_dxcom.applyFilter(dxcom_e);
		//dxcom_e_f=FilterClass_dxcom.applyFilter(dxcom);

		xcom_e_f = FilterClass_xcom.applyFilter(xcom_e);
		dxcom_e_f = FilterClass_dxcom.applyFilter(dxcom_e);
		//cout << "xcom: " << xcom << "\t dxcom: " << dxcom << endl;
		//cout << "xcom_e_f: " << xcom_e_f << "\t dxcom_e_f: " << dxcom_e_f << endl << endl;
		if (realtime <= StartTime) {
			xcop1 = 0;//Start1();
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
		logdata();
	}
}

double CyclicGaitLateral() {

	if (realtime>StartTime && store_time.back() < StartTime) {// calculate xcapture2 after start period, only once	
		xcapture2 = 0;
		//insipred by equation 4.54, to predict the next lateral foot placement
		//xcop1 = xcom + xcapture2;
		xcop1 = xcom + xcapture2;
		//ysway = xcom - xcop1;
		StepIndex += 1;
		cout <<"Stepindex"<<StepIndex<< " not in dataset"<< endl;
	}
	if (realtime > StepEndTime) {// at the beginning of new step
		StepIndex += 1;
		
		//-------------------Ideal case---------------
		//LIPM model based foot placement prediction
		double vel_tar = 0.5;//OL_sagittal.vel_target(realtime, SimTime, StepTime1);
		xcapture1 = LIPM.SagittallPos(zc, dxcom, vel_tar, StepTime2);
		
		
		//Online Estimation
		OL_sagittal.collect_walking_state(xcom - xcop1, dxcom, 1, xcapture1, StepIndex);

		OL_sagittal.StateEsimation(dxcom, 0.5, StepIndex);
		
		//-------------------with gaussian noise case---------------
		/*
		double vel_tar = OL_sagittal.vel_target(realtime, SimTime, StepTime1);
		xcapture1 = LIPM.SagittallPos(zc, dxcom_e_f, vel_tar, StepTime2);
		xcop1 = xcop_peredicted_next_step_e_f;//xcom_e_f + xcapture1;
		*/

		//-------------------Assign the global foot placement--------------
		//xcop1 = xcom + OL_sagittal.footplacement_predict;

		//-------------------check data-------------------
		if (StepIndex >= (data_set_number + 2 + 2)) {
			static const auto runOnce = [] { cout << "-------------Start online estimation----------------" << endl; return true; }();
			OL_sagittal.collect_current_walking_state(dxcom, vel_tar);
			OL_sagittal.calculate_model_coeff(OL_sagittal.dataset_past_walking_vel, OL_sagittal.dataset_next_footplacement);
			OL_sagittal.estimate_walking_state_next_step(OL_sagittal.dataset_current_walking_state, OL_sagittal.model_coeff);
			cout << "dataset_past_walking_state_stack: " << endl << OL_sagittal.dataset_past_walking_state_stack << endl << "dataset_past_walking_state: " << endl << OL_sagittal.dataset_past_walking_vel << endl;
			cout << OL_sagittal.dataset_next_footplacement << endl << "control coeff:" << OL_sagittal.model_coeff << endl << "current walking state: " << OL_sagittal.dataset_current_walking_state << endl;
			cout << "footplacement estimation: " << OL_sagittal.footplacement_predict << endl << endl;
		}
		xcop1 = xcom + OL_sagittal.footplacement_predict;
		
		//-------------------update the time-------------------
		StepEndTime += StepTime2;
		StepTime1 = StepTime2;
	}	// here the condition is time based coz the time constraint comes from sagittal plane

	return xcop1;
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
	store_fp_wrt_com.push_back(OL_sagittal.footplacement_predict);
	store_com_wrt_stance_foot.push_back(xcom - xcop1);
	store_coeff_v0.push_back(OL_sagittal.model_coeff(0, 0));
	store_coeff_vd.push_back(OL_sagittal.model_coeff(1, 0));
	store_coeff_com_offset.push_back(OL_sagittal.model_coeff(2, 0));
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
		file << store_coeff_com_offset[i] << "\n";
	};
	file.close();
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

	if (realtime < tswitch) // now update very 20 ms.
	{
		xcop1 = p_minus; //otherwise remain the CoP at current position.
	}
	else
	{
		xcop1 = p_plus;
	}
	//ddxcom = (xcom-xcop1)/zc*g;
	return xcop1;
}

double Stop(double zc, double kp, double kd, double xs, double dxs, double p_plus, double p_minus)
{
	// note that stop should be enabled in the next coming event when COM cross zero
	// so i should insert another piece of code to judge
	double ddx, maxAcc, minAcc;
	ddx = kp*(xs - xcom) + kd*(dxs - dxcom);//desired position -- ys and velocity -- dys, are 0

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



// 加到未来程序里
double LateralDifferentialDrive(double dym, double delta_dy, double ysway) {
	double dym_mod; // modified desired dym, to have differential drive

	if (abs(delta_dy) >= dym) { // delta_dy cannot be bigger than dym
		delta_dy = Sign(delta_dy)*dym;
	}

	dym_mod = dym - Sign(ysway)*delta_dy; // to summarize, only one equation
	return dym_mod;
}