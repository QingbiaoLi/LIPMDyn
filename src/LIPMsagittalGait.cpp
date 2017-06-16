/* Programmed by: Zhibin Li*/
// this program tests the LIPM dynamcis functions
// 24 April 2014
// updated: 29 May 2014
#include "windows.h"
#include <fstream>
#include <iostream>
#include "LIPMDynClass.h"
#include "FpOnlineEstimationClass.h"
#include <chrono>
#include <random>


#define DEGTORAD(x)  x*M_PI/180.0
#define RADTODEG(x)  x*180.0/M_PI
void logdata();
void savedata();
void RTControl(double SamplingTime);
double Start1(); // output cop
double CyclicGaitLateral(); // output cop
double Stop(double zc, double kp, double kd, double xs, double dxs, double p_plus, double p_minus); // output cop
double Sign(double x);
double LateralDifferentialDrive(double dym, double delta_dy, double ysway);

LIPMDynClass LIPM;
FpOnlineEstimation OL_sagittal;

vector<double> xstate(2, 0);
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
double xcop1, xcop2, ycapture1, ycapture2, ysway;
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
vector<double> store_dis_ddy(1, 0);
vector<double> store_xcop(1, 0);
vector<double> store_PredictTime(1, 0);
vector<double> store_PredictFP(1, 0);
vector<double> store_PredictVel(1, 0);
vector<double> store_PredictPos(1, 0);
double dis_ddy = 0.0;
double dis_time = 0;
double dis_duration = 0.1;
double dis_A;
bool dis_enable = 0;

//additional code for saggittal, online estimation
double gaussian_e;
int StepIndex;
double samplingTime = 0.004;
int data_set_number = 6;
double com_offset = 1.0;

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

	tswitch = LIPM.LateralSwitchTime(zc, xcom, dxcom, yf, p_plus, p_minus, StartTime); // return time of switching phase
	StepEndTime = StartTime + StepTime1; // End time of 1st periodic step.


	OL_sagittal.Init(data_set_number);
	StepIndex = 0;

	while (realtime < SimTime) {
		RTControl(samplingTime);	//  apply control at lower rate

		ddxcom = (xcom - xcop1) / zc*g + dis_ddy; // equation 4.4

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
	//dym lateral velocity it should be while cross the middle line 
	double dym = LIPM.LateralNominalVelocity(zc, p_plus, StepTime1);	//p_plus lateral foot placement (step width),
	if (realtime>StartTime && store_time.back() < StartTime) {// calculate ycapture2 after start period, only once	
		ycapture2 = 0;
		//LIPM.LateralGaitVelocityControl(zc, dym, dxcom, StepTime1); //insipred by equation 4.54, to predict the next lateral foot placement
		xcop1 = xcom + ycapture2;
		//ysway = xcom - xcop1;
		StepIndex += 1;
		cout << realtime << endl;
	}
	if (realtime>StepEndTime)	// here the condition is time based coz the time constraint comes from sagittal plane
	{
		StepIndex += 1;
		//cout << "---StepIndex ---" <<StepIndex<< endl << "xcom: " << xcom-xcop1 << "\t dxcom: " << dxcom << "\t xcop1:" << xcop1<< "\t  ycapture1:" << ycapture1 <<endl<<endl;
		//cout << "---StepIndex ---" << StepIndex << endl << realtime << endl;
		//cout << "---StepIndex ---" << StepIndex << endl << "xcom: " << xcom - xcop1 << endl;
		OL_sagittal.collect_walking_state(xcom - xcop1, dxcom, com_offset, ycapture1, StepIndex);

		OL_sagittal.StateEsimation(dxcom, 0.5, StepIndex);

		cout << "dataset_past_walking_state" << endl << OL_sagittal.dataset_past_walking_state << endl << OL_sagittal.dataset_next_footplacement << endl << "control coeff:" << OL_sagittal.control_coefficient << endl << "footplacement estimation: " << OL_sagittal.footplacement_predict << endl << endl;
		//xcop2 = OL_sagittal.footplacement_predict;

		//double xcop2_tmp = xcop1 + xcom + OL_sagittal.footplacement_predict;

		StepEndTime += StepTime2;
		StepTime1 = StepTime2; // assign the next step time to current step time
		xcop1 = xcop2;
		//cout << "---StepIndex ---" << StepIndex << endl << "xcom: " << xcom - xcop1 << endl;
	}
	//  below  apply FP control at lower rate

	double time = StepEndTime - realtime;

	/*
	double remainTime;
	double dir_drive;
	dir_drive = 0.0*sin(0.5*(realtime - StartTime));//function for walking pattern
	double dym_mod = LateralDifferentialDrive(dym, dir_drive, ysway);
	remainTime = LIPM.RemainingTimeVelocitySymmetry(zc, xcom - xcop1, dxcom, Sign(ysway)*dym_mod);//how much time left to reach a target COM position (Sign(ysway)*dym_mod)
	PredictTime = realtime + remainTime;
	*/
	xstate = LIPM.StateEvolution(zc, xcom - xcop1, dxcom, time);  // xstate is with respect to the support foot // why time is change but xstate does not changes a lot

	ycapture1 = LIPM.SagittallPos(zc, xstate[1], 0.5, StepTime2);// ideally, vel in and out are same in sagittal, 0.5 is target velocity
																 //cout <<"xstate" << xstate[1] << "\t  ycapture1:" << ycapture1 << endl << endl;
																 //ycapture2 = LIPM.LateralGaitVelocityControl(zc, dym_mod, xstate[1], StepTime2);



																 /*
																 if (StepIndex > (data_set_number + 5)) {
																 ycapture1 = LIPM.SagittallPos(zc, xstate[1], 0.7, StepTime2);// ideally, vel in and out are same in sagittal, 0.5 is target velocity
																 }
																 */
	xcop2 = xcop1 + xstate[0] + ycapture1;  // xcop2 is the global position

											/*
											if (StepIndex > (data_set_number + 5)) {
											OL_sagittal.StateEsimation(dxcom, 0.5, StepIndex);//include collect current state;
											xcop2 = xcop1 + xstate[0] + OL_sagittal.footplacement_predict;
											}
											*/
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
	store_PredictTime.push_back(PredictTime);
	store_PredictFP.push_back(xcop2);
	store_PredictVel.push_back(xstate[1]);
	store_dis_ddy.push_back(dis_ddy);
	store_PredictPos.push_back(xstate[0]);
	//11
	store_fp_wrt_com.push_back(ycapture1);
	store_com_wrt_stance_foot.push_back(xcom - xcop1);
	store_coeff_v0.push_back(OL_sagittal.control_coefficient(0, 0));
	store_coeff_vd.push_back(OL_sagittal.control_coefficient(1, 0));
	store_coeff_com_offset.push_back(OL_sagittal.control_coefficient(2, 0));
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
		file << store_PredictTime[i] << "\t";
		file << store_PredictFP[i] << "\t";
		file << store_PredictVel[i] << "\t";
		file << store_dis_ddy[i] << "\t";
		file << store_PredictPos[i] << "\t";
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

	//if (delta_dy>0) // drift to left
	//{
	//	//if(Sign(ysway)<0)// now diverge to right, so it adds abs(delta_dy) to the left in next step
	//	//{
	//	//	dym_mod=dym+abs(delta_dy);
	//	//}
	//	//else if(Sign(ysway)>0) // diverge to right
	//	//{
	//	//	dym_mod=dym-abs(delta_dy);
	//	//}
	//	dym_mod=dym-Sign(ysway)*delta_dy;
	//}
	//else // drift to right
	//{
	//	//if(Sign(ysway)<0)// now diverge to right, so next step diverge to left, it deducts abs(delta_dy)
	//	//{
	//	//	dym_mod=dym-abs(delta_dy);
	//	//}
	//	//else if(Sign(ysway)>0) // now diverge to left, next step diverge to right, adds abs(delta_dy)
	//	//{
	//	//	dym_mod=dym+abs(delta_dy);
	//	//}
	//	dym_mod=dym+Sign(ysway)*delta_dy;
	//}
	dym_mod = dym - Sign(ysway)*delta_dy; // to summarize, only one equation
	return dym_mod;
}