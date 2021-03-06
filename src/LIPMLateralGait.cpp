/* Programmed by: Zhibin Li*/
// this program tests the LIPM dynamcis functions
// 24 April 2014
// updated: 29 May 2014
#include "windows.h"
#include <fstream>
#include <iostream>
#include "LIPMDynClass.h"

#include <chrono>
#include <random>


#define DEGTORAD(x)  x*M_PI/180.0
#define RADTODEG(x)  x*180.0/M_PI
void logdata();
void savedata();
void RTControl(double SamplingTime);
double Start1(); // output cop
double CyclicGaitLateral(); // output cop
double Stop(double zc, double kp, double kd, double ys, double dys, double p_plus, double p_minus); // output cop
double Sign(double x);
double LateralDifferentialDrive(double dym, double delta_dy, double ysway);

LIPMDynClass LIPM;

vector<double> ystate(2,0);
double StepTime1=0.5;
double StepTime2=0.5;
double StartTime;
double SimTime=10.0;
double StopTime = SimTime-0.0;

double dT=0.000001;
double zc=0.5; 
double g=9.81;
double realtime=0.0;
double ycom, dycom, ddycom;
double ycop1, ycop2, ycapture1, ycapture2, ysway;
double StepEndTime;
double PredictTime; // predicted landing time

double yf=0.0;
double p_plus=0.10;// probably the lateral step size 
double p_minus=-0.10;
double tswitch;

vector<double> store_time(1,0);
vector<double> store_ycom(1,0);
vector<double> store_dycom(1,0);
vector<double> store_ddycom(1,0);
vector<double> store_dis_ddy(1,0);
vector<double> store_ycop(1,0);
vector<double> store_PredictTime(1,0);
vector<double> store_PredictFP(1,0);
vector<double> store_PredictVel(1,0);
vector<double> store_PredictPos(1,0);
double dis_ddy=0.0;
double dis_time=0;
double dis_duration=0.1;
double dis_A;
bool dis_enable=0;

double gaussian_e;

int main() {
	StartTime = 0.8; // Duration of inital step strategy. (time of initiation of lateral gait)
	ycom = 0.0; // Initial CoM position in global frame.
	dycom = 0.0; // Initial CoM velocity in global frame.
	ddycom = 0.0; // Initial CoM acceleration in global frame.	
	ycop2 = 0; // Inital predicted foot placement in global frame.

	tswitch = LIPM.LateralSwitchTime(zc, ycom, dycom, yf, p_plus, p_minus, StartTime); // return time of switching phase
	StepEndTime = StartTime+StepTime1; // End time of 1st periodic step.

	while(realtime < SimTime) {
		RTControl(0.005);	//  apply control at lower rate

		ddycom = (ycom-ycop1)/zc*g + dis_ddy; // equation 4.4

		ycom += dycom*dT + 0.5*ddycom*dT*dT;

		dycom += ddycom*dT;
		realtime += dT;
	}

	savedata();
	std::cout << "End of simulation" << endl;
	system("pause");	
    return 0;
}

void logdata() {
	store_time.push_back(realtime);
	store_ycom.push_back(ycom);
	store_dycom.push_back(dycom);
	store_ddycom.push_back(ddycom);
	store_ycop.push_back(ycop1);

	store_PredictTime.push_back(PredictTime);
	store_PredictFP.push_back(ycop2);
	store_PredictVel.push_back(ystate[1]);
	store_dis_ddy.push_back(dis_ddy);
	store_PredictPos.push_back(ystate[0]);
}

void savedata() {
	ofstream file;
	file.open("../../data/LateralGait.txt");

	for (unsigned int i=0;i<store_time.size();i++) {
		// Format for output data
		file<<store_time[i]<<"\t";
		file<<store_ycom[i]<<"\t";
		file<<store_dycom[i]<<"\t";
		file<<store_ddycom[i]<<"\t";
		file<<store_ycop[i]<<"\t";
		file<<store_PredictTime[i]<<"\t";		
		file<<store_PredictFP[i]<<"\t";
		file<<store_PredictVel[i]<<"\t";
		file<<store_dis_ddy[i]<<"\t";
		file<<store_PredictPos[i]<<"\n";
	};
	file.close();
	std::cout << "Data saved:)" <<endl;
}

void RTControl(double SamplingTime) {
	if( (int)floor(realtime/dT)%(int(SamplingTime/dT)) ==0 ) { // now update very 20 ms. measure the data out of certain dT
		if (realtime<=StartTime) {
			ycop1 = Start1();
		}
		else if ( realtime>StartTime&&realtime<=StopTime ) {
			 ycop1 = CyclicGaitLateral();
		}
		else if (realtime>StopTime) {
			double ys=0;
			double dys=0;
			double kp=50.0;
			double kd=20.0;
			ycop1=Stop(zc, kp, kd, ys, dys, p_plus, p_minus);
		}
		logdata();
	}
}

double Start1()
{
	if ( realtime>tswitch && store_time.back() < tswitch )
	{
		double t2;
		//the period that the CoM takes a trip from (ycom-p_plus,dycom) to (yf-p_plus,dxf(esitmated by  VelocityAtPositionCV))
		t2=LIPM.RemainingTimePositionSymmetry(zc, ycom-p_plus, dycom, yf-p_plus);
		cout<<t2<<"\t"<<StartTime<<"\t"<<tswitch+t2<<endl;
		StartTime=tswitch+t2;	
		StepEndTime=StartTime+StepTime2; // end time of 1st periodic step
	}

	if( realtime < tswitch ) // now update very 20 ms.
	{
		ycop1=p_minus; //otherwise remain the CoP at current position.
	}
	else
	{
		ycop1=p_plus;
	}
	//ddycom = (ycom-ycop1)/zc*g;
	return ycop1;
}

double CyclicGaitLateral() {
	//dym lateral velocity it should be while cross the middle line 
	double dym = LIPM.LateralNominalVelocity(zc, p_plus, StepTime1);	//p_plus lateral foot placement (step width),
	if (realtime>StartTime && store_time.back() < StartTime ) {
		ycapture2 = LIPM.LateralGaitVelocityControl(zc, dym, dycom, StepTime1); //insipred by equation 4.54, to predict the next lateral foot placement
		ycop1 = ycom+ycapture2;
		ysway = ycom-ycop1;
	}
	if (realtime>StepEndTime)	// here the condition is time based coz the time constraint comes from sagittal plane
	{			
		StepEndTime += StepTime2;
		StepTime1 = StepTime2; // assign the next step time to current step time
		ycop1 =  ycop2;
		ysway = ycom-ycop1;
	}
	//  below  apply FP control at lower rate

	double time = StepEndTime - realtime;
	double remainTime;

	double dir_drive;
	dir_drive = 0.0*sin(0.5*(realtime-StartTime));//function for walking pattern  

	double dym_mod=LateralDifferentialDrive(dym, dir_drive, ysway);

	remainTime = LIPM.RemainingTimeVelocitySymmetry(zc, ycom-ycop1, dycom, Sign(ysway)*dym_mod);//how much time left to reach a target COM position (Sign(ysway)*dym_mod)
	PredictTime = realtime+remainTime;
 	ystate = LIPM.StateEvolution(zc, ycom-ycop1, dycom, time);  // ystate is with respect to the support foot // why time is change but ystate does not changes a lot
	ycapture1 = LIPM.LateralGaitSymmetryTimeControl(zc, ystate[1], StepTime2); // FP position with respect to COM

	ycapture2 = LIPM.LateralGaitVelocityControl(zc, dym_mod, ystate[1], StepTime2);
	ycop2 = ycop1 + ystate[0] + ycapture2;  // ycop2 is the global position

	return ycop1;
}

double Stop(double zc, double kp, double kd, double ys, double dys, double p_plus, double p_minus)
{
	// note that stop should be enabled in the next coming event when COM cross zero
	// so i should insert another piece of code to judge
	double ddy, maxAcc, minAcc;
	ddy = kp*(ys-ycom)+kd*(dys-dycom);//desired position -- ys and velocity -- dys, are 0

	maxAcc = g*(ycom-p_minus)/zc;//equation 4.4
	minAcc = g*(ycom-p_plus)/zc;

	if (ddy>maxAcc) ddy=maxAcc;
	else if (ddy<minAcc) ddy=minAcc;

	double ycop=ycom-ddy/g*zc;

	return ycop;
}

double Sign(double x)
{
	if (x >= 0) return 1.0;
	else return -1.0;
}

// 加到未来程序里
double LateralDifferentialDrive(double dym, double delta_dy, double ysway) {
	double dym_mod; // modified desired dym, to have differential drive
	
	if (abs(delta_dy)>=dym) { // delta_dy cannot be bigger than dym
		delta_dy=Sign(delta_dy)*dym;
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
	dym_mod=dym-Sign(ysway)*delta_dy; // to summarize, only one equation
	return dym_mod;
}