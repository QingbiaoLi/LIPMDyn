/* Programmed by: Zhibin Li*/
// this program tests the LIPM dynamcis functions
// 24 April 2014
// updated: 29 May 2014
// updated: Dec 2016 by  Qingbiao Li for Online Estimation
#include "windows.h"
#include <fstream>
#include <iostream>
#include "LIPMDynClass.h"

#include "FilterClass.h"
#include <chrono>
#include <random>


#define DEGTORAD(x)  x*M_PI/180.0
#define RADTODEG(x)  x*180.0/M_PI
void logdata();
void savedata();
void RTControl(double SamplingTime);
double Start1(); // output cop

				 //double CyclicGaitLateral(); // output cop
void CyclicGaitLateral(); // output cop

double Stop(double zc, double kp, double kd, double ys, double dys, double p_plus, double p_minus); // output cop
double Sign(double x);
double LateralDifferentialDrive(double dym, double delta_dy, double ysway);

LIPMDynClass LIPM;

vector<double> ystate(2, 0);
vector<double> ystate_noise(2, 0);
vector<double> ystate_noise_filter(2, 0);
double StepTime1 = 0.5;  // Step peroid
double StepTime2 = 0.5;	// ??? 查看一下steptime1和steptime2的区别
double StartTime = 1.0; // time of initiation of lateral gait
double SimTime = 8.0; // simulation period 
double StopTime = SimTime - 0.0;// when the simulation stops

								//double dT=0.000001;
double dT = 0.0002;
double zc = 0.5; //need to add error for the height in main()
double g = 9.81;
double realtime = 0.0;
double ycom, dycom, ddycom;
double ycop1, ycop2, ycapture1, ycapture2, ysway;
//double ycop2, ycapture1, ycapture2, ysway;
double StepEndTime;
double PredictTime; // predicted landing time

double yf = 0.0;
double p_plus = 0.10;
double p_minus = -0.10;
double tswitch;

vector<double> store_time(1, 0);
vector<double> store_ycom(1, 0);
vector<double> store_dycom(1, 0);
vector<double> store_ddycom(1, 0);
vector<double> store_dis_ddy(1, 0);
vector<double> store_ycop(1, 0);
vector<double> store_PredictTime(1, 0);
vector<double> store_PredictFP(1, 0);
vector<double> store_PredictVel(1, 0);

double dis_ddy = 0.0;
double dis_time = 0;
double dis_duration = 0.1;
double dis_A;
bool dis_enable = 0;

//additional definition
double gaussian_e;
//static double mag;
//static double ycom_noise_old, ycom_noise, dycom_noise;
double ycom_noise_old, ycom_noise, dycom_noise, ycom_noise_filter_old;
double ycop_peredicted_next_step_noise;
double ddycom_noise;
double ycapture2_noise_predict;
double ycapture2_noise_observe;
vector<double> store_ycom_noise(1, 0);
vector<double> store_dycom_noise(1, 0);
vector<double> store_ycop_noise(1, 0);
vector<double> store_ycapture2_noise_predict(1, 0);
vector<double> store_ycapture2_noise_observe(1, 0);

FilterClass FilterClass_ycom, FilterClass_dycom, FilterClass_ycop, FilterClass_ycapture2, FilterClass_ystate;
//static double ycom_noise_old_filter, ycom_noise_filter, dycom_noise_filter;
double ycom_noise_old_filter, ycom_noise_filter, dycom_noise_filter;
double ycop_peredicted_next_step_noise_pure_filter;
double ycop_peredicted_next_step_noise_filter;
double ddycom_noise_filter;
double ycapture2_noise_filter_predict;
double ycapture2_noise_filter_observe;
vector<double> store_ycom_noise_filter(1, 0);
vector<double> store_dycom_noise_filter(1, 0);
vector<double> store_ycop_noise_filter(1, 0);
//vector<double> store_ycapture2_noise_filter(1,0);
vector<double> store_ycapture2_noise_filter_predict(1, 0);
vector<double> store_ycapture2_noise_filter_observe(1, 0);
vector<double> store_ycop_peredicted_next_step_noise_pure_filter(1, 0);


double step_number = 1;
double step_end = 0;
double step_end_index = 0;
double step_time;
double step_remainTime;
void logdata_ls();
void savedata_ls();
vector<double>store_time_nstep(1, 0);
vector<double>store_CoM_offset_nstep(1, 0);
vector<double>store_dycom_nstep(1, 0);
vector<double>store_ycom_nstep(1, 0);
vector<double>store_ycapture2_nstep_predict(1, 0);


vector<double>store_dycom_filter_nstep(1, 0);
vector<double>store_ycom_filter_nstep(1, 0);
vector<double>store_ycapture2_filter_nstep_predict(1, 0);
vector<double>store_ycapture2_filter_nstep_observe(1, 0);


vector<double>store_Predictpos(1, 0);
vector<double>store_ycapture2(1, 0);
vector<double>store_ystate_noise_filter(2, 0);

vector<double>store_gaussian(1, 0);
vector<double>store_gaussian_filter(1, 0);

vector<double>store_predict_pos(1, 0);
vector<double>store_predict_vel(1, 0);
vector<double>store_predict_local_footplacment(1, 0);

vector<double>store_predict_pos_noise(1, 0);
vector<double>store_predict_vel_noise(1, 0);
vector<double>store_predict_local_footplacment_noise(1, 0);

vector<double>store_predict_pos_filter(1, 0);
vector<double>store_predict_vel_filter(1, 0);
vector<double>store_predict_local_footplacment_filter(1, 0);
vector<double>store_step_end(1, 0);
vector<double>store_step_time(1, 0);
vector<double>store_step_remainTime(1, 0);

double samplingTime = 0.004;
double zc_real = zc + 0.0; //Zc error 
double CoM_offset = 0.00;//0.01 // CoM offset setting

int main()
{
	StartTime = 0.8;
	ycom = 0.0;
	dycom = 0.0;
	ddycom = 0.0;	// initial y COM condition	

	ycom_noise = 0.0;
	dycom_noise = 0.0;
	ddycom_noise = 0.0;

	ycom_noise_filter = 0.0;
	dycom_noise_filter = 0.0;
	ddycom_noise_filter = 0.0;

	ycop2 = 0;
	//ycapture=0;
	//ysway = ycom-ycop1;
	tswitch = LIPM.LateralSwitchTime(zc, ycom, dycom, yf, p_plus, p_minus, StartTime);
	double EndVel = LIPM.LateralSwitchEndVelocity(zc, ycom, dycom, yf, p_plus, p_minus, tswitch); // export for n+1 velocity????
	StepEndTime = StartTime + StepTime1; // end time of 1st periodic step


										 //setting for filter
	FilterClass_ycom.butterworth(samplingTime, 10.0, 1);
	FilterClass_dycom.butterworth(samplingTime, 10.0, 1);


	/*FilterClass_ycom.moving_average_filter(8);
	FilterClass_dycom.moving_average_filter(8);
	*/
	FilterClass_ycop.butterworth(samplingTime, 10.0, 1);
	FilterClass_ycapture2.butterworth(samplingTime, 5.0, 1);
	FilterClass_ystate.butterworth(samplingTime, 10.0, 1);


	while (realtime<SimTime)
	{
		//  apply control at lower rate  frequency of noise 200 Hz.
		RTControl(samplingTime);
		//if ( (realtime>StartTime+2.2)&&(realtime<SimTime-5.0) )
		//{		
		//	if (!dis_enable && (realtime-dis_time>0.4))
		//	{
		//		//dis_A=2.0*(double)(rand()-0.5*RAND_MAX)/(double)RAND_MAX;	//RAND_MAX 32767
		//		dis_A=Sign(rand()-0.5*RAND_MAX);
		//		dis_enable=1;
		//		dis_time=realtime;
		//	}
		//	if (dis_enable)
		//	{
		//		dis_ddy=2.0*dis_A*sin(2*M_PI*(realtime-dis_time)/(2*dis_duration));
		//		if (realtime>(dis_time+dis_duration))
		//		{
		//			dis_enable=0;
		//			dis_time=realtime;
		//		}
		//	}
		//}

		//---------------------------------------------------------------------???????
		ddycom = (ycom - ycop1) / zc_real*g + dis_ddy; // equation 4.4 zc_real is kind of way to cheat a real Zc,dis_ddy = 0;

													   // Original code
		ycom += dycom*dT + 0.5*ddycom*dT*dT;
		dycom += ddycom*dT;

		//dycom_noise+= ddycom_noise*dT;
		realtime += dT;
		//Sleep(0.5);

		//dycom_noise_filter=FilterClass_dycom.applyFilter(dycom_noise);

	}

	savedata();
	savedata_ls();
	cout << "End of simulation" << endl;
	//Sleep(1000);
	system("pause");
	return 1;
}

void RTControl(double SamplingTime)
{

	if ((int)floor(realtime / dT) % (int(SamplingTime / dT)) == 0) // now update very 20 ms. // samplingTime=0.002, dT=0.0002
	{

		//---------------------------Add gaussian noise in ycom.-------------------
		std::random_device rd;
		std::mt19937 generator(rd());

		// values near the mean are the most likely
		// standard deviation affects the dispersion of generated values from the mean
		static double noise_std = 0.0002;// If std=0.0001, then 65% of your values will be in the range [-0.0001, 0.0001]. 
		static double mag = 1.0 / 5.0; // if mag =1/4 or 1/5, nearly all of the value would be in [-std,+std]
		std::normal_distribution<> distribution(0.0, (mag*noise_std));

		gaussian_e = distribution(generator);


		ycom_noise_old = ycom_noise;
		ycom_noise = ycom + gaussian_e + CoM_offset;// it is a realistic measurement of CoM with gaussian noise and offset.
		dycom_noise = (ycom_noise - ycom_noise_old) / SamplingTime;

		//ddycom_noise = (ycom_noise-ycop_peredicted_next_step_noise)/zc*g + dis_ddy;
		//ycom_noise_old_filter = ycom_noise_filter;
		//dycom_noise_filter=(ycom_noise_filter-ycom_noise_old_filter)/SamplingTime;
		//dycom_noise_filter=FilterClass_dycom.applyFilter(dycom_noise);
		//dycom_noise_filter=FilterClass_dycom.applyFilter(dycom);


		ycom_noise_filter = FilterClass_ycom.applyFilter(ycom_noise);
		dycom_noise_filter = FilterClass_dycom.applyFilter(dycom_noise);
		step_end_index = step_end_index + 1;
		/*	dycom_noise_filter= (ycom_noise_filter-ycom_noise_filter_old)/SamplingTime;
		ycom_noise_filter_old=ycom_noise_filter;*/
		//---------------------------------------------------------------------

		if (realtime <= StartTime)
		{
			ycop1 = Start1();
		}
		else if (realtime>StartTime&&realtime <= StopTime)
		{

			//vector<double> ycop_temp(2,0);
			CyclicGaitLateral();
			//ycop1 = ycop_temp[0];
			//ycop_peredicted_next_step_noise=ycop_temp[1];
			//cout<< ycop1 << "\t" << ycop_peredicted_next_step_noise << endl;
			//ycop1 = CyclicGaitLateral_e();
		}
		else if (realtime>StopTime)
		{
			double ys = 0;
			double dys = 0;
			double kp = 50.0;
			double kd = 20.0;
			ycop1 = Stop(zc, kp, kd, ys, dys, p_plus, p_minus);
		}

		logdata();

	}
	//else {
	//ycom_noise = FilterClass_ycom.applyFilter(ycom);
	//dycom_noise = FilterClass_dycom.applyFilter(dycom);

	//}
	/*ycom_noise_filter=FilterClass_ycom.applyFilter(ycom_noise);
	dycom_noise_filter=FilterClass_dycom.applyFilter(dycom_noise);*/
}



double Start1()
{
	if (realtime>tswitch && store_time.back() < tswitch)
	{
		double t2;
		t2 = LIPM.RemainingTimePositionSymmetry(zc, ycom - p_plus, dycom, yf - p_plus);
		cout << t2 << "\t" << StartTime << "\t" << tswitch + t2 << endl;
		StartTime = tswitch + t2;
		StepEndTime = StartTime + StepTime2; // end time of 1st periodic step
	}

	if (realtime < tswitch) // now update very 20 ms.
	{
		ycop1 = p_minus;
		ycop2 = -ycop1;
		ycop_peredicted_next_step_noise_filter = -ycop1;

	}
	else
	{
		ycop1 = p_plus;
		ycop2 = -ycop1;
		ycop_peredicted_next_step_noise_filter = -ycop1;


	}
	//ddycom = (ycom-ycop1)/zc*g;
	return ycop1;
}

void CyclicGaitLateral()
{
	double dym = LIPM.LateralNominalVelocity(zc, p_plus, StepTime1);//lateral velocity it should be while cross the middle line 
																	//LateralNominalSway(double zc, double stepwidth, double StepTime)
																	//double vel = LIPM.LateralNominalVelocity(zc, 0.2, 0.1);
																	//double x = LIPM.RemainingTime(zc, -p_plus, vel, -p_plus); // use ysway here to keep the gait symetry

	if (realtime>StartTime && store_time.back() < StartTime)
	{
		//double remainTime;
		//ycapture1 = LIPM.LateralGaitSymmetryTimeControl(zc, dycom, StepTime1);
		//ycop1 = ycom+ycapture1;
		//ysway = ycom-ycop1;
		//remainTime = LIPM.RemainingTime(zc, ycom-ycop1, dycom, ysway);
		ycapture2 = LIPM.LateralGaitVelocityControl(zc, dym, dycom, StepTime1);
		ycop1 = ycom_noise_filter + ycapture2;
		ysway = ycom_noise_filter - ycop1;
		//remainTime = LIPM.RemainingTime(zc, ycom-ycop1, dycom, ysway);	
		//cout<<endl;
		//// add noise
		//ycapture2_noise = LIPM.LateralGaitVelocityControl(zc, dym, dycom_noise, StepTime1);
		//ycop_peredicted_next_step_noise = ycom_noise+ycapture2_noise;

		//ycop[0]=ycop1;
		//ycop[1]=ycop_peredicted_next_step_noise;
		//cout<< ycop[0] << "\t" << ycop[1] << endl;


	}
	if (realtime>StepEndTime)	// here the condition is time based coz the time constraint comes from sagittal plane
	{
		StepEndTime += StepTime2;
		StepTime1 = StepTime2; // assign the next step time to current step time
		logdata_ls();
		//ycop1 =  ycop2;
		ycop1 = ycop_peredicted_next_step_noise_filter;
		ysway = ycom - ycop1;
		//cout<<"Sign ysway:"<<Sign(ysway)<<endl;
		//Sleep(100);
		//double verr = dycom-Sign(dycom)*dym;
		//cout<<"verr:"<<verr<<endl;

		//log the velocity--dycom, center of mass position -- ycom, foot placement--ycaputure of last loop
		step_number += 1;
		step_end = step_end_index;
	}

	/////  below  apply FP control at lower rate

	double time = StepEndTime - realtime;
	double remainTime;
	step_time = time;

	double dir_drive;
	dir_drive = 0.0*sin(0.5*(realtime - StartTime)); // sinware for lateral position
													 //dir_drive=-0.05;

	double dym_mod = LateralDifferentialDrive(dym, dir_drive, ysway);
	//double dym_mod =LateralDifferentialD`rive(dym_noise, dir_drive, ysway);

	//cout<<"dym_mod:"<<dym_mod<<"\t"<<Sign(ysway)<<endl; //Sleep(10);

	//remainTime = LIPM.RemainingTimePositionSymmetry(zc, ycom-ycop1, dycom, ysway); // use ysway here to keep the gait symetry
	remainTime = LIPM.RemainingTimeVelocitySymmetry(zc, ycom - ycop1, dycom, Sign(ysway)*dym_mod);    // update remain time according to velocity symmetry
	step_remainTime = remainTime;
	//StepEndTime=realtime+remainTime;
	PredictTime = realtime + remainTime;
	ystate = LIPM.StateEvolution(zc, ycom - ycop1, dycom, time);  // ystate is with respect to the support foot,[0] for position,[1] for veolicty  where ' ycom-ycop1 = x - p* 'is to transfer global coordinate into local coordinatesssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
	cout << "ystate[0]" << ystate[0] << endl;
	ycapture1 = LIPM.LateralGaitSymmetryTimeControl(zc, ystate[1], StepTime2); // FP position with respect to COM(local coordinate system??)

	ycapture2 = LIPM.LateralGaitVelocityControl(zc, dym_mod, ystate[1], StepTime2);//lateral foot placement

	double c = 0.0;
	//ycop2 = ycop1 + ystate[0]+1.0122*ycapture1+0*(1.0-c)*ycapture2;  // ycop2 is the global position
	ycop2 = ycop1 + ystate[0] + c*ycapture1 + (1.0 - c)*ycapture2;  // ycop2 is the global position of the predicted foot placement for the next step


//----------------------------------------------------------------------------------------------------------------------

//-----------------------------get filtered signal -----------------------------

//Use ideal equation --logically more reasonable

	ystate_noise_filter = LIPM.StateEvolution(zc, ycom_noise_filter - ycop1, dycom_noise_filter, time);  // ystate is with respect to the support foot
																										 //ystate_noise_filter = LIPM.StateEvolution(zc, ycom_noise_filter-ycop1, dycom, time); 
	ycapture2_noise_filter_predict = LIPM.LateralGaitVelocityControl(zc, dym_mod, ystate_noise_filter[1], StepTime2);
	//ycapture2_noise_filter_observe = ycop2-ycom_noise_filter;
	//ycapture2_noise_filter_observe = ycop2-ycom_noise;
	//ycapture2_noise_filter_observe = ystate[0]+ycapture2-ystate_noise_filter[0];
	ycop_peredicted_next_step_noise_filter = ycop1 + ystate_noise_filter[0] + ycapture2_noise_filter_predict;

	//Sleep(100);
	//return ycop1,ycop_peredicted_next_step_noise;
}


double Stop(double zc, double kp, double kd, double ys, double dys, double p_plus, double p_minus)
{
	// note that stop should be enabled in the next coming event when COM cross zero
	// so i should insert another piece of code to judge
	double ddy, maxAcc, minAcc;
	ddy = kp*(ys - ycom) + kd*(dys - dycom);

	maxAcc = g*(ycom - p_minus) / zc;
	minAcc = g*(ycom - p_plus) / zc;

	if (ddy>maxAcc) ddy = maxAcc;
	else if (ddy<minAcc) ddy = minAcc;

	double ycop = ycom - ddy / g*zc;

	return ycop;
}


double Sign(double x)
{
	if (x >= 0) return 1.0;
	else return -1.0;
}

// 加到未来程序里
double LateralDifferentialDrive(double dym, double delta_dy, double ysway)
{
	double dym_mod; // modified desired dym, to have differential drive

	if (abs(delta_dy) >= dym) // delta_dy cannot be bigger than dym
	{
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


void logdata()
{
	store_time.push_back(realtime);
	store_ycom.push_back(ycom);
	store_dycom.push_back(dycom);
	store_ddycom.push_back(ddycom);
	store_ycop.push_back(ycop1);// global current foot placement in ideal simulation
								//store_ycop.push_back(ycop[0]);

	store_PredictTime.push_back(PredictTime);
	store_PredictFP.push_back(ycop2);// global predicted foot placement in ideal simulation
	store_PredictVel.push_back(ystate[1]);
	store_dis_ddy.push_back(dis_ddy);

	//additional data
	//--------------------Log noise data---------------------
	store_ycom_noise.push_back(ycom_noise);
	store_dycom_noise.push_back(dycom_noise);
	store_ycop_noise.push_back(ycop_peredicted_next_step_noise); // global foot placement in simulation with noise
	store_ycapture2_noise_predict.push_back(ycapture2_noise_predict);
	store_ycapture2_noise_observe.push_back(ycapture2_noise_observe);

	//--------------------Log filtered data---------------------
	store_ycom_noise_filter.push_back(ycom_noise_filter);
	//store_dycom_noise_filter.push_back(ystate_noise_filter[0]);//dycom_noise_filter
	store_dycom_noise_filter.push_back(dycom_noise_filter);
	store_ycop_noise_filter.push_back(ycop_peredicted_next_step_noise_filter);// global foot placement in simulation with noise after filter
																			  // local foot placement in simulation with noise after filter
	store_ycapture2_noise_filter_predict.push_back(ycapture2_noise_filter_predict);
	store_ycapture2_noise_filter_observe.push_back(ycapture2_noise_filter_observe);


	//log the gassuian noise data
	store_gaussian.push_back(gaussian_e);
	//add filter for the gassuian noise
	store_gaussian_filter.push_back(FilterClass_ycapture2.applyFilter(gaussian_e));
	store_ycop_peredicted_next_step_noise_pure_filter.push_back(ycop_peredicted_next_step_noise_pure_filter);

	store_predict_pos.push_back(ystate[0]);
	store_predict_vel.push_back(ystate[1]);
	store_predict_local_footplacment.push_back(ycapture2);

	store_predict_pos_noise.push_back(ystate_noise[0]);
	store_predict_vel_noise.push_back(ystate_noise[1]);
	store_predict_local_footplacment_noise.push_back(ycapture2_noise_predict);

	store_predict_pos_filter.push_back(ystate_noise_filter[0]);
	store_predict_vel_filter.push_back(ystate_noise_filter[1]);
	store_predict_local_footplacment_filter.push_back(ycapture2_noise_filter_predict);

	store_step_time.push_back(step_time);
	store_step_remainTime.push_back(step_remainTime);


}


void savedata()
{
	ofstream file;
	//file.open("../Data/LateralGait.txt");
	file.open("../Data/LateralGait_with_Gerror_change_loop.txt");
	for (unsigned int i = 0; i<store_time.size(); i++)
	{

		// Format for output data
		file << store_time[i] << "\t";
		file << store_ycom[i] << "\t";
		file << store_dycom[i] << "\t";
		file << store_ddycom[i] << "\t";
		file << store_ycop[i] << "\t";
		file << store_PredictTime[i] << "\t";
		file << store_PredictFP[i] << "\t";
		file << store_PredictVel[i] << "\t";
		file << store_dis_ddy[i] << "\t";

		//additional data
		//--------------------Log noise data---------------------
		file << store_ycom_noise[i] << "\t";
		file << store_dycom_noise[i] << "\t";
		file << store_ycop_noise[i] << "\t";
		file << store_ycapture2_noise_predict[i] << "\t";
		file << store_ycapture2_noise_observe[i] << "\t";


		//--------------------Log filtered data---------------------
		file << store_ycom_noise_filter[i] << "\t";
		file << store_dycom_noise_filter[i] << "\t";
		file << store_ycop_noise_filter[i] << "\t";
		file << store_ycapture2_noise_filter_predict[i] << "\t";
		file << store_ycapture2_noise_filter_observe[i] << "\t";

		//log the gassuian noise data 
		file << store_gaussian[i] << "\t";
		file << store_gaussian_filter[i] << "\t";
		file << store_ycop_peredicted_next_step_noise_pure_filter[i] << "\t";


		file << store_predict_pos[i] << "\t";
		file << store_predict_vel[i] << "\t";
		file << store_predict_local_footplacment[i] << "\t";

		file << store_predict_pos_noise[i] << "\t";
		file << store_predict_vel_noise[i] << "\t";
		file << store_predict_local_footplacment_noise[i] << "\t";

		file << store_predict_pos_filter[i] << "\t";
		file << store_predict_vel_filter[i] << "\t";
		file << store_predict_local_footplacment_filter[i] << "\t";
		file << store_step_time[i] << "\t";
		file << store_step_remainTime[i] << "\n";

	};
	file.close();
	cout << "Data saved:)" << endl;
}

//--------------------------- Export data for Online Esimtaiton -------------------------------------------
void logdata_ls()
{
	store_time_nstep.push_back(realtime);
	//store_CoM_offset_nstep.push_back(CoM_offset);
	store_dycom_filter_nstep.push_back(dycom_noise_filter);// velocity of Ceter of mass with noise after fileter
	store_ycom_filter_nstep.push_back(ycom_noise_filter - ycop1);// position of Ceter of mass to local coordinate system with noise after fileter
	store_ycapture2_filter_nstep_predict.push_back(ycapture2_noise_filter_predict);// local foot placement
	store_ycapture2_filter_nstep_observe.push_back(ycapture2_noise_filter_observe);

	store_ycom_nstep.push_back(ycom - ycop1);
	store_dycom_nstep.push_back(dycom);
	store_ycapture2_nstep_predict.push_back(ycapture2);
	store_step_end.push_back(step_end);

	//store_ycapture2_nstep.push_back(ycop_peredicted_next_step_noise_filter);
}


void savedata_ls()
{
	ofstream file;
	//file.open("../Data/LateralGait.txt");
	file.open("../Data/OnlineEstimation.txt");
	for (unsigned int j = 0; j<step_number; j++)
	{

		// Format for output data
		file << store_time_nstep[j] << "\t";
		//file<<store_CoM_offset_nstep[j]<<"\t";
		file << store_ycom_filter_nstep[j] << "\t";
		file << store_dycom_filter_nstep[j] << "\t";
		file << store_ycapture2_filter_nstep_predict[j] << "\t";
		file << store_ycapture2_filter_nstep_observe[j] << "\t";

		file << store_ycom_nstep[j] << "\t";
		file << store_dycom_nstep[j] << "\t";
		file << store_ycapture2_nstep_predict[j] << "\t";
		file << store_step_end[j] << "\n";
	};
	file.close();
	cout << "Online Estimation Data saved:)" << endl;
}
//---------------------------------------------------------------------------------------
