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
#include "RobotParameter.h"


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
		xcom_e = xcom + gaussian_e + com_offset;// it is a realistic measurement of CoM with gaussian noise and offset.
		dxcom_e = (xcom_e - xcom_e_old) / SamplingTime;

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
		xcop1 = xcom + xcapture2;
		StepIndex += 1;
		cout <<"Stepindex"<<StepIndex<< " not in dataset"<< endl;
	}
	if (realtime > StepEndTime) {
		StepIndex += 1;
		// origin
		//xcapture1 = LIPM.SagittallPos(zc, dxcom, 0.5, StepTime2);
		//OL_sagittal.collect_walking_state(xcom - xcop1, dxcom, 0, xcapture1, StepIndex);
		//OL_sagittal.StateEsimation(dxcom, 0.5, StepIndex);
		//xcapture1 = OL_sagittal.footplacement_predict;
		//xcop1 = xcop_peredicted_next_step_e_f;// xcom + xcapture1;
		// add noise and filter	
		//OL_sagittal.StateEsimation(dxcom_e_f, 0.5, StepIndex);
		//xcapture1 = LIPM.SagittallPos(zc, dxcom_e_f, 0.5, StepTime2);	
		//xcop1 = xcom_e_f + xcapture1;
		//cout << "xcom: "<<xcom-xcop1<<"\t xcom_e_f: "<<xcom_e_f - xcop1 <<"\t xcop2: " << xcop2 <<"\t xcop_peredicted_next_step_e_f: "<< xcop_peredicted_next_step_e_f <<endl<<endl;
		//cout << "xcapture1: " << xcapture1 << "\t xcapture1_e_f_observe: " << xcapture1_e_f_observe << endl << endl;
		//-------------------online estimation--------------------
		/*
		//cout << "---xcop_peredicted_next_step_e_f: " << xcop_peredicted_next_step_e_f << endl;
		//OL_sagittal.collect_walking_state(xcom_e_f - xcop1, dxcom_e_f, 1, xcapture2_e_f_predict, StepIndex);
		cout << "---StepIndex ---" << StepIndex << endl << "xcom: " << xcom_e_f << "\t " << xcom - xcop1 << endl;
		//OL_sagittal.StateEsimation(dxcom_e_f, 0.5, StepIndex);
		if (StepIndex >= (data_set_number + 2 + 1)) {
			//xcop_peredicted_next_step_e_f = xcom_e_f +OL_sagittal.footplacement_predict;
			//xcop_peredicted_next_step_e_f = xcom + OL_sagittal.footplacement_predict;
			cout << "dataset_past_walking_state_stack" << OL_sagittal.dataset_past_walking_state_stack << "dataset_past_walking_state" << endl << OL_sagittal.dataset_past_walking_state << endl << OL_sagittal.dataset_next_footplacement << endl << "control coeff:" << OL_sagittal.control_coefficient << endl << "current walking state" << OL_sagittal.dataset_current_walking_state << endl << "footplacement estimation: " << OL_sagittal.footplacement_predict << endl << endl;
			//xcop1 = xcop_peredicted_next_step_e_f; // xcop2;//
			//cout << "fp_wrt_com: " << OL_sagittal.dataset_next_footplacement << endl;
		}
		*/

		//-------------------Assign the global foot placement--------------
		xcop1 = xcop2;
		StepEndTime += StepTime2;
		StepTime1 = StepTime2; // assign the next step time to current step time

		 //-------------------check data-------------------
		//cout << "x accel: " << ddxcom*zc / (xcom - xcop1) << endl;
		//cout << "Stepindex" << StepIndex << "  in dataset" << endl;
		//cout << "----xcom: " << xcom << "\t dxcom: " << dxcom << endl;
		//cout << "----xcom_e_f: " << xcom_e_f << "\t dxcom_e_f: " << dxcom_e_f << endl << endl;
	}

	//  below  apply FP control at lower rate
	double time = StepEndTime - realtime;
	double vel_tar = 0.5;// OL_sagittal.vel_target(realtime, SimTime, StepTime1);

	xstate = LIPM.StateEvolution(zc, xcom - xcop1, dxcom, time);  // xstate is with respect to the support foot // why time is change but xstate does not changes a lot
	xcapture1 = LIPM.SagittallPos(zc, xstate[1], vel_tar, StepTime2);// ideally, vel in and out are same in sagittal, 0.5 is target velocity
	xcop2 = xcop1 + xstate[0] + xcapture1;  // xcop2 is the global position

	//cout << xcop2 << endl;

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

	// traj file
	store_traj_LHip.push_back(1);
	store_traj_LKnee.push_back(2);
	store_traj_RHip.push_back(3);
	store_traj_RKnee.push_back(4);
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

	file.open("../../data/tarj.txt");
	for (unsigned int i = 0; i<store_time.size(); i++) {
		// Format for output data
		file << store_traj_LHip[i] << "\t";
		file << store_traj_LKnee[i] << "\t";
		file << store_traj_RHip[i] << "\t";
		file << store_traj_RKnee[i] << "\n";
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

