/* Programmed by: Zhibin Li*/
// this program tests the LIPM dynamcis functions
// 24 April 2014
#include "windows.h"
#include <fstream>
#include <iostream>
#include "LIPMDynClass.h"

#define DEGTORAD(x)  x*M_PI/180.0
#define RADTODEG(x)  x*180.0/M_PI
void logdata();
void savedata();

LIPMDynClass LIPM;

vector<double> ystate(2,0);
double StepTime1=0.5;
double StepTime2=0.5;
double SimTime=3.0;
double dT=0.00002;
double zc=0.5; 
double g=9.81;
double realtime=0;
double ycom, dycom, ddycom;
double ycop1, ycop2, ycapture, ysway;
double StepEndTime;
double PredictTime; // predicted landing time

vector<double> store_time;
vector<double> store_ycom;
vector<double> store_dycom;
vector<double> store_ycop;
vector<double> store_PredictFP;
vector<double> store_PredictTime;

void main()
{
	ycom=0.0;
	dycom=0.0;
	ddycom=0;	// initial condition
	
	ycop1 = LIPM.LateralGaitSymmetryTimeControl(zc, dycom, StepTime1);
	ycop2=0;
	ycapture=0;
	ysway = ycom-ycop1;
	StepEndTime=StepTime1; // end time of current step


	while(realtime<SimTime)
	{	

		if (realtime>StepEndTime)	// here the condition is time based coz the time constraint comes from sagittal plane
		{			
			StepEndTime += StepTime2;
			StepTime1 = StepTime2; // assign the next step time to current step time
			ycop1 =  ycop2;
			ysway = ycom-ycop1;
		}

		/////  below  apply FP control at lower rate
		if( (int)floor(realtime/dT)%(int(0.005/dT)) ==0 ) // now update very 20 ms.
		{
			double time = StepEndTime - realtime;
			PredictTime = realtime+LIPM.RemainingTime(zc, ycom-ycop1, dycom, ysway); // use ysway here to keep the gait symetry
			ystate = LIPM.StateEvolution(zc, ycom-ycop1, dycom, time);  // ystate is with respect to the support foot
			ycapture = LIPM.LateralGaitSymmetryTimeControl(zc, ystate[1], StepTime2);
			ycop2 = ycop1 + ystate[0]+ycapture;  // ycop2 is the global position
			//cout<< ycop2 <<endl;
			logdata();
			//system("pause");
		}

		/////  above apply FP control at lower rate

		ddycom = (ycom-ycop1)/zc*g;
		ycom+= dycom*dT + 0.5*ddycom*dT*dT;
		dycom+= ddycom*dT;
		realtime += dT;
		//Sleep(0.5);
	}

	savedata();
	cout<< "End of simulation" <<endl;
	//system("pause");	
}

void logdata()
{
	store_time.push_back(realtime);
	store_ycom.push_back(ycom);
	store_ycop.push_back(ycop1);

	store_PredictTime.push_back(PredictTime);
	store_PredictFP.push_back(ycop2);
}


void savedata()
{
	ofstream file;
	file.open("../Data/FootPlacement.txt");
	for (int i=0;i<store_time.size();i++)
	{
		file<<store_time[i]<<"\t";
		file<<store_ycom[i]<<"\t";
		file<<store_ycop[i]<<"\t";
		file<<store_PredictTime[i]<<"\t";
		file<<store_PredictFP[i]<<"\n";
	};
	file.close();
	cout<< "Data saved:)" <<endl;
}


