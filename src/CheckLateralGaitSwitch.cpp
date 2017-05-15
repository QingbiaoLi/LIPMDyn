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
void Sim1();
void Sim2();

vector<double> store_time;
vector<double> store_ycom;
vector<double> store_dycom;
vector<double> store_ddycom;
vector<double> store_ycop;
vector<double> store_tswitch;
vector<double> store_PredictVel;

LIPMDynClass LIPM;

vector<double> ystate(2,0);
double StepTime1=0.5;
double StepTime2=0.5;
double SimTime=2;
double dT=0.00002;
double zc=0.5; 
double g=9.81;
double realtime=0;
double ycom, dycom, ddycom;
double ycop1, ycop2, ycapture, ysway;
double StepEndTime;
double PredictTime; // predicted landing time
double PredictVel;

double yy0=0.0;
double dy0=0.0;
double yf=-0.05;
double p_plus=0.10;
double p_minus=-0.10;
double tswitch;
double y0switch;
double dy0switch;
int sampleSet=9;
int sampleNum=0;
double tscan=0.2;
double tend;

/*conclusion: there is not big difference in the switching time if the total time is very small, 
0.5s, or larger (>1.5s, which is very similar to infinite solution already), the swithch time is 
always between 0.1s to 0.15s. so actually we can consider this number as a general/universal
easy to use number to initiate the lateral gait, and leave the fine adjustment of sway time by 
the support foot COP control in reality. So humans may learn this know-how already, so when we
initiate lateral gait, we apply a EMG impulse that generates muscle power that lasts for about 
0.15s, and it works in most cases already! */

bool method=0;

void main()
{
	while(sampleNum<sampleSet)
	{
		tend=sampleNum*tscan+0.4;
		SimTime = tend;
		realtime=0;
		yy0=0;
		dy0=-0.15;
		ycom=yy0;
		dycom=dycom;
		ddycom=0;
		tswitch = LIPM.LateralSwitchTime(zc, yy0, dy0, yf, p_plus, p_minus, tend);
		PredictVel = LIPM.LateralSwitchEndVelocity(zc, yy0, dy0, yf, p_plus, p_minus, tswitch);
		store_tswitch.push_back(tend);
		store_tswitch.push_back(tswitch);
		

		while(realtime<=SimTime)
		{	
			if (method)
			{
				Sim1();	
			}else
			{
				Sim2();
			}
			logdata();
			realtime += dT;
		}
		store_PredictVel.push_back(PredictVel);
		store_PredictVel.push_back(dycom);

		sampleNum++;
	}
	savedata();
	cout<< "End of simulation" <<endl;
	//system("pause");	
}

void logdata()
{
	store_time.push_back(realtime);
	store_ycom.push_back(ycom);	
	store_dycom.push_back(dycom);
	store_ddycom.push_back(ddycom);
	store_ycop.push_back(ycop1);

	//store_PredictTime.push_back(PredictTime);
	//store_PredictFP.push_back(ycop2);
}


void savedata()
{
	ofstream file;
	file.open("../Data/LSData1.txt");
	for (int i=0;i<store_time.size();i++)
	{
		file<<store_time[i]<<"\t";
		file<<store_ycom[i]<<"\t";
		file<<store_dycom[i]<<"\t";
		file<<store_ddycom[i]<<"\t";
		file<<store_ycop[i]<<"\n";		
	};
	file.close();

	file.open("../Data/LSSwitch1.txt");
	for (int i=0;i<store_tswitch.size();i++)
	{
		file<<store_tswitch[i]<<"\n";		
	};
	file.close();	

	file.open("../Data/LSVel1.txt");
	for (int i=0;i<store_PredictVel.size();i++)
	{
		file<<store_PredictVel[i]<<"\n";		
	};
	file.close();

	
	cout<< "Data saved:)" <<endl;
}

void Sim1()
{
		if ( realtime>tswitch && store_time.back() < tswitch )
		{
			double t2;
			t2=LIPM.RemainingTime(zc, ycom-p_plus, dycom, yf-p_plus);
			cout<<t2<<endl;
			//SimTime=tswitch+t2;			
		}

		if( realtime < tswitch ) // now update very 20 ms.
		{
			ycop1=p_minus;
		}
		else
		{
			ycop1=p_plus;
		}
		ddycom = (ycom-ycop1)/zc*g;
		ycom += dycom*dT + 0.5*ddycom*dT*dT;
		dycom += ddycom*dT;
}


void Sim2()
{
	if ( realtime>tswitch && store_time.back() < tswitch )
	{
		y0switch=ycom;
		dy0switch=dycom;
		double t2;
		t2=LIPM.RemainingTime(zc, ycom-p_plus, dycom, yf-p_plus);				
		cout<<t2<<"\t"<<SimTime<<"\t"<<tswitch+t2<<endl;
		SimTime=tswitch+t2;
	}

	if( realtime < tswitch ) // now update very 20 ms.
	{
		ycop1=p_minus;
		ystate = LIPM.StateEvolution(zc, yy0-p_minus, dy0, realtime);
		ycom=ystate[0]+p_minus;
		dycom=ystate[1];
	}
	else
	{
		ycop1=p_plus;
		ystate = LIPM.StateEvolution(zc, y0switch-p_plus, dy0switch, realtime-tswitch);
		ycom=ystate[0]+p_plus;
		dycom=ystate[1];
	}		
}