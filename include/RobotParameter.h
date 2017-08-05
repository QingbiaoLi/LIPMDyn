#ifndef _robotPara_H// header guards
#define _robotPara_H
#include <cmath>


//#include <math.h>       /* atan */
const double pi = 3.14159265358979323846264;
/*pls be awared that this self defined funtion deg2rad() doesnt support
any arithmetic calculation(+_x/), no matter inside or outside bracket*/
#define deg2rad(x) x*pi/180;
#define rad2deg(x) x*180/pi;
const double g = 9.81;
/*
struct robotlimb{
double max;// max length
double min;// min length
double length;
bool contact;// stance leg/contact arm or swing leg/idle arm
//bool aerial;// swing leg or idle arm
};
robotlimb legR, legL, armR, armL;
*/

/* robot dimension , mass, inertia */
const double upperleg =  0.342;//0.5;//
const double lowerleg =  0.355; //0.55;//
const double fullleg = upperleg + lowerleg;// full leg length
const double hip_offset = 0.07260;// from pelvis center to hip joint


double r_max ;//maximum leg extension, also it is leg length of new swing leg in when heelstrike

double l_swing_min;
double r_strike;// min stance leg length at heel-strike: r_strike->r_max
double r_stance; // 
double r_swing;// min swing leg length: r_max->r_min->r_strike
short int leg_state = 3;// 0->00, 1->10, 2->01, 3->11; left stance=1; right stance=1
double rp = 0.30; // pendulum length from COM to the support center(FT sensor center), seems not used

double theta_swing;//swing leg length in polar coordinate
double theta_swing0;//swing leg length in polar coordinate
double theta_swingf;//swing leg length in polar coordinate

double theta_stance;//stance leg length in polar coordinate
double theta_stance0;//stance leg length in polar coordinate
double theta_stancef;//stance leg length in polar coordinate

double l_stance0;
double l_stancef;

double l_swing0;
double l_swingf;


double r_left;//left leg length in polar coordinate
double theta_left;//left leg angle in polar coordinate
double r_right;//right leg length in polar coordinate
double theta_right;//right leg angle in polar coordinate

/*mass, inertia*/
const double M_tot = 17.0;
double I_tot = M_tot*rp*rp;// approximated inertia in single mass model
						   //const double  k_approx=0.98;

double Wc = sqrt(g / r_max);// natural frequency
double Tc = 1.0 / Wc;

/*gait parameters*/
double phi;//inter leg angle
double theta_slope;


double steplength; // step length
double lift;// 5 cm ground clearance
double ratio;// r_strike/r_max
double Tstep;// period of one step
double Tstep_des;// period of one step
int numberofsteps = 0;

bool initial_step;
double angleL[3];
double angleR[3];
//r_min=fullleg-lift;// 5 cm less

// LIPMsagittalGait origin header

void logdata();
void savedata();
void RTControl(double SamplingTime);
double Start1(); // output cop
double CyclicGaitLateral(); // output cop
double Stop(double zc, double kp, double kd, double xs, double dxs, double p_plus, double p_minus); // output cop
double Sign(double x);

LIPMDynClass LIPM;
OnlineEstimation OL_sagittal;

double time_cyclic;
vector<double> xstate(2, 0);
vector<double> xstate_e_f(2, 0);
double StepTime1 = 0.45;
double StepTime2 = 0.45;
double StartTime;
double SimTime = 10.0;
double StopTime = SimTime - 0.0;

double dT = 0.0002;
double zc = 0.8; //1.0;//
double hip2ground = 0.65; //1.0;
double angle_bottom;
//double g = 9.81;
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

vector<double> store_traj_LHip(1, 0);
vector<double> store_traj_LKnee(1, 0);
vector<double> store_traj_RHip(1, 0);
vector<double> store_traj_RKnee(1, 0);


vector<double> store_theta_stance(1, 0);
vector<double> store_l_stance(1, 0);
vector<double> store_theta_swing(1, 0);
vector<double> store_l_swing(1, 0);

vector<double> store_theta_left(1, 0);
vector<double> store_r_left(1, 0);
vector<double> store_theta_right(1, 0);
vector<double> store_r_right(1, 0);

vector<double> store_xcom_e(1, 0);
vector<double> store_dxcom_e(1, 0);
double dis_time = 0;
double dis_duration = 0.1;
double dis_A;
bool dis_enable = 0;

//additional code for saggittal, online estimation

int StepIndex;
double samplingTime = 0.004;
int data_set_number = 6;
double com_offset = 0.00;
double fp_wrt_com;

//add noise
double gaussian_e;
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

double vel_f_predict;
double pos_f_predict;
#endif