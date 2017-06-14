#include "LIPMDynClass.h"
#include <iostream>

LIPMDynClass::LIPMDynClass() :
	g(9.81)
{
	;
}

LIPMDynClass::~LIPMDynClass()
{
}


double LIPMDynClass::Sign(double x)
{
	if (x >= 0) return 1;
	else return -1;
}


//vector<double> LIPMDynClass::foo()
//{
//	vector<double> a(3,0);
//	return a;
//}

vector<double> LIPMDynClass::StateEvolution(double zc, double x0, double dx0, double time)
{
	vector<double> statef(2, 0);
	if (time<0)
	{
		cout << "time should not be negative" << endl;
		return statef;
	}
	// it applies to all cases, no matter diverging or converging
	double Tc, xf, dxf;
	Tc = sqrt(zc / g);
	xf = x0*cosh(time / Tc) + dx0*Tc*sinh(time / Tc);
	dxf = x0 / Tc*sinh(time / Tc) + dx0*cosh(time / Tc);

	statef[0] = xf;
	statef[1] = dxf;

	return statef;
}

double LIPMDynClass::FPCritical(double zc, double dxf)
{
	double Tc, CP; // CP is one step capture point
	Tc = sqrt(zc / g);
	CP = Tc*dxf;
	return CP;	// foot placement with respect to the COM
}


double LIPMDynClass::VelocityAtPositionDV(double zc, double x0, double dx0, double xf)
{
	// this function only works for the E>0 (cross apex) 单调函数 Monotonic function in phase plane
	// where xf>x0, dx>0; or xf<x0, dx<0;
	// answers what the velocity is at the future position
	double Tc, dxf, E;
	dxf = 0; // default value
	Tc = sqrt(zc / g);
	E = dx0*dx0 - x0*x0 / (Tc*Tc);

	if (E>0)
	{
		if ((xf >= x0&&dx0>0) || (xf <= x0&&dx0<0))
		{
			dxf = Sign(dx0)*sqrt(dx0*dx0 + g / zc*(xf*xf - x0*x0));
		}
	}
	else // if E<0 then it is NOT monotonic function in phase plane
	{
		cout << "cannot cross Apex, function is not programmed for this case" << endl;
	}
	return dxf;	// 
}


double LIPMDynClass::VelocityAtPositionCV(double zc, double x0, double dx0, double xf)
{
	// this function only works for the ConVerging (CV) phase, Non-Monotonic function in phase plane
	// x, dx 同号 diverging section; x dx 异号, converging phase, dxf 和 dx0 反号
	// answers what the velocity is at the future position
	double Tc, dxf, E, temp;
	dxf = 0; // default value
	Tc = sqrt(zc / g);
	E = dx0*dx0 - x0*x0 / (Tc*Tc);

	if (E<0)
	{
		temp = abs(dx0*dx0 + g / zc*(xf*xf - x0*x0));
		//if (temp<0)
		//{
		//	temp=-temp;
		//}
		if (Sign(x0) != Sign(dx0)) // converging phase
		{
			dxf = -Sign(dx0)*sqrt(temp);
		}
		else if (Sign(x0) == Sign(dx0)) // diverging phase
		{
			dxf = Sign(dx0)*sqrt(temp);
		}
	}
	else // if E>0 then it is cross apex case 
	{
		cout << "cross Apex, function is not programmed for this case" << endl;
	}
	return dxf;	// 
}


double LIPMDynClass::PositionAtVelocityCV(double zc, double x0, double dx0, double dxf)
{
	// this function 
	double xf, E;
	xf = 0; // default value
	E = dx0*dx0 - x0*x0*g / zc;

	if (E<0) // ONLY for converging phase
	{
		xf = Sign(x0)*sqrt(x0*x0 + zc / g*(dxf*dxf - dx0*dx0));
	}
	else // if E>0 then it's in diverging phase, and it's NOT monotonic function in phase plane
	{
		cout << "It can cross Apex, function is not programmed for this case" << endl;
	}
	return xf;	// 
}


double LIPMDynClass::LateralNominalVelocity(double zc, double stepwidth, double StepTime)
{
	// given the step time, and nominal lateral foot placement (step width),
	// calculate what the lateral velocity it should be while cross the middle line 
	double Tc, dym, Tstep, pwidth;
	Tc = sqrt(zc / g);
	Tstep = abs(StepTime);
	pwidth = abs(stepwidth);
	if (Tstep<0.01) // step period 0.01s already doesnt make sense
	{
		Tstep = 0.01;
	}
	//dym = pwidth/Tc*(exp(Tstep/Tc)-exp(-Tstep/Tc))/(exp(Tstep/Tc)+exp(-Tstep/Tc)+2.0);
	dym = pwidth / Tc*(exp(Tstep / Tc) - 1.0) / (exp(Tstep / Tc) + 1.0);
	return dym;
}

double LIPMDynClass::LateralNominalSway(double zc, double stepwidth, double StepTime)
{
	// given the nominal step width (COP), and step period Tstep
	// this function calculates what the maximum sway distance is
	// from the neutral center
	double Tc, ys, Tstep, pwidth;
	Tc = sqrt(zc / g);
	pwidth = abs(stepwidth);
	Tstep = abs(StepTime);
	if (Tstep>87 * Tc)
	{
		Tstep = 87 * Tc; // prevent numerical overflow
	}
	ys = pwidth*(1.0 - 2.0*exp(Tstep / Tc / 2) / (exp(Tstep / Tc) + 1.0));
	return ys;
}

double LIPMDynClass::LateralGaitVelocityControl(double zc, double dym, double dy0, double StepTime)
{
	// given the step time, target (nominal) velocity dym and initial lateral velocity dy0 
	// cross the mid line,
	// calculate what the lateral foot placement should be for recovering velocity
	double Tc, dyf, yfp, Tstep;
	Tc = sqrt(zc / g);
	Tstep = abs(StepTime);
	dyf = -Sign(dy0)*dym;
	if (Tstep<0.01) // step period 0.01s already doesnt make sense
	{
		Tstep = 0.01;
	}
	yfp = (Tc*dy0*(exp(Tstep / Tc) + exp(-Tstep / Tc)) - 2.0*Tc*dyf) / (exp(Tstep / Tc) - exp(-Tstep / Tc));
	return yfp;
}



double LIPMDynClass::LateralGaitSymmetryTimeControl(double zc, double dy0, double StepTime)
{
	// for lateral foot placement control, definitely converging phase
	// by default the entry position is the same as exit position.
	double Tc, xfp;
	Tc = sqrt(zc / g);
	xfp = Tc*dy0*(exp(StepTime / Tc) + 1.0) / (exp(StepTime / Tc) - 1.0);
	return xfp;
}


double LIPMDynClass::RemainingTimePositionSymmetry(double zc, double x0, double dx0, double xf)
{
	// this function uses existing functions to rule out how much time left to reach a target COM position
	double t, dxf;
	dxf = VelocityAtPositionCV(zc, x0, dx0, xf); // this is a general function		

	t = StateTransitionTime(zc, x0, dx0, xf, dxf);

	return t;	// foot placement with respect to the COM
}

double LIPMDynClass::RemainingTimeVelocitySymmetry(double zc, double x0, double dx0, double dxf)
{
	// this function uses existing functions to rule out how much time left to reach a target COM position
	double t, xf;
	xf = PositionAtVelocityCV(zc, x0, dx0, dxf); // this is a general function		

	t = StateTransitionTime(zc, x0, dx0, xf, dxf);

	return t;	// foot placement with respect to the COM
}


double LIPMDynClass::StopApexPositionControl(double zc, double dxf, double x_apex)
{
	/*	here the velocity at Apex is assumed to be zero
	zc:		COM height
	dxf:	the end velocity before the new stance foot
	x_apex: the angle where the new stance leg shall stop at, typically within 30 deg
	*/
	double Tc, xfp;
	Tc = sqrt(zc / g);
	xfp = Sign(dxf)*Tc*sqrt(dxf*dxf + x_apex*x_apex*g / zc);
	return xfp;	// with respect to the COM position
				/*this answers where you should put your leg if you want to stop at
				a specific apex position with respect to new support foot, this is a general equation
				*/
}


double LIPMDynClass::CrossApexVelocityControl(double zc, double dxf, double v_apex)
{
	/*	for sagittal gait control
	dxf:	end velocity before impact
	xf:		end position before impact
	zc:		nominal COM height, natural frequency wn is replaced by sqrt(g/zc)
	v_apex:	desired angular velocity at apex point
	range of v_apex:  0=<v_apex<=dxf
	*/
	if (abs(v_apex)>abs(dxf))	// desired v_apex cannot be bigger than dxf
	{
		v_apex = abs(dxf)*v_apex / abs(v_apex);	// scale down
	}
	double Tc, xfp;
	Tc = sqrt(zc / g);
	xfp = Sign(dxf)*Tc*sqrt(dxf*dxf - v_apex*v_apex);
	return xfp;
	//xfp minimum is 0, when v_apex is the same as dxf. if v_apex
	//is large than dxf then the equation is not valid, since xfp
	//is the same sign as xf(n) and there doesnt exist apex anymore
	// 公式仅在 dxf=<v_apex<=0的情况下成立(能量不能为负)
}

double LIPMDynClass::StopApexTimeControl(double zc, double dxf, double tau)
{
	/*	Tc:  sqrt(zc/g)
	dxf: final velocity of stance leg before a new step
	tau: the time since new stance to the apex
	range of tau: 0=<tau<oo
	*/
	double Tc, exponent, D, xfp;
	D = 0;
	Tc = sqrt(zc / g);
	exponent = 2 * abs(tau) / Tc;

	if (exponent>87) //exp(87) is the max that a double number can represent
	{
		D = 1.0;
	}
	else
	{
		D = (exp(exponent) + 1.0) / (exp(exponent) - 1.0);
	}
	xfp = D*Tc*dxf;
	return xfp;
	/*this answers where you should put your leg if you want a specific time to reach the apex
	and this FE is always beyond the critical FE so the COM always returns, motion is bounded
	*/
}

double LIPMDynClass::CrossApexTimeControl(double zc, double dxf, double tau)
{
	/*	for sagittal gait control
	Tc:  sqrt(zc/g)
	dxf: final angular velocity of stance leg
	xf: final angle of stance leg
	tau: the time since new stance to the apex
	range of tau: 0=<tau<oo
	*/
	double Tc, exponent, D, xfp;
	D = 0;
	Tc = sqrt(zc / g);
	exponent = 2 * abs(tau) / Tc;

	if (exponent>87) //exp(87) is the max that a double number can represent
	{
		D = 1.0;
	}
	else
	{
		D = (exp(exponent) - 1.0) / (exp(exponent) + 1.0);
	}
	xfp = D*Tc*dxf;
	return xfp;

	/*this answers where you should put your leg if you want a specific time to cross the apex
	and this FE is always behind the critical FE so the COM always passes apex
	*/
}


double LIPMDynClass::StateTransitionTime(double zc, double x0, double dx0, double xf, double dxf)
{
	double a1, a2, b1, b2, c1, c2, Tc, t;
	Tc = sqrt(zc / g);
	a1 = x0 + Tc*dx0;
	b1 = xf + Tc*dxf;
	a2 = xf - Tc*dxf;
	b2 = x0 - Tc*dx0;
	t = 0;

	if (a1 != 0 && a2 == 0)
	{
		c1 = b1 / a1;
		if (c1>0)
		{
			t = Tc*log(c1);
		}
		else
		{
			cout << "error: c1<0" << endl;
			system("pause");
		}
	}
	else if (a2 != 0 && a1 == 0)
	{
		c2 = b2 / a2;
		if (c2>0)
		{
			t = Tc*log(c2);
		}
		else
		{
			cout << "error: c2<0" << endl;
			system("pause");
		}
	}
	else if (a1 != 0 && a2 != 0)
	{
		c1 = b1 / a1;
		c2 = b2 / a2;
		if (c1>0)
		{
			t = Tc*log(c1);
			//printf("\n%f\n",t);
		}
		else
		{
			cout << "error: c1<0" << endl;
			system("pause");
		}
		if (c2>0)
		{
			t = Tc*log(c2);
			//printf("\n%f\n",t);	
		}
		else
		{
			cout << "error: c2<0" << endl;
			system("pause");
		}
	}
	else
	{
		cout << "error:both c1 c2 zero" << endl;
		system("pause");
	}

	return t;
}

double LIPMDynClass::SagittalNominalVelocity(double zc, double step, double Tstep)
{
	double Tc, Tabs, dx;
	Tc = sqrt(zc / g);
	Tabs = abs(Tstep);
	if (Tabs<0.01) // step period 0.01s already doesnt make sense
	{
		Tabs = 0.01;
	}
	dx = (2 * step / Tabs)*exp(Tabs / Tc) / (exp(2 * Tabs / Tc) - 1.0);
	return dx; // the velocity already contains the sign depending on step length's sign
}

vector<double> LIPMDynClass::MaxAcceleration(double zc, double x0, double dx0, double xf, double p)
{
	// given the initial state, and final position allowed, before hitting the friction cone
	// what is the max speed you can get, and what is the min time spent
	// p is the max position you can shift the COP backward to accelerate 
	double x0new, xfnew, dxfmax, tmin;
	vector<double> boundary(2, 0);
	x0new = x0 - p;
	xfnew = xf - p;

	dxfmax = Sign(x0new)*sqrt(dx0*dx0 + g / zc*(xfnew*xfnew - x0new*x0new));
	tmin = StateTransitionTime(zc, x0new, dx0, xfnew, dxfmax);

	boundary[0] = dxfmax;
	boundary[1] = tmin;
	return boundary; // return the max speed you can get at specified position, and this takes the minimum time
}


double LIPMDynClass::LateralSwitchTime(double zc, double y0, double dy0, double yf, double p_plus, double p_minus, double t)
{
	// given the initial state, and final position to be reached
	// the max and min COP to accelerate, what is the switching time
	// default direction bool dir=1, move to left. t 是输入的总时间
	double Tc, x0, dx0, xf, tswitch, C, D, t2; // x is just the unknown
	Tc = sqrt(zc / g);
	x0 = y0 - p_minus;
	xf = yf - p_plus;
	dx0 = dy0;
	double temp;
	temp = t / Tc;	// actually if t>1.5s or 1.8s, the COM almost stops at support foot,
	if (temp>80)
	{
		temp = 80;
		cout << "overflow of exp" << endl;
		system("pause");
	}

	C = ((x0 + Tc*dx0)*exp(temp) + (x0 - Tc*dx0)*exp(-temp) - 2.0*xf) / (p_plus - p_minus);
	D = C / 2.0;
	if (D<1.0)
	{
		D = 1.0;
		cout << "COM pass over apex, obital energy>0. Precondition violated!" << endl;
	}
	t2 = Tc * log(D + sqrt(D*D - 1.0)); // 反双曲函数, 第二相的过程时间
	tswitch = t - t2;   // x is the duration of 2nd phase

	return tswitch; // return time of switching phase
}


double LIPMDynClass::LateralSwitchEndVelocity(double zc, double y0, double dy0, double yf, double p_plus, double p_minus, double tswitch)
{
	/*this function follows up LateralSwitchTime, to predict what is the end velocity in the end*/
	double x0, dx0, xf, dxf; // x is just the unknown
	vector<double> ystate(2, 0);
	//Tc=sqrt(zc/g);
	x0 = y0 - p_minus;
	dx0 = dy0;
	xf = yf - p_plus;

	ystate = StateEvolution(zc, x0, dx0, tswitch);

	dxf = VelocityAtPositionCV(zc, ystate[0] + (p_minus - p_plus), ystate[1], xf);

	return dxf; // return end velocity at the end of 2nd phase, validated, very accurate
}



//Vector3d LIPMDynClass::calcAngularVelocity(double dT, Vector3d ri0,Vector3d ri1)
//{
//	// ri0 and ri1 are the previous and current position vector
//	static double a, b, c, w_norm;
//	static Vector3d C, w, axis;
//	C = ri1-ri0;
//	a = ri0.norm();
//	b = ri1.norm();
//	c = C.norm();
//
//	w_norm = acos( (a*a+b*b-c*c)/(2*a*b) ) / dT;  //余弦定理
//	axis = ri0.cross(ri1); // axis of angular velocity
//	w = axis/axis.norm() * w_norm;   // angular velocity vector
//
//	return w;
//}

//-----------------  chengxu added 20141105
double LIPMDynClass::LateralNominalVelSSend(double zc, double stepwidth, double dym, double StepTime)
{
	// given the step time, and nominal lateral foot placement (step width),
	// calculate what the lateral velocity it should be at SS end 
	double Tc, Tstep, pwidth, dymSS;
	Tc = sqrt(zc / g);
	Tstep = abs(StepTime);
	pwidth = abs(stepwidth);
	if (Tstep<0.01) // step period 0.01s already doesnt make sense
	{
		Tstep = 0.01;
	}
	dymSS = pwidth / Tc*sinh(-0.1*Tstep / Tc) + dym*cosh(-0.1*Tstep / Tc);
	return dymSS;
}

double LIPMDynClass::deltaLateralEndVel(double vd0, double vdf, double deltaV0)
{
	double deltaVf;
	if (vdf >= 0) {
		deltaVf = sqrt(4 * vdf*vdf - 4 * (deltaV0*deltaV0)) - 2 * vdf;
		//deltaVf = sqrt(4*vdf*vdf+4*(2*vd0*deltaV0+deltaV0*deltaV0))-2*vdf;
	}
	else {
		deltaVf = -sqrt(4 * vdf*vdf - 4 * (deltaV0*deltaV0)) - 2 * vdf;
		//deltaVf = -sqrt(4*vdf*vdf+4*(2*vd0*deltaV0+deltaV0*deltaV0))-2*vdf;
	}
	deltaVf = abs(0.5*deltaVf);
	//if(vd0*deltaV0 >= 0){
	//	deltaVf = Sign(vdf)*deltaVf;
	//}
	//else{
	//	deltaVf = -Sign(vdf)*deltaVf;
	//}
	return deltaVf;
}

double LIPMDynClass::deltaLateralPos(double zc, double deltaVel, double StepTime)
{
	double Tc, deltaPos;
	Tc = sqrt(zc / g);
	// footplacement error = deltaVel*Tc* coth(\tau_{s})  
	deltaPos = deltaVel*Tc*(exp(2 * StepTime / Tc) + 1.0) / (exp(2 * StepTime / Tc) - 1.0);
	return deltaPos;
}

double LIPMDynClass::SagittallPos(double zc, double dx0, double dxm, double StepTime)
{
	//the use of SagittallPos, given the current measured velocity (dx0), to answer how much you adjust your footplacement (relative change of the foot placement position,deltaPos), 
	//  how fast or within what time (StepTime) for the COM velocity to reach the desired the one (dxm)
	double Tc, deltaPos;
	Tc = sqrt(zc / g);

	deltaPos = 2 * Tc*dxm / (exp(-StepTime / Tc) - exp(StepTime / Tc)) + dx0*Tc*(exp(2 * StepTime / Tc) + 1.0) / (exp(2 * StepTime / Tc) - 1.0);
	return deltaPos;
}


//Vector3d LIPMDynClass::calcAngularVelocity(double dT, Vector3d ri0,Vector3d ri1)
//{
//	// ri0 and ri1 are the previous and current position vector
//	static double a, b, c, w_norm;
//	static Vector3d C, w, axis;
//	C = ri1-ri0;
//	a = ri0.norm();
//	b = ri1.norm();
//	c = C.norm();
//
//	w_norm = acos( (a*a+b*b-c*c)/(2*a*b) ) / dT;  //余弦定理
//	axis = ri0.cross(ri1); // axis of angular velocity
//	w = axis/axis.norm() * w_norm;   // angular velocity vector
//
//	return w;
//}
