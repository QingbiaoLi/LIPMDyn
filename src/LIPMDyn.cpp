/* Programmed by: Zhibin Li*/
// this program tests the LIPM dynamcis functions
#include <iostream>
#include "LIPMDynClass.h"

#define DEGTORAD(x)  x*M_PI/180.0
#define RADTODEG(x)  x*180.0/M_PI

LIPMDynClass LIPM;

void main()
{
	float ln, dthetaf, thetaf, theta_apex, tau, w_apex, x, t;
	ln=0.5; 
	dthetaf=1.0;
	thetaf=M_PI/6.0; //30 degree
	w_apex=0.5/ln;

	//system("pause");

	x=LIPM.CPCritical(ln,  dthetaf);
	cout<< RADTODEG(x) <<endl;

	theta_apex=DEGTORAD(-10);
	x=LIPM.StopApexFPControl(ln, dthetaf, thetaf, theta_apex);
	cout<< RADTODEG(x) <<"<"<< RADTODEG(theta_apex) <<endl;
	//system("pause");

	tau=100;
	x=LIPM.StopApexTimeControl(ln, dthetaf, thetaf, tau);
	cout<< RADTODEG(x) <<endl;
	//system("pause");

	tau=100;
	x=LIPM.CrossApexTimeControl(ln, dthetaf,thetaf, tau);
	cout<< RADTODEG(x) <<endl;

	w_apex = 0.0*dthetaf;
	x=LIPM.CrossApexVelocityControl(ln, dthetaf,thetaf, w_apex);
	cout<< RADTODEG(x) <<endl;

	t=LIPM.TransitionTime(ln, 0, dthetaf, thetaf,  1.5*dthetaf);
	cout<< endl;
	cout<< t <<endl;

	system("pause");


	
}

