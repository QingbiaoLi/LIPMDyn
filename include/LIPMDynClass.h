#ifndef _LIPMDyn_H// header guards
#define _LIPMDyn_H
// 基于线性倒立摆模型的动力学参数预测控制

#include "MatrixVector.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

using namespace std;

class LIPMDynClass
{
public:
	LIPMDynClass();
	virtual ~LIPMDynClass();
	double FPCritical(double zc, double dxf);

	double StopApexPositionControl(double zc, double dxf, double x_apex);
	double StopApexTimeControl(double zc, double dxf, double tau);
	double CrossApexVelocityControl(double zc, double dxf, double v_apex);	
	double CrossApexTimeControl(double zc, double dxf, double tau);
	
	double StateTransitionTime(double zc, double x0, double dx0, double xf, double dxf);
	
	double VelocityAtPositionCV(double zc, double x0, double dx0, double xf);
	double VelocityAtPositionDV(double zc, double x0, double dx0, double xf);
	double PositionAtVelocityCV(double zc, double x0, double dx0, double dxf);

	double LateralGaitSymmetryTimeControl(double zc, double dy0, double StepTime);

	double LateralGaitVelocityControl(double zc, double dym, double dy0, double StepTime);

	double RemainingTimePositionSymmetry(double zc, double x0, double dx0, double xf);
	double RemainingTimeVelocitySymmetry(double zc, double x0, double dx0, double dxf);
	vector<double> StateEvolution(double zc, double x0, double dx0, double time);
	
	double SagittalNominalVelocity(double zc, double step,double Tstep);
	double LateralNominalVelocity(double zc, double stepwidth, double StepTime);

	double LateralNominalSway(double zc, double stepwidth, double StepTime);

	vector<double> MaxAcceleration(double zc, double x0, double dx0, double xf, double p);

	double LateralSwitchTime(double zc, double y0, double dy0, double yf, double p_plus, double p_minus, double t);
	double LateralSwitchEndVelocity(double zc, double y0, double dy0, double yf, double p_plus, double p_minus, double tswitch);

	//Vector3d calcAngularVelocity(double dT, Vector3d ri0,Vector3d ri1);
	//vector<double> foo();

private:
	const double g;
	double Sign(double x);

};
#endif