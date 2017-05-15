#include <Eigen/Dense>
using namespace Eigen;

extern Matrix3d Rx(double);
extern Matrix3d Ry(double);
extern Matrix3d Rz(double);
extern Matrix3d Rodrigues(Vector3d w,double theta);