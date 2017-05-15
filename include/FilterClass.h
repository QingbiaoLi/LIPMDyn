#ifndef _Filter_H// header guards
#define _Filter_H

#define MAX_FILTER_LENGTH 8
//typedef double double; //number type
#define _USE_MATH_DEFINES
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

class FilterClass 
{

public:
	//Define Variables
	//************************************************************************************
	int Ncoeff; 					// Number of coefficients
	double x[MAX_FILTER_LENGTH]; 	// Circular buffer of incoming values to be filtered
	double y[MAX_FILTER_LENGTH];		// Result and accumulators
	double a[MAX_FILTER_LENGTH];		// Y coefficient
	double b[MAX_FILTER_LENGTH];		// X coefficient
	FilterClass();
	/*End Variables definition
	************************************************************************************
	************************************************************************************
	//Begin function heders
	************************************************************************************/
	void clear_filter();  //set values to 0
	void int_diff_filter(double X);
	/*Filter Functions
	T is the sampleTime in second	
	cutoff is the cutoff frequency in Hz
	N is the order of the filter*/

	void least_squares_filter(double T, int N); 
	void moving_average_filter(int N);
	void butterworth(double T, double cutoff, int N);
	void butterDifferentiator(double T, double cutoff, int N);

	//Use Filter Function
	std::vector<double> applyFilter(const std::vector<double>& X);
	Eigen::Vector3d ApplyButterworth(double T, double cutoff, int N, Eigen::Vector3d raw_data);
	std::vector<double> ApplyButterworth(double T, double cutoff, int N, const std::vector<double> &raw_data);

	double applyFilter(double input);
	double differentiator(double input);
	//end Functions

    inline void Init(){InitF=false;};

private:
	bool InitF;
	std::vector<std::vector<double> > Xvec;
	std::vector<std::vector<double> > Yvec;
};
#endif
