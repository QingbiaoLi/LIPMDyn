/*****************************************************************************
AlexFilterClass.h

Description:	Header file of AlexFilterClass

@Version:	1.0
@Author:	Chengxu Zhou (chengxu.zhou@iit.it)
@Release:	2014/08/20
@Update:	2014/10/07
*****************************************************************************/
#ifndef ALEXFILTER_CLASS_H// header guards
#define ALEXFILTER_CLASS_H

#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>

class AlexFilterClass
{
public:
	AlexFilterClass();
	// ~AlexFilterClass();

	std::vector<double> applyFilter(int filter_size, const std::vector<double> &raw_data);
	std::vector<double> returnBuffer(int index); // index from front to back, 0 is the newest element
	int returnBufferSize();
	std::vector<double> applyMeanFilter(int filter_size, const std::vector<double> &raw_data);

    Eigen::Vector3d applyFilter(int filter_size, const Eigen::Vector3d &raw_data);
    Eigen::Vector3d applyMeanFilter(int filter_size, const Eigen::Vector3d &raw_data);

private:
	std::vector<std::vector<double> > F_buffer; 
	bool initF;
	// std::vector<std::vector<double> >::iterator it_buffer;	
};
#endif

