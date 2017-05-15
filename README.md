# LIPMDyn
To run this program, you need to:
1. Download [Eigen Library]( http://eigen.tuxfamily.org/index.php?title=Main_Page).

2. Add this package into additinal include directrise,for example:
..\..\eigen; ..\..\header

3. Stucture of the program

* /include:	header files for this program. 

	"MatrixVector.h": user defined matrix based on Eigen library

	"LIPMDynClass.h": key function based on Linear inverted pendulum (LIP) model;

	Optional:

	"AlexFilterClass.h" and "FilterClass.h": Filter file for filtering noise;


* /src: 		header files for this program. 

	"MatrixVector.cpp": user defined matrix based on Eigen library

	"LIPMDynClass.cpp": a class contains key function based on Linear inverted pendulum (LIP) model.

	"LIPMLateralGait.cpp": main function of the simulation of lateral CoM motion and foot placment based on LIP model

	Optional: 

	"AlexFilterClass.cpp" and "FilterClass.cpp": Filter file for filtering noise;
	
* * Most of these are include in LIPMDynClass, ignore them unless you are interested to understand individual function.
	LateralGaitStop.cpp

	LIPMLateralFP.cpp

	LIPM3D.cpp

	CheckLateralGaitSwitch.cpp

	LIPMDyn.cpp

* /data: 		folder to save the log out data, temporary called as "LateralGait.txt" 

  "checkLateralGait.m": Matlab file to read "LateralGait.txt" and plot the data
