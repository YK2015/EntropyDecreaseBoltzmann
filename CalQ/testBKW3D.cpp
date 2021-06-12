/**
* Copyright (C) 2021 All rights reserved.
* @file testBKW3D.cpp
* @brief  
* @author Yang Kuang (matkuan@nus.edu.sg)
* @date 2021-04-26
* @return 
*/


#include "BKW.h"
#include <string>

typedef double  r_type;
typedef std::complex<double>  c_type;
r_type R = 6.62132;
	
int main(int argc, char ** argv)
{
	int N = std::atoi(argv[1]);
	r_type dt = std::atof(argv[2]);
	r_type TEnd = std::atof(argv[3]);

	bool FixEntropy = true;

	std::cout << "##################################################\n"
		<< "# Entropy fixed method for BKW model\n" 
		<< "# spatial domain  :\t[-"<<R<<", "<<R<< "]\n" 
		<< "# number of grid N:\t" << N << "^3\n" 
		<< "# time interval   :\t[0, "<<TEnd <<"]\n" 
		<< "# time step     dt:\t"<<dt<< "\n";
	if(FixEntropy)
		std::cout << "# Entropy fixing stratege is adopted.\n";
	else
		std::cout << "# Entropy fixing stratege is not adopted.\n";

	std::cout << "##################################################\n";

	BKW3D<c_type,r_type> bkw3d(N,R);
	bkw3d.setFixEntropy(FixEntropy);
	//bkw3d.RK2(0,TEnd,dt);
	//bkw3d.CFE(0,TEnd,dt,beta);
	bkw3d.FEuler(0,TEnd,dt);
	std::string filename = "N"+std::to_string(N)+".dt"+std::to_string(dt)+".dat";
	bkw3d.outV1Direction(filename);

	return 0;
}




