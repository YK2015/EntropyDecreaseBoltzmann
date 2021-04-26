/**
* Copyright (C) 2021 All rights reserved.
* @file testBKW3D.cpp
* @brief  
* @author Yang Kuang (matkuan@nus.edu.sg)
* @date 2021-04-26
* @return 
*/


#include "BKW.h"

typedef double  r_type;
typedef std::complex<double>  c_type;
r_type R = 6.62132;
	
int main(int argc, char ** argv)
{
	int N = std::atoi(argv[1]);

	BKW3D<c_type,r_type> bkw3d(N,R);
	bkw3d.RK2(0,0.1,0.001);

	return 0;
}




