/**
* Copyright (C) 2021 All rights reserved.
* @file BKW.h
* @brief  
* @author Yang Kuang (matkuan@nus.edu.sg)
* @date 2021-04-26
* @return 
*/

#ifndef __BKW_H__
#define __BKW_H__

#include "QCalculator.h"
#include "fft3d.h"
#include <type_traits>

template <class c_type, class r_type>
class BKW3D
{
	private:
		int N0;	/// number of grid in one direction
		int N3d;	/// number of grid in whole space

		r_type R;	/// computational domain [-R,R]^3

		c_type* F;	/// function value in real space (currently we use ctype for convenience)
		c_type* Fhat;	/// function value in reciprocal space
		c_type* Q;	/// kernel Q


	public:
		BKW3D(){};
		BKW3D(int _N0, r_type _R):R(_R){initialize(_N0);};
		~BKW3D(){delete[] F; delete[] Fhat;delete [] Q;};

	public:
		void initialize(int);
		void reinit(int);
		void resetR(r_type);
		r_type funExact(r_type t, r_type x, r_type y, r_type z);	/// reference function (used for taking the initial value)
		r_type calEntropy();

		void setInitialValue(r_type t);
		void fft3(bool isInverse=0); /// normally fft3, if isInverse then ifft3

		void RK2(r_type T0, r_type TEnd, r_type dt);		/// run one step RK2
};

#define TEMPLATE template<class c_type, class r_type>
#define THIS BKW3D<c_type,r_type>
TEMPLATE
void THIS::initialize(int _N0)
{
	N0 = _N0;
	N3d = N0*N0*N0;
  F = new c_type[N3d];
  Fhat = new c_type[N3d];
  Q = new c_type[N3d];
};

TEMPLATE
void THIS::reinit(int _N0)
{
	delete [] F;
	delete [] Fhat;
	delete [] Q;
	initialize(_N0);
}


TEMPLATE
void THIS::resetR(r_type _R)
{
	R = _R;
	std::cout << "# Reset computatational domain to ["<<-R << ", "<<R<<"]^3\n";
}

TEMPLATE
r_type THIS::funExact(r_type t, r_type x, r_type y, r_type z)
{
	r_type r = sqrt((x*x+y*y+z*z))/(R*sqrt(3));
	r_type val = 2;
	for(int j = 1; j < 39; ++j)
		val += r_type(j)/741.  * sin(j*r*M_PI);

	return val;
}

TEMPLATE
void THIS::setInitialValue(r_type t)
{
	r_type h = 2*R/r_type(N0);
	int l = 0;
	r_type x,y,z;
	for(int i = 0; i < N0; ++i){
		x = -R + i*h;
		for(int j = 0; j < N0; ++j){
			y = -R + j*h;
			for(int k = 0; k < N0; ++k){
				z = -R + k*h;
				F[l] = c_type(funExact(t,x,y,z)); /// real to complex
				l++;
			}
		}
	}
}

TEMPLATE
r_type THIS::calEntropy()
{
	r_type entropy = 0;
	for(int i = 0; i < N3d; ++i)
		entropy += std::real(F[i]) * log(std::real(F[i]));
	entropy *= std::pow(2*R,3.)/r_type(N3d);
	return entropy;
}


TEMPLATE
void THIS::fft3(bool isInverse)
{
	fft3d f3d;
	if(isInverse)
		f3d.ifft(N0,N0,N0,reinterpret_cast<fftw_complex*>(F),reinterpret_cast<fftw_complex*>(Fhat));
	else
		f3d.fft(N0,N0,N0,reinterpret_cast<fftw_complex*>(F),reinterpret_cast<fftw_complex*>(Fhat));
}

TEMPLATE
void THIS::RK2(r_type T0, r_type TEnd, r_type dt)
{
	bool isInverse = 1;

	std::cout << "# Initialization...\n";
	setInitialValue(T0);
	fft3();

	std::cout << "# QCalculator calculate B...\n";
	QCalculator<c_type,r_type> qcal(N0,0);
	qcal.genB();

	c_type  F1[N3d];	/// tmp vaule in reciprocal space

	std::cout << "# Start the evolution...\n";
	r_type  t = T0, entropy, entropy0 = 1.;
	do{
		t += dt;
		qcal.calQ(Q,Fhat);
		for(int l = 0; l < N0; ++l)
			F1[l] = Fhat[l] + 0.5*dt*Q[l];
		qcal.calQ(Q,F1);
		for(int l = 0; l < N0; ++l)
			Fhat[l] = Fhat[l] + dt*Q[l];

		/// if calculate entropy
		fft3(isInverse);
		entropy = calEntropy();
		std::cout <<std::fixed <<std::setprecision(14)<< t << "\t" << entropy << "\t" << entropy0-entropy  <<  "\n";
		entropy0 = entropy;
	}while(t < TEnd);

}

#undef TEMPLATE
#undef THIS



#endif ///BKW_H__

