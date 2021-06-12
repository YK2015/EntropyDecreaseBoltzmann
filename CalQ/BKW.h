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
#include "TimeWatch.h"
#include "EntropyFix.h"
#include <cmath>
#include <iomanip>
#include <string>
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

		bool FixEntropy=false;

	public:
		BKW3D(){};
		BKW3D(int _N0, r_type _R):R(_R){initialize(_N0);};
		~BKW3D(){delete[] F; delete[] Fhat;delete [] Q;};

	public:
		void initialize(int);
		void reinit(int);
		void resetR(r_type);
		void setFixEntropy(bool _fixEntropy=true){FixEntropy=_fixEntropy;};
		r_type funExact(r_type t, r_type x, r_type y, r_type z);	/// reference function (used for taking the initial value)
		r_type calEntropy();
		r_type flogf(r_type);

		void setInitialValue(r_type t);
		void fft3(bool isInverse=0); /// normally fft3, if isInverse then ifft3

		void FEuler(r_type T0, r_type TEnd, r_type dt);		/// run one step RK2
		void RK2(r_type T0, r_type TEnd, r_type dt);		/// run one step RK2
		void SSP3(r_type T0, r_type TEnd, r_type dt);		/// run one step RK2
		void RK4(r_type T0, r_type TEnd, r_type dt);		/// run one step RK2
		void DUFORT(r_type T0, r_type TEnd, r_type dt, r_type beta);		
		void CFE(r_type T0, r_type TEnd, r_type dt, r_type beta);		

		void outV1Direction(std::string);
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
#if 1
	r_type r1 = (x+R)/(2.*R);/// r1 \in [0,1]
	r_type r2 = (y+R)/(2.*R);
	r_type r3 = (z+R)/(2.*R);
	r_type r = sqrt((r1-0.5)*(r1-0.5) +(r2-0.5)*(r2-0.5) +(r3-0.5)*(r3-0.5))/sqrt(3);//+ r2*r2 + r3*r3);
	r_type val = 0.0000001;
	
//	if(r < 0.25)
	{
		val = 3.2;
		//int j = 11.;
		for(int j = 0; j <11; ++j)
			//val += 0.01 * r_type(rand () % 5);
			val += std::abs(r_type(j))/55 * sin(j*(r1-0.5)*M_PI) + std::abs(r_type(j))/55.  * sin(j*(r2-0.5)*M_PI) + std::abs(r_type(j))/55.  * sin(j*(r3-0.5)*M_PI);//* cos(j*r2*M_PI)* sin(j*r3*M_PI);
	}

	return val;
#elif 
	r_type v2 = (x*x+y*y+z*z);
	r_type s = 1.-0.4*exp(-t/6.);
	return 1./std::pow(2.*M_PI*s,1.5)*exp(-v2/(2.*s)) * ((5.*s-3.)/(2.*s)+ (1.-s)/(2.*s*s)*v2);
#endif

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
	for(int i = 0; i < N3d; ++i){
	//	if(std::real(F[i]) > 1.0e-12)
	//	entropy += std::real(F[i]) * log(std::real(F[i]))-std::real(F[i]);
	entropy += flogf(std::real(F[i]));
		
	}
	entropy *= std::pow(2*R,3.)/r_type(N3d);
	return entropy;
}
TEMPLATE
r_type THIS::flogf(r_type f)
{
	r_type val = 0.;
	if(std::abs(f) < 1.0e-6) val = 0.;
	else val = f*log(f) - f;

	return val;
}

TEMPLATE
void THIS::outV1Direction(std::string filename)
{
	r_type x;
	int l ;
	std::ofstream write(filename);
	double h = 2.*R/double(N0);
	write << "# x\tf(t)\tf_0\n";
	r_type fval;
	for(int i = 0; i < N0; ++i){
		x = -R + i *h; 
		fval = funExact(0,x,h/2,h/2);
		l =  (i*N0+N0/2)*N0 + N0/2;
		write << std::fixed << std::setprecision(8) << x  << "\t" << std::real(F[l]) << "\t" << fval << "\n";
	}
	write.close();

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
void THIS::FEuler(r_type T0, r_type TEnd, r_type dt)
{
	bool isInverse = 1;

	std::cout << "# Initialization...\n";
	setInitialValue(T0);
	fft3();

	std::cout << "# QCalculator calculate B...\n";
	QCalculator<c_type,r_type> qcal(N0,0);
	{
		TimeWatch<r_type> tw;
		qcal.genB();
	}

	c_type  F1[N3d];	/// tmp vaule in reciprocal space
	r_type * Ftmp; Ftmp = new r_type [N3d]; 
	EntropyFix<int,r_type> efix(N3d); 

	std::cout << "# Start the evolution using Forward Euler...\n";
	std::cout << "# Step\tt\tentropy_{n}\tentropy_{n-1}-entropy_{n}\tfixedEntropy\tfixedDiff\n";
	TimeWatch<r_type> tw;
	r_type  t = T0, entropy, entropy0 = calEntropy();
	std::cout <<std::fixed <<std::setprecision(8) <<   "00\t"<< t << "\t" << entropy0 << "\t-\t-\t-\t\n";
	int step = 1;
	std::ofstream writeF("F.dat");
	do{
		t += dt;
		qcal.calQ(Q,Fhat);
		for(int l = 0; l < N3d; ++l)
			Fhat[l] = Fhat[l] + dt*Q[l];

		/// if calculate entropy
		fft3(isInverse);
		entropy = calEntropy();
		if(step>1){
			std::cout <<std::fixed <<std::setprecision(8)<< step << "\t" << t << "\t" << entropy << "\t" << entropy0-entropy  <<  "\t";
			if(FixEntropy && entropy - entropy0 > 1.0e-10){
				for(int i =0;i<N3d;++i) Ftmp[i]=std::real(F[i]);
				efix.FixSol(N3d,Ftmp,entropy0/(std::pow(2*R,3)));
				for(int i =0;i<N3d;++i) F[i]=c_type{Ftmp[i]};
				entropy = calEntropy();
				std::cout << entropy << "\t" << entropy0 - entropy;
			}
			else 
				std::cout << entropy << "\t" << entropy0 - entropy;
			std::cout << std::endl;
		}
		else
			std::cout <<std::fixed <<std::setprecision(8) << step << "\t"<< t << "\t" << entropy << "\t-\t-\t-\t\n";
		entropy0 = entropy; step ++;
		
		writeF << t << "\t";
		for(int i =0;i < N3d;++i) 
			writeF  << std::real(F[i]) << "\t";
		writeF << "\n";
	}while(t < TEnd);
	writeF.close();

}


TEMPLATE
void THIS::CFE(r_type T0, r_type TEnd, r_type dt, r_type beta)
{
	bool isInverse = 1;

	std::cout << "# Initialization...\n";
	setInitialValue(T0);
	fft3();

	std::cout << "# QCalculator calculate B...\n";
	QCalculator<c_type,r_type> qcal(N0,0);
	{
		TimeWatch<r_type> tw;
		qcal.genB();
	}

	c_type  F1[N3d];	/// tmp vaule in reciprocal space
	r_type * Ftmp; Ftmp = new r_type [N3d]; 
	EntropyFix<int,r_type> efix(N3d); 

	std::cout << "# Start the evolution using C-Forward Euler...\n";
	std::cout << "# t\taverage\tentropy_{n}\tentropy_{n-1}-entropy_{n}\n";
	TimeWatch<r_type> tw;
	r_type  t = T0, entropy, entropy0 = 1.;
	int step = 0;
	do{
		t += dt;
		qcal.calQ(Q,Fhat);
		for(int l = 0; l < N3d; ++l)
			Fhat[l] = Fhat[l] + (1.-exp(-beta*dt))/beta * Q[l];

		/// if calculate entropy
		fft3(isInverse);

		r_type sum = 0.;
		for(int i =0; i < N3d; ++i) sum += std::real(F[i]);
		sum /= N3d;
		entropy = calEntropy();
		if(step){
			std::cout <<std::fixed <<std::setprecision(8)<< t  << "\t" << sum << "\t" << entropy << "\t" << entropy0-entropy  <<  "\t";
			if(FixEntropy&&entropy - entropy0 > 1.0e-10){
				for(int i =0;i<N3d;++i) Ftmp[i]=std::real(F[i]);
				efix.FixSol(N3d,Ftmp,entropy0/(std::pow(2*R,3)));
				for(int i =0;i<N3d;++i) F[i]=c_type{Ftmp[i]};

				sum=0.;
				for(int i =0; i < N3d; ++i) sum += std::real(F[i]);
				sum /= N3d;

				entropy = calEntropy();
				std::cout << entropy << "\t" << entropy0 - entropy << "\t" << sum;
			}
			std::cout << std::endl;
		}
		else
			std::cout <<std::fixed <<std::setprecision(8)<< t << "\t" << entropy << "\t-\n";
		entropy0 = entropy; step ++;
	}while(t < TEnd);

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
	{
		TimeWatch<r_type> tw;
		qcal.genB();
	}

	c_type  F1[N3d];	/// tmp vaule in reciprocal space
	r_type * Ftmp; Ftmp = new r_type [N3d]; 
	EntropyFix<int,r_type> efix(N3d); 

	std::cout << "# Start the evolution using RK2...\n";
	std::cout << "# t\tentropy_{n}\tentropy_{n-1}-entropy_{n}\n";
	TimeWatch<r_type> tw;
	r_type  t = T0, entropy, entropy0 = 1.;
	int step = 0;
	do{
		t += dt;
		qcal.calQ(Q,Fhat);
		for(int l = 0; l < N3d; ++l)
			F1[l] = Fhat[l] + 0.5*dt*Q[l];
		qcal.calQ(Q,F1);
		for(int l = 0; l < N3d; ++l)
			Fhat[l] = Fhat[l] + dt*Q[l];

		/// if calculate entropy
		fft3(isInverse);

		r_type sum = 0.;
		for(int i =0; i < N3d; ++i) sum += std::real(F[i]);
		sum /= N3d;
		entropy = calEntropy();
		if(step){
			std::cout <<std::fixed <<std::setprecision(8)<< t  << "\t" << sum<< "\t" << entropy << "\t" << entropy0-entropy  <<  "\t";
//			if(entropy - entropy0 > 1.0e-10){
			if(0){

				for(int i =0;i<N3d;++i) Ftmp[i]=std::real(F[i]);
				r_type sum = 0.;
				for(int i =0; i < N3d; ++i) sum += std::real(F[i]);
				sum /= N3d;
				//std::cout << "\t\tBeforeFix:: aver = " << sum << "\n";

				efix.FixSol(N0*N0*N0,Ftmp,entropy0/(std::pow(2*R,3)));
				for(int i =0;i<N3d;++i) F[i]=c_type{Ftmp[i]};

				sum = 0.;
				for(int i =0; i < N3d; ++i) sum += std::real(F[i]);
				sum /= N3d;

				entropy = calEntropy();
				fft3();
				std::cout << entropy << "\t" << entropy0 - entropy << "\t" << sum;
			}
			std::cout << std::endl;
		}
		else
			std::cout <<std::fixed <<std::setprecision(8)<< t << "\t" << sum<< "\t" << entropy << "\t-\n";
		entropy0 = entropy; step ++;
	}while(t < TEnd);

}


TEMPLATE
void THIS::SSP3(r_type T0, r_type TEnd, r_type dt)
{
	bool isInverse = 1;

	std::cout << "# Initialization...\n";
	setInitialValue(T0);
	fft3();

	std::cout << "# QCalculator calculate B...\n";
	QCalculator<c_type,r_type> qcal(N0,0);
	{
		TimeWatch<r_type> tw;
		qcal.genB();
	}


	EntropyFix<int,r_type> efix(N3d); 
		
	c_type  F1[N3d],F2[N3d];	/// tmp vaule in reciprocal space
	r_type * Ftmp; Ftmp = new r_type [N3d]; 

	std::cout << "# Start the evolution using SSP3...\n";
	std::cout << "# t\tentropy_{n}\tentropy_{n-1}-entropy_{n}\tNew entropy\tnew entropy_{n-1}-entropy_{n}\n";
	TimeWatch<r_type> tw;
	r_type  t = T0, entropy, entropy0 = 1.;
	int step = 0;
	do{
		t += dt;
		qcal.calQ(Q,Fhat);
		for(int l = 0; l < N3d; ++l)
			F1[l] = Fhat[l] + dt*Q[l];
		qcal.calQ(Q,F1);
		for(int l = 0; l < N3d; ++l)
			F2[l] = 0.75*Fhat[l] + 0.25*F1[l] + 0.25*dt*Q[l];
		qcal.calQ(Q,F2);
		for(int l = 0; l < N3d; ++l)
			Fhat[l] = Fhat[l]/3. + 2.*F2[l]/3. + 2./3.*dt*Q[l];

		/// if calculate entropy
		fft3(isInverse);

		entropy = calEntropy();
		if(step){
			std::cout <<std::fixed <<std::setprecision(8)<< t << "\t" << entropy << "\t" << entropy0-entropy  <<  "\t";
			if(entropy - entropy0 > 1.0e-10){
				
					for(int i =0;i<N3d;++i) Ftmp[i]=std::real(F[i]);
					efix.FixSol(N0*N0*N0,Ftmp,entropy0/(std::pow(2*R,3)));
					for(int i =0;i<N3d;++i) F[i]=c_type{Ftmp[i]};
				
				entropy = calEntropy();
				std::cout << entropy << "\t" << entropy0 - entropy;
			}
			std::cout << std::endl;
		}
		else
			std::cout <<std::fixed <<std::setprecision(8)<< t << "\t" << entropy << "\t-\n";
		entropy0 = entropy; step ++;
	}while(t < TEnd);

	delete [] Ftmp;
}


TEMPLATE
void THIS::RK4(r_type T0, r_type TEnd, r_type dt)
{
	bool isInverse = 1;

	std::cout << "# Initialization...\n";
	setInitialValue(T0);
	fft3();

	std::cout << "# QCalculator calculate B...\n";
	QCalculator<c_type,r_type> qcal(N0,0);
	{
		TimeWatch<r_type> tw;
		qcal.genB();
	}


	EntropyFix<int,r_type> efix(N3d); 
		
	c_type  F1[N3d],F2[N3d],F3[N3d];	/// tmp vaule in reciprocal space
	c_type  Q1[N3d],Q2[N3d],Q3[N3d];	/// tmp vaule in reciprocal space
	r_type * Ftmp; Ftmp = new r_type [N3d]; 

	std::cout << "# Start the evolution using RK4...\n";
	std::cout << "# t\tentropy_{n}\tentropy_{n-1}-entropy_{n}\tNew entropy\tnew entropy_{n-1}-entropy_{n}\n";
	TimeWatch<r_type> tw;
	r_type  t = T0, entropy, entropy0 = 1.;
	int step = 0;
	do{
		t += dt;
		qcal.calQ(Q,Fhat);
		for(int l = 0; l < N3d; ++l)
			F1[l] = Fhat[l] + 0.5*dt*Q[l];
		qcal.calQ(Q1,F1);
		for(int l = 0; l < N3d; ++l)
			F2[l] = Fhat[l] +  0.5*dt*Q1[l];
		qcal.calQ(Q2,F2);
		for(int l = 0; l < N3d; ++l)
			F3[l] = Fhat[l] +  dt*Q2[l];
		qcal.calQ(Q3,F3);

		for(int l = 0; l < N3d; ++l)
			Fhat[l] = Fhat[l] + dt/6. * (Q[l]+2.*Q1[l]+2.*Q2[l]+Q3[l]);

		/// if calculate entropy
		fft3(isInverse);

		entropy = calEntropy();
		if(step){
			std::cout <<std::fixed <<std::setprecision(8)<< t << "\t" << entropy << "\t" << entropy0-entropy  <<  "\t";
			if(entropy - entropy0 > 1.0e-10){
				
					for(int i =0;i<N3d;++i) Ftmp[i]=std::real(F[i]);
					efix.FixSol(N0*N0*N0,Ftmp,entropy0/(std::pow(2*R,3)));
					for(int i =0;i<N3d;++i) F[i]=c_type{Ftmp[i]};
				
				entropy = calEntropy();
				std::cout << entropy << "\t" << entropy0 - entropy;
			}
			std::cout << std::endl;
		}
		else
			std::cout <<std::fixed <<std::setprecision(8)<< t << "\t" << entropy << "\t-\n";
		entropy0 = entropy; step ++;
	}while(t < TEnd);

	delete [] Ftmp;
}


TEMPLATE
void THIS::DUFORT(r_type T0, r_type TEnd, r_type dt, r_type beta)
{
	bool isInverse = 1;

	std::cout << "# Initialization...\n";
	setInitialValue(T0);
	fft3();

	std::cout << "# QCalculator calculate B...\n";
	QCalculator<c_type,r_type> qcal(N0,0);
	{
		TimeWatch<r_type> tw;
		qcal.genB();
	}

	c_type  F1[N3d], F_old[N3d];	/// tmp vaule in reciprocal space
	for(int i = 0;i < N3d; ++i) F_old[i] = Fhat[i];	/// at step n-1
	r_type * Ftmp; Ftmp = new r_type [N3d]; 
	EntropyFix<int,r_type> efix(N3d); 

	std::cout << "# Start the evolution using RK2...\n";
	std::cout << "# t\tentropy_{n}\tentropy_{n-1}-entropy_{n}\n";
	TimeWatch<r_type> tw;
	r_type  t = T0, entropy, entropy0 = 1.;
	r_type dtRK2 = dt/8;
	int step = 0;
	do{
		t += dtRK2;
		qcal.calQ(Q,Fhat);
		for(int l = 0; l < N3d; ++l)
			F1[l] = Fhat[l] + 0.5*dt*Q[l];
		qcal.calQ(Q,F1);
		for(int l = 0; l < N3d; ++l)
			Fhat[l] = Fhat[l] + dt*Q[l];

		/// if calculate entropy
		fft3(isInverse);

		r_type sum = 0.;
		for(int i =0; i < N3d; ++i) sum += std::real(F[i]);
		sum /= N3d;
		entropy = calEntropy();
		if(step){
			std::cout <<std::fixed <<std::setprecision(8)<< t  << "\t" << sum<< "\t" << entropy << "\t" << entropy0-entropy  <<  "\t";
//			if(entropy - entropy0 > 1.0e-10){
			if(0){

				for(int i =0;i<N3d;++i) Ftmp[i]=std::real(F[i]);
				efix.FixSol(N0*N0*N0,Ftmp,entropy0/(std::pow(2*R,3)));
				for(int i =0;i<N3d;++i) F[i]=c_type{Ftmp[i]};

				sum = 0.;
				for(int i =0; i < N3d; ++i) sum += std::real(F[i]);
				sum /= N3d;

				entropy = calEntropy();
				fft3();
				std::cout << entropy << "\t" << entropy0 - entropy << "\t" << sum;
			}
			std::cout << std::endl;
		}
		else
			std::cout <<std::fixed <<std::setprecision(8)<< t << "\t" << sum<< "\t" << entropy << "\t-\n";
		entropy0 = entropy; step ++;
	}while(step < 8);


	std::cout << "# Start the evolution using DUFORT...\n";
	std::cout << "# t\tentropy_{n}\tentropy_{n-1}-entropy_{n}\n";
	step = 0;
	do{
		t += dt;
		qcal.calQ(Q,Fhat);
		for(int l = 0; l < N3d; ++l) F1[l] = Fhat[l];
		for(int l = 0; l < N3d; ++l)
			Fhat[l] = 1./(1.+beta*dt) *(F1[l]+dt*Q[l]+beta*dt*(F1[l]-F_old[l]));

		/// if calculate entropy
		fft3(isInverse);

		r_type sum = 0.;
		for(int i =0; i < N3d; ++i) sum += std::real(F[i]);
		sum /= N3d;
		entropy = calEntropy();
		if(step){
			std::cout <<std::fixed <<std::setprecision(8)<< t  << "\t" << sum<< "\t" << entropy << "\t" << entropy0-entropy  <<  "\t";
//			if(entropy - entropy0 > 1.0e-10){
			if(0){

				for(int i =0;i<N3d;++i) Ftmp[i]=std::real(F[i]);
				efix.FixSol(N0*N0*N0,Ftmp,entropy0/(std::pow(2*R,3)));
				for(int i =0;i<N3d;++i) F[i]=c_type{Ftmp[i]};

				sum = 0.;
				for(int i =0; i < N3d; ++i) sum += std::real(F[i]);
				sum /= N3d;

				entropy = calEntropy();
				fft3();
				std::cout << entropy << "\t" << entropy0 - entropy << "\t" << sum;
			}
			std::cout << std::endl;
		}
		else
			std::cout <<std::fixed <<std::setprecision(8)<< t << "\t" << sum<< "\t" << entropy << "\t-\n";
		entropy0 = entropy; step ++;
		
		for(int i = 0; i< N3d; ++i) F_old[i] = F1[i];
	}while(t < TEnd);

}


#undef TEMPLATE
#undef THIS



#endif ///BKW_H__

