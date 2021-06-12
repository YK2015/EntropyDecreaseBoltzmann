#ifndef _ENTROPYFIX_H_
#define _ENTROPYFIX_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <assert.h>
//#include <complex>
//#include "mkl.h"


template <class i_type, class r_type>
class EntropyFix
{
	private:
		i_type Nt;						///	total number of dof in all directions
		r_type eps;

		r_type Averagef = 1.0;			/// constant state of f
		r_type *LogInput, *LogOutput;	/// storage of function to take log function in MKL
		bool   Average_flag = false;	/// indicator of the average of f is set or not

	private:
		EntropyFix(){	assert(("Default constructor is not allowed in EntropyFix class", 0 == 1));	};

	public:
		~EntropyFix(){ delete [] LogInput;	delete [] LogOutput; };
		EntropyFix(i_type _Nt, r_type _eps = 1.0e-10) : Nt(_Nt), eps(_eps) {
			LogInput = new r_type[Nt];
			LogOutput = new r_type[Nt];
		}

	public:

		r_type EntropyCal(const std::vector<r_type> & fin);			/// Calculate the entropy for given fin
		r_type EntropyCal(i_type Nf, r_type * fin);			/// Calculate the entropy for given fin

		void setAverage(const std::vector<r_type> & fin);			/// Calculate the average
		void setAverage(i_type Nf, r_type* fin);			/// Calculate the average

		void FixSol(std::vector<r_type> & fin, r_type AimEntropy);	/// Bisection method to find the entropic solution
		void FixSol(i_type Nf, r_type* fin, r_type AimEntropy);	/// Bisection method to find the entropic solution

		void FixSol(const std::vector<r_type> &fin, std::vector<r_type> &fout, r_type AimEntropy);	/// If we need the input function unchanged. Future work.
};

#define TEMPLATE template <class i_type, class r_type> 
#define THIS EntropyFix<i_type, r_type>

	TEMPLATE
void THIS::setAverage(const std::vector<r_type> & fin)
{
	i_type Nf = fin.size();
	assert(Nf == Nt);

	Average_flag = true;
	//Averagef = cblas_dasum(Nf, &(fin[0]), 1);
	Averagef = 0.;
	for(int i = 0; i < Nf; ++i) Averagef += fin[i];
	Averagef /= (r_type)Nf;
	std::cout << "\nAverafe f\t" << Averagef << "\n";
}

	TEMPLATE
void THIS::setAverage(i_type Nf, r_type* fin)
{
	assert(Nf == Nt);

	Average_flag = true;
//	Averagef = cblas_dasum(Nf, &(fin[0]), 1);
	Averagef = 0.;
	for(int i = 0; i < Nf; ++i) Averagef += fin[i];

	Averagef /= (r_type)Nf;
	std::cout << "\nAverafe f\t" << Averagef << "\n";
//	std::cout << "\n\t\t#Before FIxing:: " << Averagef << "\n";
}

	TEMPLATE
r_type THIS::EntropyCal(const std::vector<r_type> & fin)
{
	i_type Nf = fin.size();
	assert(Nf == Nt);

	r_type res = 0.0;

	if (!Average_flag) {
		setAverage(fin);
	}

	/// I will write two versions. The following one is standard one, meaning no MKL.
	/// The standard one is faster, compared to the MKL_sequential.
#ifndef _MKL_H
	for (i_type i=0; i<Nf; i++) {
		if (fin[i] > eps) {
			res += fin[i] * log(fin[i]);
		}
	}
#elif
	
	/// Second version with MKL, r_type should be double. Otherwise use v?Ln

	//for (i_type i=0; i<Nf; i++) {
	i_type Nnonzero=0;	/// number of non-zero components of fin
	for (auto x : fin)	{
		if (x > eps) {
			LogInput[Nnonzero++] = x; 
		}
		/*	else
				{
				LogInput[i++] = 1.0;
				}	*/
	}
	vdLn( Nnonzero, (const r_type*)LogInput, &(LogOutput[0]) );
	res = cblas_ddot(Nnonzero, &(LogInput[0]), 1, &(LogOutput[0]), 1);

#endif

	res = res / (r_type)Nf - Averagef;
	return res;
}

	TEMPLATE
r_type THIS::EntropyCal(i_type Nf, r_type *fin)
{
	assert(Nf == Nt);

	r_type res = 0.0;

	if (!Average_flag) {
		setAverage(Nf, fin);
	}

	/// I will write two versions. The following one is standard one, meaning no MKL.
	/// The standard one is faster, compared to the MKL_sequential.
	
#ifndef _MKL_H
	for (i_type i=0; i<Nf; i++) {
		if (fin[i] > eps) {
			res += fin[i] * log(fin[i]);
		}
	}
#elif

	/// Second version with MKL, r_type should be double. Otherwise use v?Ln

	//for (i_type i=0; i<Nf; i++) {
	i_type Nnonzero=0;	/// number of non-zero components of fin
	for(i_type i=0; i<Nf; i++){
		if(fin[i]>eps)
			LogInput[Nnonzero++] = fin[i]; 
	}

	vdLn( Nnonzero, (const r_type*)LogInput, &(LogOutput[0]) );
	res = cblas_ddot(Nnonzero, &(LogInput[0]), 1, &(LogOutput[0]), 1);
#endif

	res = res / (r_type)Nf - Averagef;
	return res;
}


	TEMPLATE
void THIS::FixSol(std::vector<r_type> & fin, r_type AimEntropy)
{
	i_type Nf = fin.size();
	assert(Nf == Nt);

	std::vector<r_type> fMid(Nt);

	r_type betaL = 0.0, betaR = 1.0, betaMid, entMid, betaFbar;
	r_type BiEps = eps;					///	this epsilon could be different to the epsilon in calculating the entropy

	while ((betaR - betaL) > BiEps) {
		betaMid = 0.5 * (betaL + betaR);

		betaFbar = betaMid * Averagef;
		for (auto &x : fMid) {
			x = betaFbar;
		}

//		cblas_daxpy(Nf, 1.0 - betaMid, &(fin[0]), 1, &(fMid[0]), 1);
		for(int i = 0; i< Nf; ++i) fMid[i] += (1.-betaMid)*fin[i];

		entMid = EntropyCal(fMid);
		if ( (entMid < AimEntropy) && ((AimEntropy - entMid)<eps) ) {
			//cblas_dcopy(Nf, &(fMid[0]), 1, &(fin[0]), 1);
			for(int i =0;i < Nf;++i) fin[i] = fMid[i];
			return;
		}

		if (entMid < AimEntropy) {
			betaR = betaMid;
		}
		else
		{
			betaL = betaMid;
		}
	}

	betaFbar = betaR * Averagef;
	//cblas_dscal(Nf, 1.0 - betaR, &(fin[0]), 1);
	for(auto &x :fin) x *= (1.0-betaR);
	for (auto &x : fin) {
		x += betaFbar;
	}
}


TEMPLATE
void THIS::FixSol(i_type Nf, r_type * fin, r_type AimEntropy)
{
	assert(Nf == Nt);

	//std::vector<r_type> fMid(Nt);
	r_type fMid[Nt];

	r_type betaL = 0.0, betaR = 1.0, betaMid, entMid, betaFbar;
	r_type BiEps = eps;					///	this epsilon could be different to the epsilon in calculating the entropy

	while ((betaR - betaL) > BiEps) {
		betaMid = 0.5 * (betaL + betaR);

		if(!Average_flag) setAverage(Nf, 	fin);
		betaFbar = betaMid * Averagef;
		for(i_type i=0; i<Nt; i++)
			fMid[i] = betaFbar;

		//cblas_daxpy(Nf, 1.0 - betaMid, fin, 1, fMid, 1);
		for(int i = 0; i< Nf; ++i) fMid[i] += (1.-betaMid)*fin[i];
		entMid = EntropyCal(Nt,fMid);
		if ( (entMid < AimEntropy) && ((AimEntropy - entMid)<eps) ) {
			//cblas_dcopy(Nf, fMid, 1, fin, 1);
			for(int i =0;i < Nf;++i) fin[i] = fMid[i];
			return;
		}

		if (entMid < AimEntropy) {
			betaR = betaMid;
		}
		else
		{
			betaL = betaMid;
		}
	}

	betaFbar = betaR * Averagef;
//	cblas_dscal(Nf, 1.0 - betaR,fin, 1);
	for(i_type i=0; i < Nf; ++i) fin[i] *= (1.0-betaR);
	for(i_type i=0; i<Nf; i++)
		fin[i] += betaFbar;
}
#undef TEMPLATE
#undef THIS


#endif

