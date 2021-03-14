#ifndef _CALBBOX1D_H_
#define _CALBBOX1D_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <assert.h>


template <class value_type>
class QCalculator
{
	public:
		int N;	///	number of dof in  each direction, only support odd N now
		int n;	///	number of positive exponents 

		std::vector<value_type> Q;	/// stored in the standard way 
		std::vector<value_type> f;	/// stored in fftw way
		std::vector<value_type> B;	/// stored in a specified manner for fast evaluation of Q


		int alpha;	/// the type of the kernel B; alpha = 0 associated with Maxwell molecules
								/// alpha = 1 correpsonds to the hard spheres

public:
		QCalculator(){};
		~QCalculator(){};
		QCalculator(int _N, int _alpha) : alpha(_alpha) {initialize(_N); };

public:
		void initialize(int _N) {
			N = _N; n = N/2;

			int Nf = 2*N;
			int NB = N;
			f.resize(Nf);	
			Q.resize(NB);
			B.resize(NB*NB);
		};	

		value_type funF(value_type, value_type);

	  void genB();	/// generate B
	 	void calQ(std::vector<value_type> & Qout, 
				const std::vector<value_type> & fin);
		void transf(const std::vector<value_type> & fin);
	
};

#define TEMPLATE template <class value_type>
#define THIS QCalculator<value_type> 
TEMPLATE
value_type THIS::funF(value_type xi, value_type eta)
{
		value_type p = xi + eta, q = xi - eta;
		value_type res = value_type(0);
		if(alpha == 0){
			if(eta < 1.0e-10){
				if(xi < 1.0e-10) res = value_type(1./3.);
				else res = -(xi*cos(xi)+sin(xi))/(xi*xi*xi);
			}
			else{
				if(xi < 1.0e-10) res = -(eta*cos(eta)+sin(eta))/(eta*eta*eta);
				else if (abs(xi-eta) < 1.0e-10) res = (xi-cos(xi)*sin(xi))/(2.*xi*xi*xi);
				else	res = (p*sin(q) - q*sin(p))/(2.*xi*eta*p*q);
			}
				
		}
		else if (alpha == 1){
			if(eta < 1.0e-10){
				if(xi < 1.0e-10) res = value_type(1./4.);
				else res = (-2.-(-2.+xi*xi)*cos(xi)+2.*xi*sin(xi))/(xi*xi*xi*xi);
			}
			else{
				if(xi < 1.0e-10) res = (-2.-(-2.+eta*eta)*cos(eta)+2.*eta*sin(eta))/(eta*eta*eta*eta);
				else if (abs(xi-eta) < 1.0e-10) res = -(-1.-2.*xi*xi+cos(2.*xi)+2.*xi*sin(2.*xi))/(xi*xi*xi*xi);
				else	res = (q*sin(q)+cos(q)-p*sin(q)-cos(p))/(2.*xi*eta*q*q) - 2./(p*p*q*q);
			}
		}
		else res = value_type(0);
		return res;
}

TEMPLATE 
void THIS::transf(const std::vector<value_type> & fin)
{
	for(int i = 0; i < N; ++i){
		f[i] = fin[i];
		f[i+N] = fin[i];
	}
}

TEMPLATE
void THIS::genB()
{
	for(auto& x : B) x = 0.;
	value_type xi, eta;
	value_type lambda = value_type(2./(3.+sqrt(2)) * 4.*atan(1.));

	int m,ib=0;
	for(int l = 0; l < N; ++l)
		for(int k = 0; k < N; ++k){
				m = n-l+k;
				if(m < 0) m += N;
				else if (m > N-1) m -= N;

				xi = fabs(l-n+m-n) * lambda;	/// note l corresponds to l -n
				eta = fabs(l-m) * lambda;
				B[ib] = funF(xi,eta);
				xi = fabs(m-n+m-n) * lambda;
				eta = 0;
				B[ib++] -= funF(xi,eta);
			//	B[ib++] = m;
		}
}

TEMPLATE
void THIS::calQ(std::vector<value_type > & Qout,
		const std::vector<value_type> & fin)
{
	for(auto& x : Qout) x = 0;
	transf(fin);

	int m = 0;
	std::cout << "n = " << n << "\n";
	int iQ;
	for(int l = n+1; l < N+n+1; ++l){
		iQ = 0;
		for(int k = 3*n+2-l; k < 5*n+3-l; ++k){
			Qout[iQ++] += f[l] * f[k] * B[m++];
		}

	}

}

#undef TEMPLATE
#undef THIS


template<class value_type>
class directEvaluation
{
	public:
		int N, n;
		std::vector<value_type> f; 
		std::vector<value_type> Q;
		std::vector<std::vector<value_type> > B;

		int alpha;	/// the type of the kernel B; alpha = 0 associated with Maxwell molecules
								/// alpha = 1 correpsonds to the hard spheres

	public:
		directEvaluation(){};
		~directEvaluation(){};
		directEvaluation(int _N, int _alpha) : alpha(_alpha) {initialize(_N); };

public:
		void initialize(int _N) {
			N = _N; n = N/2;
			f.resize(N);	
			Q.resize(N);
			B.resize(N,std::vector<value_type>(N));
		};	

		value_type funF(value_type, value_type);

	  void genB();	/// generate B
	 	virtual	void calQ(std::vector<value_type> & Qout, 
				const std::vector<value_type> & fin);
		void transf(const std::vector<value_type> & fin);
};

#define TEMPLATE template <class value_type>
#define THIS directEvaluation<value_type> 
TEMPLATE
value_type THIS::funF(value_type xi, value_type eta)
{
		value_type p = xi + eta, q = xi - eta;
		value_type res = value_type(0);
		if(alpha == 0){
			if(eta < 1.0e-10){
				if(xi < 1.0e-10) res = value_type(1./3.);
				else res = -(xi*cos(xi)+sin(xi))/(xi*xi*xi);
			}
			else{
				if(xi < 1.0e-10) res = -(eta*cos(eta)+sin(eta))/(eta*eta*eta);
				else if (abs(xi-eta) < 1.0e-10) res = (xi-cos(xi)*sin(xi))/(2.*xi*xi*xi);
				else	res = (p*sin(q) - q*sin(p))/(2.*xi*eta*p*q);
			}
				
		}
		else if (alpha == 1){
			if(eta < 1.0e-10){
				if(xi < 1.0e-10) res = value_type(1./4.);
				else res = (-2.-(-2.+xi*xi)*cos(xi)+2.*xi*sin(xi))/(xi*xi*xi*xi);
			}
			else{
				if(xi < 1.0e-10) res = (-2.-(-2.+eta*eta)*cos(eta)+2.*eta*sin(eta))/(eta*eta*eta*eta);
				else if (abs(xi-eta) < 1.0e-10) res = -(-1.-2.*xi*xi+cos(2.*xi)+2.*xi*sin(2.*xi))/(xi*xi*xi*xi);
				else	res = (q*sin(q)+cos(q)-p*sin(q)-cos(p))/(2.*xi*eta*q*q) - 2./(p*p*q*q);
			}
		}
		else res = value_type(0);
		return res;
}

TEMPLATE
void THIS::genB()
{
	for(int i = 0; i < N; ++i)
		for(int j = 0; j < N; ++j) 
			B[i][j] = 0.;
	value_type xi, eta;
	value_type lambda = value_type(2./(3.+sqrt(2)) * 4.*atan(1.));
	for(int l = 0; l < N; ++l){
		for(int m = 0; m < N; ++m){
			xi = fabs(l-n+m-n) * lambda;	/// note l corresponds to l -n
			eta = fabs(l-m) * lambda;
			B[l][m] = funF(xi,eta);
			xi = fabs(m-n+m-n) * lambda;
			eta = 0;
			B[l][m] -= funF(xi,eta);
		}
	}
}

TEMPLATE
void THIS::calQ(std::vector<value_type> & Qout, 
				const std::vector<value_type> & fin)
{
	for(auto & x : Qout) x = 0.;

	transf(fin);
	for(int l = 0; l < N; ++l){
		for(int k = 0; k < N; ++k){
			int m;
			if(l < k-n) m = k-l-n-1;
			else if (l > k+n) m = k-l+3*n+1;
			else m = k-l+n;
			Qout[k] += f[l]*f[m]*B[l][m];
		}
	}

}

TEMPLATE
void THIS::transf(const std::vector<value_type> & fin)
{
	for(int i = 0; i < n; ++i)
		f[i] = fin[i+n+1];
	for(int i = 0; i < n+1; ++i)
		f[i+n] = fin[i];
}

#undef TEMPLATE
#undef THIS

#endif
