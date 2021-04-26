#ifndef _QCALCULATOR_H_
#define _QCALCULATOR_H_

#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <complex>

//const int N = 11;
//const int NQ =  N*N*N;
//const int Nf = 8*NQ;
//const int NB = NQ*NQ;

template <class c_type, class r_type>
class QCalculator
{
	private:
		int N;
		int NB;
		int NQ;
		int Nf;
		int n;	///	number of positive exponents 

		c_type * Q;	/// stored in the standard way, complex
		c_type * f;	/// stored in fftw way, complex*2
		r_type * B;	/// stored in a specified manner for fast evaluation of Q, real type

		int alpha;	/// the type of the kernel B; alpha = 0 associated with Maxwell molecules
								/// alpha = 1 correpsonds to the hard spheres
								
		bool IS_FILTER = 1;	/// adopt filtering (default setting)

public:
		QCalculator(){};
		~QCalculator(){delete[] Q;delete[] f; delete[] B;};
		QCalculator(int _N, int _alpha) : N(_N), alpha(_alpha) {initialize();};

public:
		
		void initialize() {
			n = (N-1)/2;

			Nf = std::pow(2*N,3);
			NQ = std::pow(N,3);
			NB = NQ*NQ;
			f = new c_type[Nf];	
			Q = new c_type[NQ];
			B = new r_type[NQ*NQ];
		};	
		
		void closeFilter(){IS_FILTER = 0;};

		void genB();	/// generate B
		r_type funF(r_type, r_type);
		r_type filtering(int);
	 	virtual	void calQ(c_type * Qout, c_type* fin);
		void transf(c_type * fin); /// input f
		void convertQ(c_type * Qout); ///  output Q
};

#define TEMPLATE template <class c_type, class r_type> 
#define THIS QCalculator<c_type, r_type>

TEMPLATE
r_type THIS::filtering(int beta)
{
	r_type tmp = r_type(n+1);
	r_type PI = r_type(4.*atan(1.));
	r_type b = fabs(r_type(beta));
	return ((tmp-b)*cos(PI*b/tmp)+sin(PI*b/tmp)/tan(PI/tmp))/tmp;
}

TEMPLATE
r_type THIS::funF(r_type xi, r_type eta)
{
		r_type p = xi + eta, q = xi - eta;
		r_type res = r_type(0);
		if(alpha == 0){
			if(eta < 1.0e-10){
				if(xi < 1.0e-10) res = r_type(1./3.);
				else res = (-xi*cos(xi)+sin(xi))/(xi*xi*xi);
			}
			else{
				if(xi < 1.0e-10) res = (-eta*cos(eta)+sin(eta))/(eta*eta*eta);
				else if (std::abs(xi-eta) < 1.0e-10) res = (xi-cos(xi)*sin(xi))/(2.*xi*xi*xi);
				else	res = (p*sin(q) - q*sin(p))/(2.*xi*eta*p*q);
			}
				
		}
		else if (alpha == 1){
			if(eta < 1.0e-10){
				if(xi < 1.0e-10) res = r_type(1./4.);
				else res = (-2.-(-2.+xi*xi)*cos(xi)+2.*xi*sin(xi))/(xi*xi*xi*xi);
			}
			else{
				if(xi < 1.0e-10) res = (-2.-(-2.+eta*eta)*cos(eta)+2.*eta*sin(eta))/(eta*eta*eta*eta);
				else if (std::abs(xi-eta) < 1.0e-10) res = -(-1.-2.*xi*xi+cos(2.*xi)+2.*xi*sin(2.*xi))/(8.*xi*xi*xi*xi);
				else	res = (q*sin(q)+cos(q))/(2.*xi*eta*q*q)-(p*sin(p)+cos(p))/(2.*xi*eta*p*p) - 2./(p*p*q*q);
			}
		}
		else res = r_type(0);
		return res;
}

TEMPLATE
void THIS::transf(c_type * fin)
{
	int k = 0, kin1, kin2, kin3, kin;
	for(int k1 = 0; k1 < 2*N; ++k1){
		kin1 = (k1 < N) ? k1 : k1-N;
		for(int k2 = 0; k2 < 2*N; ++k2){
			kin2 = (k2 < N) ? k2 : k2-N;
			for(int k3 = 0; k3 < 2*N; ++k3){
				kin3 = (k3 < N) ? k3 : k3-N;
				kin = (kin1*N+kin2)*N + kin3;
				f[k++] = fin[kin];					
			}
		}
	}
}

TEMPLATE
void THIS::convertQ(c_type* Qout)
{
	/*
  int k = 0, kout1, kout2, kout3, kout;
	
	for(int k1 = 0; k1 < N; ++k1){				/// is this a bug? number of elements in Q should be 2N: No, since Q is of size N
		kout1 = (k1 < n) ? k1+n+1 : k1-n;
		for(int k2 = 0; k2 < N; ++k2){
			kout2 = (k2 < n) ? k2+n+1 : k2-n;
			for(int k3 = 0; k3 < N; ++k3){
				kout3 = (k3 < n) ? k3+n+1 : k3-n;
				kout = (kout1*N+kout2)*N + kout3;
				Qout[kout] = Q[k++];					
			}
		}
	}
	*/

	int k = 0, kq1, kq2, kq3, kq;
	for(int k1 = 0; k1 < N; ++k1){
		kq1 = (k1 < n+1) ? k1+n : k1-(n+1);
		for(int k2 = 0; k2 < N; ++k2){
			kq2 = (k2 < n+1) ? k2+n : k2-(n+1);
			for(int k3 = 0; k3 < N; ++k3){
				kq3 = (k3 < n+1) ? k3+n : (k3-n-1);
				kq = (kq1*N + kq2)*N + kq3;
				Qout[k++] = Q[kq];
			}
		}
	}
}

TEMPLATE
void THIS::genB()
{
	for(int i = 0; i < NB; ++i) B[i] = 0.;
	r_type xi, eta;
	r_type lambda = r_type(2./(3.+sqrt(2.)) * M_PI);

	auto calXi = [&](int x1, int x2, int x3){
		return r_type(sqrt(x1*x1+x2*x2+x3*x3))*lambda;
	};
	auto calFil = [&](int x1, int x2, int x3){
		return filtering(x1)*filtering(x2)*filtering(x3);
	};

	int m1,m2,m3,m, ib =0;
	int xi1,xi2,xi3, eta1,eta2,eta3;
	int ibm = 0;

	r_type rval = (4.*M_PI) * std::pow(2.*lambda/M_PI * 6.62132/N,3+alpha);

	r_type fil_val = r_type(1.);
	std::vector<r_type> Bmm(N*N*N);
	std::ofstream writeB("diagB.dat");
	for(int k1 = -n; k1 < n+1; ++k1)
		for(int k2 = -n; k2 < n+1; ++k2)
			for(int k3 = -n; k3 < n+1; ++k3){
				xi = calXi(2*k1,2*k2,2*k3);
				//xi = r_type(sqrt(xi1+xi2+xi3))*lambda;	 /// it looks like a bug
				eta = r_type(0);
				if(IS_FILTER) fil_val = calFil(k1,k2,k3);
				Bmm[ibm++] = funF(xi,eta)*fil_val*fil_val*rval;
				writeB << k1 << "\t" << k2 << "\t" << k3 << "\t" << Bmm[ibm-1]<< "\n";


			}
		
	writeB.close();
	

  r_type sigmal = r_type(1);
  r_type sigmam = r_type(1);
	for(int l1 = 0; l1 < N; ++l1)
		for(int l2 = 0; l2 < N; ++l2)
			for(int l3 = 0; l3 < N; ++l3){
				if(IS_FILTER) sigmal = calFil(l1-n,l2-n,l3-n);
				for(int k1 = 0; k1 < N; ++k1){
					m1 = k1-l1+n;	/// x1-index in the array
					if(m1 < 0) m1 += N;
					else if (m1 > N-1) m1 -= N; /// move the outside element back
					xi1 = l1-n+m1-n;  eta1 = l1-m1; /// xi = l+m, eta = l-m

					for(int k2 = 0; k2 < N; ++k2){
						m2 = k2-l2+n;
						if(m2 < 0) m2 += N;
						else if (m2 > N-1) m2 -= N;
						xi2 = l2-n+m2-n; eta2 = l2-m2; 

						for(int k3 = 0; k3 < N; ++k3){
							m3 = k3-l3+n;
							if(m3 < 0) m3 += N;
							else if (m3 > N-1) m3 -= N;
							m = (m1*N+m2)*N+m3;
							xi3 = l3-n+m3-n; eta3 = l3-m3; 
							xi = calXi(xi1,xi2,xi3);
							eta = calXi(eta1,eta2,eta3);
							if(IS_FILTER) sigmam = calFil(m1-n,m2-n,m3-n);
							B[ib++] = funF(xi,eta)*sigmal*sigmam*rval - Bmm[m];
						}
					}
				}
			}
}

TEMPLATE
void THIS::calQ(c_type * Qout, c_type* fin)
{
	for(int i = 0; i < NQ; ++i) Q[i] = c_type(0);
	transf(fin);

	int N2 = 2*N;
	int m = 0,k,l,k_l;
	for(int l1 = n+1; l1 < N+n+1; ++l1){
		for(int l2 = n+1; l2 < N+n+1; ++l2){
			for(int l3 = n+1; l3 < N+n+1; ++l3){
				k = 0; l = (l1*N2+l2)*N2 + l3;
				for(int k1 = 3*n+2-l1; k1 < 5*n+3-l1; ++k1){
					for(int k2 = 3*n+2-l2; k2 < 5*n+3-l2; ++k2){
						for(int k3 = 3*n+2-l3; k3 < 5*n+3-l3; ++k3){
							k_l = (k1*N2+k2)*N2 + k3;
							Q[k++] += f[l] * f[k_l] * B[m++];
						}
					}
				}
			}
		}
	}

	convertQ(Qout);

}
#undef TEMPLATE
#undef THIS


#endif
