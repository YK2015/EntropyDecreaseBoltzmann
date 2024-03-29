#ifndef _CALBBOX_H_
#define _CALBBOX_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <complex>

#define DEBUG 1

template <class c_type, class r_type>
class QCalculator
{
	public:
		int N;	///	number of dof in  each direction, only support odd N now
		int n;	///	number of positive exponents 

		std::vector<c_type> Q;	/// stored in the standard way, complex
		std::vector<c_type> f;	/// stored in fftw way, complex
		std::vector<r_type> B;	/// stored in a specified manner for fast evaluation of Q, real type


		int alpha;	/// the type of the kernel B; alpha = 0 associated with Maxwell molecules
								/// alpha = 1 correpsonds to the hard spheres

public:
		QCalculator(){};
		~QCalculator(){};
		QCalculator(int _N, int _alpha) : alpha(_alpha) {initialize(_N); };

public:
		void initialize(int _N) {
			N = _N; n = N/2;

			int Nf = std::pow(2*N,3);
			int NB = std::pow(N,3);
			f.resize(Nf);	
			Q.resize(NB);
			B.resize(NB*NB);
		};	

		void genB();	/// generate B
		r_type funF(r_type, r_type);
	 	virtual	void calQ(std::vector<c_type> & Qout, 
				const std::vector<c_type> & fin);
		void transf(const std::vector<c_type> & fin); /// input f
		void convertQ(std::vector<c_type> & Qout); ///  output Q
};

#define TEMPLATE template <class c_type, class r_type> 
#define THIS QCalculator<c_type, r_type>

TEMPLATE
r_type THIS::funF(r_type xi, r_type eta)
{
		r_type p = xi + eta, q = xi - eta;
		r_type res = r_type(0);
		if(alpha == 0){
			if(eta < 1.0e-10){
				if(xi < 1.0e-10) res = r_type(1./3.);
				else res = -(xi*cos(xi)+sin(xi))/(xi*xi*xi);
			}
			else{
				if(xi < 1.0e-10) res = -(eta*cos(eta)+sin(eta))/(eta*eta*eta);
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
				else if (std::abs(xi-eta) < 1.0e-10) res = -(-1.-2.*xi*xi+cos(2.*xi)+2.*xi*sin(2.*xi))/(xi*xi*xi*xi);
				else	res = (q*sin(q)+cos(q)-p*sin(q)-cos(p))/(2.*xi*eta*q*q) - 2./(p*p*q*q);
			}
		}
		else res = r_type(0);
		return res;
}

TEMPLATE
void THIS::transf(const std::vector<c_type> & fin)
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
void THIS::convertQ(std::vector<c_type> & Qout)
{
  int k = 0, kout1, kout2, kout3, kout;
	for(int k1 = 0; k1 < N; ++k1){
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

}

TEMPLATE
void THIS::genB()
{
	for(auto& x : B) x = 0.;
	r_type xi, eta;
	r_type lambda = r_type(2./(3.+sqrt(2)) * 4.*atan(1.));

	int m1,m2,m3,m, ib =0;
	int xi1,xi2,xi3, eta1,eta2,eta3;
	int ibm = 0;

	std::vector<r_type> Bmm(N*N*N);
	for(int k1 = 0; k1 < N; ++k1){
		xi1 = k1+k1-2*n;  xi1 *= xi1;
		for(int k2 = 0; k2 < N; ++k2){
			xi2 = k2+k2-2*n; xi2 *= xi2;
			for(int k3 = 0; k3 < N; ++k3){
				xi3 = k3+k3-2*n; xi3 *= xi3;
				xi = r_type(sqrt(xi1+xi2+xi3)*lambda);	
				eta = r_type(0);
				Bmm[ibm++] = funF(xi,eta);
			}
		}
	}



	for(int l1 = 0; l1 < N; ++l1)
		for(int l2 = 0; l2 < N; ++l2)
			for(int l3 = 0; l3 < N; ++l3)
				for(int k1 = 0; k1 < N; ++k1){
					m1 = n-l1+k1;
					if(m1 < 0) m1 += N;
					else if (m1 > N-1) m1 -= N;
					xi1 = l1-n+m1-n; xi1 *= xi1; eta1 = l1-m1; eta1 *= eta1;

					for(int k2 = 0; k2 < N; ++k2){
						m2 = n-l2+k2;
						if(m2 < 0) m2 += N;
						else if (m2 > N-1) m2 -= N;
						xi2 = l2-n+m2-n; xi2 *= xi2; eta2 = l2-m2; eta2 *= eta2;

						for(int k3 = 0; k3 < N; ++k3){
							m3 = n-l3+k3;
							if(m3 < 0) m3 += N;
							else if (m3 > N-1) m3 -= N;
							m = (m1*N+m2)*N+m3;
							xi3 = l3-n+m3-n; xi3 *= xi3; eta3 = l3-m3; eta3 *= eta3;
							xi = r_type(sqrt(xi1+xi2+xi3)) * lambda;	
							eta = r_type(sqrt(eta1+eta2+eta3)) * lambda;	
							B[ib++] = funF(xi,eta) - Bmm[m];
						}
					}
				}
				
}

TEMPLATE
void THIS::calQ(std::vector<c_type > & Qout,
		const std::vector<c_type> & fin)
{
	for(auto& x : Q) x = c_type(0);
	transf(fin);

	int Nf = 2*N;
	int m = 0,k,l,k_l;
	for(int l1 = n+1; l1 < N+n+1; ++l1){
		for(int l2 = n+1; l2 < N+n+1; ++l2){
			for(int l3 = n+1; l3 < N+n+1; ++l3){
				k = 0; l = (l1*Nf+l2)*Nf + l3;
				for(int k1 = 3*n+2-l1; k1 < 5*n+3-l1; ++k1){
					for(int k2 = 3*n+2-l2; k2 < 5*n+3-l2; ++k2){
						for(int k3 = 3*n+2-l3; k3 < 5*n+3-l3; ++k3){
							k_l = (k1*Nf+k2)*Nf + k3;
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

template<class c_type, class r_type>
class directEvaluation
{
	public:
		int N, n;
		std::vector<c_type> f; 
		std::vector<c_type> Q;
		std::vector<r_type> B;

		int alpha;	/// the type of the kernel B; alpha = 0 associated with Maxwell molecules
								/// alpha = 1 correpsonds to the hard spheres

	public:
		directEvaluation(){};
		~directEvaluation(){};
		directEvaluation(int _N, int _alpha) : alpha(_alpha) {initialize(_N); };

public:
		void initialize(int _N) {
			N = _N; n = N/2;
			int Nf = 8*N*N*N;
			int NQ = N*N*N;
			int NB = NQ*NQ;
			f.resize(Nf);	
			Q.resize(NQ);
			B.resize(NB);
			std::cout << "#\tN = " << N << "\n#\tsize(f) = " << Nf 
				<< "\n#\tsize(Q) = " << NQ
				<< "\n#\tsize(B) = " << NB << "\n";
		};	

		r_type funF(r_type, r_type);

	  void genB();	/// generate B
	 	virtual	void calQ(std::vector<c_type> & Qout, 
				const std::vector<c_type> & fin);
		void transf(const std::vector<c_type> & fin);
		void convertQ(std::vector<c_type> & Qout); ///  output Q
};

#define TEMPLATE template <class c_type, class r_type>
#define THIS directEvaluation<c_type, r_type> 
TEMPLATE
r_type THIS::funF(r_type xi, r_type eta)
{
		r_type p = xi + eta, q = xi - eta;
		r_type res = r_type(0);
		if(alpha == 0){
			if(eta < 1.0e-10){
				if(xi < 1.0e-10) res = r_type(1./3.);
				else res = -(xi*cos(xi)+sin(xi))/(xi*xi*xi);
			}
			else{
				if(xi < 1.0e-10) res = -(eta*cos(eta)+sin(eta))/(eta*eta*eta);
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
				else if (std::abs(xi-eta) < 1.0e-10) res = -(-1.-2.*xi*xi+cos(2.*xi)+2.*xi*sin(2.*xi))/(xi*xi*xi*xi);
				else	res = (q*sin(q)+cos(q)-p*sin(q)-cos(p))/(2.*xi*eta*q*q) - 2./(p*p*q*q);
			}
		}
		else res = r_type(0);
		return res;
}

TEMPLATE
void THIS::genB()
{
	int NB = B.size();
	for(int i = 0; i < NB; ++i)
			B[i] = 0.;

	r_type xi, eta;
	r_type lambda = r_type(2./(3.+sqrt(2)) * 4.*atan(1.));
	int xi1,xi2,xi3,eta1,eta2,eta3, ib=0;

	std::vector<r_type> Bmm(N*N*N);
	for(int k1 = 0; k1 < N; ++k1){
		xi1 = k1+k1-2*n;  xi1 *= xi1;
		for(int k2 = 0; k2 < N; ++k2){
			xi2 = k2+k2-2*n; xi2 *= xi2;
			for(int k3 = 0; k3 < N; ++k3){
				xi3 = k3+k3-2*n; xi3 *= xi3;
				xi = r_type(sqrt(xi1+xi2+xi3)*lambda);	
				eta = r_type(0);
				Bmm[ib++] = funF(xi,eta);
			}
		}
	}

	ib = 0; int ibm = 0;
	for(int l1 = 0; l1 < N; ++l1)
		for(int l2 = 0; l2 < N; ++l2)
			for(int l3 = 0; l3 < N; ++l3){
				ibm = 0;	
				for(int m1 = 0; m1 < N; ++m1){
					xi1 = l1-n+m1-n; xi1 *= xi1; eta1 = l1-m1; eta1 *= eta1;
					for(int m2 = 0; m2 < N; ++m2){
						xi2 = l2-n+m2-n; xi2 *= xi2; eta2 = l2-m2; eta2 *= eta2;
						for(int m3 = 0; m3 < N; ++m3){
							xi3 = l3-n+m3-n; xi3 *= xi3; eta3 = l3-m3; eta3 *= eta3;
							xi = r_type(sqrt(xi1+xi2+xi3)) * lambda;	
							eta = r_type(sqrt(eta1+eta2+eta3)) * lambda;	
							B[ib++] = funF(xi,eta) - Bmm[ibm++];
						}
					}
				}
			}
		
	
}

TEMPLATE
void THIS::calQ(std::vector<c_type> & Qout, 
				const std::vector<c_type> & fin)
{
	for(auto & x : Q) x = 0.;

	int Nb = N*N*N;

	transf(fin);
	int m1, m2, m3, k,m,l;
	for(int l1 = 0; l1 < N; ++l1){
		for(int l2 = 0; l2 < N; ++l2){
			for(int l3 = 0; l3 < N; ++l3){
				l = (l1*N+l2)*N+l3; k = 0;

				for(int k1 = 0; k1 < N; ++k1){
					if(l1 < k1-n) m1 = k1-l1-n-1;
					else if (l1 > k1+n) m1 = k1-l1+3*n+1;
					else m1 = k1-l1+n;

					for(int k2 = 0; k2 < N; ++k2){
						if(l2 < k2-n) m2 = k2-l2-n-1;
						else if (l2 > k2+n) m2 = k2-l2+3*n+1;
						else m2 = k2-l2+n;

						for(int k3 = 0; k3 < N; ++k3){
							if(l3 < k3-n) m3 = k3-l3-n-1;
							else if (l3 > k3+n) m3 = k3-l3+3*n+1;
							else m3 = k3-l3+n;
							m = (m1*N+m2)*N+m3;
							Q[k++] += f[l]*f[m]*B[l*Nb+m];
						}
					}
				}
			}
		}
	}

	convertQ(Qout);

}

TEMPLATE
void THIS::transf(const std::vector<c_type> & fin)
{
	int k = 0, kin1, kin2, kin3, kin;
	for(int k1 = 0; k1 < N; ++k1){
		kin1 = (k1 < n) ? k1+n+1 : k1-n;
		for(int k2 = 0; k2 < N; ++k2){
			kin2 = (k2 < n) ? k2+n+1 : k2-n;
			for(int k3 = 0; k3 < N; ++k3){
				kin3 = (k3 < n) ? k3+n+1 : k3-n;
				kin = (kin1*N+kin2)*N + kin3;
				f[k++] = fin[kin];	
			}
		}
	}
}
TEMPLATE
void THIS::convertQ(std::vector<c_type> & Qout)
{
  int k = 0, kout1, kout2, kout3, kout;
	for(int k1 = 0; k1 < N; ++k1){
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

}
#undef TEMPLATE
#undef THIS


#undef DEBUG

#endif
