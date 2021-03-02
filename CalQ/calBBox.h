#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <assert.h>

#define DEBUG 1

template <class value_type, int DIM>
class QCalculator
{
private:
	int N;	///  number of dof in  each direction, only support odd N now
	int n;

	std::vector<value_type> Q;
	std::vector<value_type> f;
	std::vector<value_type> B;

public:
		QCalculator(){};
		~QCalculator(){};
		QCalculator(int _N) {initialize(_N); };

public:
		void initialize(int _N) {
			N = _N; n = N/2;

			int Nf = std::pow(2*N,DIM);
			int NB = std::pow(N,DIM);
			f.resize(Nf);	
			Q.resize(NB);
			B.resize(NB*NB);
		};	

	  void genB(int alpha);	/// generate B
		void calQ(std::vector<value_type> & Qout, 
				const std::vector<value_type> & fin);

#if DEBUG
		void displayQ()
		{
			int Nb = std::pow(N,DIM);
			std::cout << "# size of Q: " << Nb  << "\n\n";
			for(int i = 0; i < Nb; ++i){
				std::cout << std::fixed << std::setprecision(4) << Q[i] << "\t";
				std::cout << "\n";
			}

		}
		void displayB()
		{
			int Nb = std::pow(N,DIM);
			std::cout << "# size of B: " << Nb*Nb << " = " << Nb << " * " << Nb << "\n\n";
			for(int i = 0; i < Nb; ++i){
				for(int j = 0; j < Nb; ++j){
					std::cout << std::fixed << std::setprecision(4) << B[j*Nb + i] << "\t";
				}
				std::cout << "\n";
			}

		}
#endif
	
};


#define TEMPLATE template <class value_type, int DIM> 
#define THIS QCalculator<value_type,DIM>

TEMPLATE
void THIS::genB(int alpha)
{
	for (auto &x : B) x = 0.;

	/// the expression for the integral 
	auto Func = [&](value_type xi, value_type eta) {
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
	};

	int Nb = std::pow(N,DIM);
	if(DIM == 3){
		int kb=0,xi1,xi2,xi3,eta1,eta2,eta3;
		value_type xi, eta;
		double lambda = 2./(3.+sqrt(2)) * 4.*atan(1.);
		std::vector<value_type> Bmm(Nb);
		for(int k1 = 0; k1 < N; ++k1){
			xi1 = k1+k1-2*n;  xi1 *= xi1;
			for(int k2 = 0; k2 < N; ++k2){
				xi2 = k2+k2-2*n; xi2 *= xi2;
				for(int k3 = 0; k3 < N; ++k3){
					xi3 = k3+k3-2*n; xi3 *= xi3;
					xi = value_type(sqrt(xi1+xi2+xi3)*lambda);	
					eta = value_type(0);
					Bmm[kb++] = Func(xi,eta);
				}
			}
		}

		kb = 0; int kbm = 0;
		for(int k1 = 0; k1 < N; ++k1){
			for(int k2 = 0; k2 < N; ++k2){
				for(int k3 = 0; k3 < N; ++k3){
					kbm = 0;
					for(int l1 = 0; l1 < N; ++l1){
						xi1 = k1+l1-2*n; eta1 = k1-l1;
						xi1 *= xi1; eta1 *= eta1;
						for(int l2 = 0; l2 < N; ++l2){
							xi2 = k2+l2-2*n; eta2 = k2-l2;
							xi2 *= xi2; eta2 *= eta2;
							for(int l3 = 0; l3 < N; ++l3){
								xi3 = k3+l3-2*n; eta3 = k3-l3;
								xi3 *= xi3; eta3 *= eta3;
								xi = value_type(sqrt(xi1+xi2+xi3) * lambda); 
								eta = value_type(sqrt(eta1+eta2+eta3) * lambda); 
								B[kb++] = Func(xi,eta) - Bmm[kbm++];
							}
						}
					}
				}
			}
		}
	}

}


TEMPLATE
void THIS::calQ(std::vector<value_type>& Qout, const std::vector<value_type> & fin)
{
	int factor = std::pow(2,DIM);
	assert(factor*fin.size() == f.size());

	
	if(DIM == 3){
		int kin1, kin2, kin3, kin, k,l,m, m1, m2, m3;
		
		/// assign values of f
		/// f(k1,k2,k3) --index--> (k1*N+k2)*N+k3
		k = 0;
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

		/// calculate Q
		for(auto & x : Q) x = value_type(0);
		int N2 = N*2, mb1, mb2, mb3, mb,lb;
		for(int l1 = 0; l1 < N; ++l1){
			for(int l2 = 0; l2 < N; ++l2){
				for(int l3 = 0; l3 < N; ++l3){
					l = ((l1+n+1)*N2+l2+n+1)*N2 + l3+n+1;
					lb = (l1*N+l2)*N + l3;
					k = 0; 
					for(int k1 = 0; k1 < N; ++k1){
						m1 = k1-l1+N; if(m1 < n+1) mb1 = m1+n; else if (m1 < 3*n+2) mb1 = m1 - n-1; else mb1 = m1-3*n-2;
						for(int k2 = 0; k2 < N; ++k2){
							m2 = k2-l2+N; if(m2 < n+1) mb2 = m2+n; else if (m2 < 3*n+2) mb2 = m2 - n-1; else mb2 = m2-3*n-2;
							for(int k3 = 0; k3 < N; ++k3){
								m3 = k3-l3+N;
								if(m3 < n+1) mb3 = m3+n; else if (m3 < 3*n+2) mb3 = m3 - n-1; else mb3 = m3-3*n-2;
								m = (m1*N2+m2)*N2+m3; mb = (mb1*N+mb2)*N+mb3;
								Q[k++]	 +=  f[l]*f[m]*B[lb*N*N*N+mb];
							}
						}
					}
					
				}
			}
		}

	}

}

#undef THIS
#undef TEMPLATE
#undef DEBUG
