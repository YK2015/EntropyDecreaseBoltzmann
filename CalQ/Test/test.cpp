/**
* Copyright (C) 2021 All rights reserved.
* @file test.cpp
* @brief  one d for symmetric storage
* @author Yang Kuang (matkuan@nus.edu.sg)
* @date 2021-01-26
* @return 
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "TimeWatch.h"
using namespace std;

void refUpdateQ(int N, const vector<double>  &, const vector<double>&,
		vector<double>&);

void UpdateQ(int N, const vector<double > &, const vector<double>&,
		vector<double>&);

void UpdateQLongF(int N, const vector<double > &, const vector<double>&,
		vector<double>&);

double l2error(const vector<double>&, const vector<double>& );

int main(int argc, char * argv[])
{
	int N0 = atoi(argv[1]);
	int N = N0*N0*N0;
	vector<double> f(N);
	vector<double> lf(27*N);
	vector<double> B(N*N);
	for(int i = 0;i < N; ++i){
		f[i] =  i* 0.1; 
		for(int j = 0; j < N; ++j)
			B[i*N+j] = (i+j)* 0.1;
	}

	for(int j = 0; j < 3*N0; ++j){
		int j0 = j % N0;
		for(int i = 0; i < 3*N0; ++i){
			int i0 = i % N0;
			for(int k = 0; k < 3*N0; ++k){
				int k0 = k % N0;
				lf[(j*(3*N0)+i)*3*N0+k]	 = f[(j0*N0+i0)*N0+k0];
			}
		}
	}

	vector<double> Q(N);
	cout << "# Classical method \n";
	refUpdateQ(N0,B,f,Q);

	vector<double> Q2(N);
	cout << "# Update in a new rule \n";
	UpdateQ(N0,B,f,Q2);

	
	vector<double> Q3(N);
	cout << "# Update with Trippled f \n";
	UpdateQLongF(N0,B,lf,Q3);


	auto err = l2error(Q2,Q3);
	cout << "l2 err = " << err << "\n";

	return 0;
}

void refUpdateQ(int N, const vector<double>  & B, const vector<double>& f,
		vector<double>& Q)
{
	TimeWatch<float> tw;
	int N0 = Q.size(); int n = N/2;
	cout << "n = " << n << "\n";
	for(double & x : Q) 	x  = 0.;

#if 1
	int k,l, m, m1, m2,m3;
	for(int k3 = 0; k3 < N; ++k3){
		for(int k2 = 0; k2 < N; ++k2){
			for(int k1 = 0; k1 < N; ++k1){
				k = (k3*N+k2)*N+k1;
				for(int l3 = 0; l3 < N; ++l3){
					if(l3 > k3+n) m3 = k3-l3+3*n+1;
					else if(l3 < k3-n) m3 = k3-l3-n-1;
					else m3 = k3-l3+n;

					for(int l2 = 0; l2 < N; ++l2){
						if(l2 > k2+n) m2 = k2-l2+3*n+1;
						else if(l2 < k2-n) m2 = k2-l2-n-1;
						else m2 = k2-l2+n;

						for(int l1 = 0; l1 < N; ++l1){
							if(l1 > k1+n) m1 = k1-l1+3*n+1;
							else if (l1 < k1-n) m1 = k1-l1-n-1;
							else m1 = k1-l1+n;
									
							l = (l3*N+l2)*N+l1;	m = (m3*N+m2)*N+m1;

							Q[k] += B[k*N0+l] * f[l] * f[m];
						}
					}
				}
			}
		}
	}
#else
	for(int k = 0; k < n; ++k){
			for(int l = 0; l < k+n+1; ++l)
				Q[k] += B[k*N+l] * f[l] * f[k-l+n];
			for(int l = k+n+1; l < N; ++l)
				Q[k] += B[k*N+l] * f[l] * f[k-l+3*n+1];
	}
	for(int k = n; k < N; ++k){
		for(int l = 0; l < k-n; ++l)
			Q[k] += B[k*N+l] * f[l] * f[k-l-n-1];
		for(int l = k-n; l < N; ++l)
			Q[k] += B[k*N+l] * f[l] * f[k-l+n];
	}
#endif
}

void UpdateQ(int N, const vector<double > & B, const vector<double>& f,
		vector<double>& Q)
{
	TimeWatch<float> tw;
	int N0 = Q.size();  int n = N/2;
	cout << "n = " << n << "\n";
	for(double & x : Q) 	x  = 0.;

#if 1
	int k, l, m, m1, m2,m3;
	for(int l3 = 0; l3 < N; ++l3){
	for(int l2 = 0; l2 < N; ++l2){
		for(int l1 = 0; l1 < N; ++l1){
			l = (l3*N+l2)*N+l1;
			for(int k3 = 0; k3 < N; ++k3){
				if(k3 < l3-n) m3 = k3-l3+3*n+1;
				else if(k3 > l3+n) m3 = k3-l3-n-1;
				else m3 = k3-l3+n;
			for(int k2 = 0; k2 < N; ++k2){
				if(k2 < l2-n) m2 = k2-l2+3*n+1;
				else if(k2 > l2+n) m2 = k2-l2-n-1;
				else m2 = k2-l2+n;

				for(int k1 = 0; k1 < N; ++k1){
					if(k1 < l1-n) m1 = k1-l1+3*n+1;
					else if (k1 > l1+n) m1 = k1-l1-n-1;
					else m1 = k1-l1+n;
									
					k = (k3*N+k2)*N + k1; m = (m3*N+m2)*N + m1;
					Q[k] += B[l*N0+k] * f[l] * f[m];
				}
				}
			}
		}
	}
	}
#else
	int m, k1,k2,l1,l2,m1,m2;
	for(int l = 0; l < N0; ++l){
		l2 = l/N; l1 = l%N;
		for(int k = 0; k < N0; ++k){
			k2 = k/N; k1 = k%N;
				if(k2 < l2-n) m2 = k2-l2+3*n+1;
				else if(k2 > l2+n) m2 = k2-l2-n-1;
				else m2 = k2-l2+n;

				if(k1 < l1-n) m1 = k1-l1+3*n+1;
				else if (k1 > l1+n) m1 = k1-l1-n-1;
				else m1 = k1-l1+n;
			 	m	= m2*N + m1;

				Q[k] += B[l*N0+k] * f[l] * f[m];
		}
	}
#endif
}

void UpdateQLongF(int N, const vector<double > & B, const vector<double>& f,
		vector<double>& Q)
{
	TimeWatch<float> tw;
	int N0 = Q.size(); int n = N/2;
	cout << "n = " << n << "\n";
	for(double & x : Q) 	x  = 0.;

	int k, l, l_f, m, m1, m2, m3;
	int kb = 0;
	for(int l3 = 0; l3 < N; ++l3){
	for(int l2 = 0; l2 < N; ++l2){
		for(int l1 = 0; l1 < N; ++l1){
			l = (l3*N+l2)*N+l1;
			l_f = (l3*3*N+l2)*3*N+l1;

			k = 0;
			for(int k3 = 0; k3 < N; ++k3){
				 m3 = k3-l3+n+N;
			for(int k2 = 0; k2 < N; ++k2){
				 m2 = k2-l2+n+N;
				 m = (m3*3*N+m2)*(3*N)-l1+n+N;
				for(int k1 = 0; k1 < N; ++k1){
					//m1 = k1-l1+n+N;
				 	//m = (m3*3*N+m2)*(3*N) + m1;
#if 0
					k = (k3*N+k2)*N + k1;
					Q[k] += B[l*N0+k] * f[l_f] * f[m];
#else
					Q[k] += B[kb] * f[l_f] * f[m];
					kb++, k++,m++;
#endif
				}
			}
		}
	}
	}
	}
}



double l2error(const vector<double>& v1, const vector<double>& v2)
{
	double err = 0;
	for(u_int i = 0; i < v1.size(); ++i)
		err += (v1[i]-v2[i])*(v1[i]-v2[i]);
	return sqrt(err);
}
