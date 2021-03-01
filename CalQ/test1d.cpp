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

void refUpdateQ(const vector<double>  &, const vector<double>&,
		vector<double>&);

void UpdateQ(const vector<double > &, const vector<double>&,
		vector<double>&);

void UpdateQLongF(const vector<double > &, const vector<double>&,
		vector<double>&);

double l2error(const vector<double>&, const vector<double>& );

int main(int argc, char * argv[])
{
	u_int N = atoi(argv[1]);
	vector<double> f(N);
	vector<double> lf(3*N);
	vector<double >  B(N*N);
	for(u_int i = 0;i < N; ++i){
		f[i] =  i* 0.1; 
		lf[i] = f[i]; 
		lf[i+N] = f[i]; 
		lf[i+2*N] = f[i]; 
		for(u_int j = 0; j < N; ++j)
			B[i+j*N] = (i+j)* 0.1;
	}

	vector<double> Q(N);
	cout << "# Classical method \n";
	refUpdateQ(B,f,Q);

	vector<double> Q2(N);
	cout << "# Update in a new rule \n";
	UpdateQ(B,f,Q2);

	vector<double> Q3(N);
	cout << "# Update with Trippled f \n";
	UpdateQLongF(B,lf,Q3);


	auto err = l2error(Q2,Q3);
	cout << "l2 err = " << err << "\n";

	return 0;
}

void refUpdateQ(const vector<double>  & B, const vector<double>& f,
		vector<double>& Q)
{
	TimeWatch<float> tw;
	int N = Q.size();int n = N/2;
	cout << "n = " << n << "\n";
	for(double & x : Q) 	x  = 0.;

#if 1
	int m;
	for(int k = 0; k < N; ++k){
		for(int l = 0; l < N; ++l){
			if(l > k+n) m = k-l+3*n+1;
			else if (l < k-n) m = k-l-n-1;
			else m = k-l+n;
			Q[k] += B[k*N+l] * f[l] * f[m];
			/*
			if(l > k+n)
				Q[k] += B[k*N+l] * f[l] * f[k-l+3*n+1];
			else if(l < k-n)
				Q[k] += B[k*N+l] * f[l] * f[k-l-n-1];
			else
				Q[k] += B[k*N+l] * f[l] * f[k-l+n];
				*/
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

void UpdateQ(const vector<double > & B, const vector<double>& f,
		vector<double>& Q)
{
	TimeWatch<float> tw;
	int N = Q.size(); int n = N/2;
	for(double & x : Q) 	x  = 0.;
#if 0
	for(int l = 0; l < N; ++l){
		for(int k = 0; k < N; ++k){
			if(k < l-n)
				Q[k] += B[l*N+k] * f[l] * f[k-l+3*n+1];
			else if(k > l+n)
				Q[k] += B[l*N+k] * f[l] * f[k-l-n-1];
			else
				Q[k] += B[l*N+k] * f[l] * f[k-l+n];
		}
	}
#else
	for(int l = 0; l < n; ++l){
		for(int k = 0; k < l+n+1; ++k)
			Q[k] += B[l*N+k] * f[l] * f[k-l+n];
		for(int k = l+n+1; k < N; ++k)
			Q[k] += B[l*N+k] * f[l] * f[k-l-n-1];
	} 
	for(int l = n; l < N; ++l){
		for(int k = 0; k < l-n; ++k)
			Q[k] += B[l*N+k] * f[l] * f[k-l+3*n+1];
		for(int k = l-n; k < N; ++k)
			Q[k] += B[l*N+k] * f[l] * f[k-l+n];
	}
#endif
}

void UpdateQLongF(const vector<double > & B, const vector<double>& f,
		vector<double>& Q)
{
	TimeWatch<float> tw;
	int N = Q.size(); int n = N/2;
	for(double & x : Q) 	x  = 0.;
#if 1
	for(int l = 0; l < N; ++l){
		for(int k = 0; k < N; ++k){
			Q[k] += B[l*N+k] * f[l] * f[k - l +n+N];
		}
	}
#else
	for(u_int l = 0; l < n; ++l){
		for(u_int k = 0; k < l+n+1; ++k)
			Q[k] += B[l*N+k] * f[l] * f[k-l+n];
		for(u_int k = l+n+1; k < N; ++k)
			Q[k] += B[l*N+k] * f[l] * f[k-l-n-1];
	} 
	for(u_int l = n; l < N; ++l){
		for(u_int k = 0; k < l-n; ++k)
			Q[k] += B[l*N+k] * f[l] * f[k-l+3*n+1];
		for(u_int k = l-n; k < N; ++k)
			Q[k] += B[l*N+k] * f[l] * f[k-l+n];
	}
#endif

}



double l2error(const vector<double>& v1, const vector<double>& v2)
{
	double err = 0;
	for(u_int i = 0; i < v1.size(); ++i)
		err += (v1[i]-v2[i])*(v1[i]-v2[i]);
	return sqrt(err);
}
