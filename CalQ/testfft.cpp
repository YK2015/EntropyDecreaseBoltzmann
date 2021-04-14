#include "QCalculator.h"
#include "fft3d.h"
#include "TimeWatch.h"
#include <stdlib.h>
#include <random>

typedef double  r_type;
typedef std::complex<double>  c_type;
	
/// domain [-6,6]^N
r_type fun(r_type t, r_type x, r_type y, r_type z)
{
	r_type v2 = (x*x+y*y+z*z);
	r_type s = 1.-2.*exp(-t/6.)/5.;
	return 1./std::pow(2.*M_PI*s,3./2.)*exp(-v2/(2.*s)) * ((5.*s-3.)/(2.*s)+ (1.-s)/(2.*s*s)*v2);
}


int main(int argc, char** argv)
{

	int N0 = std::atoi(argv[1]);

	int N = N0*N0*N0;

	fftw_complex F[N], Fhat[N], Ftest[N];
	
	double h = 12./double(N0-1);
	int l = 0;
	double x,y,z;
	for(int i = 0; i < N0; ++i){
		x = -6. + i*h;
		for(int j = 0; j < N0; ++j){
			y = -6. + j*h;
			for(int k = 0; k < N0; ++k){
				z = -6. + k*h;
				F[l] = fun(0,x,y,z);
				l++;
			}
		}
	}

	/// fft
	fft3d f3d;
	f3d.fft(N0,N0,N0,F,Fhat);
	f3d.ifft(N0,N0,N0,Ftest,Fhat);

  l = 0;
	for(int i = 0; i < N0; ++i)
		for(int j = 0; j < N0; ++j)
			for(int k = 0; k < N0; ++k){
				std::cout << "("<<i<<","<<j<<","<<k<<")\t";
				std::cout << Fhat[l] << "\t";
				std::cout << F[l] << "\t" << Ftest[l] << "\n";
				l ++;
			}




	return 0;
}

