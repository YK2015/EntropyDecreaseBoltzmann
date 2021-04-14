/**
* Copyright (C) 2021 All rights reserved.
* @file fft3d.h
* @brief  input real data, output frequency data, and vice visa
* @author Yang Kuang (matkuan@nus.edu.sg)
* @date 2021-04-13
* @return 
*/

#ifndef __FFT3D_H__
#define __FFT3D_H__

#include <cmath>
#include <complex>
#include <fftw3.h>

class fft3d
{
	public:
		fft3d(){};
		~fft3d(){};

	public:
		/// fftw_complex = std::complex<double> 
		void fft(int Nx, int Ny, int Nz, fftw_complex* re , fftw_complex* fre){
			fftw_plan p;
			p = fftw_plan_dft_3d(Nx,Ny,Nz, re, fre, FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(p);
			fftw_destroy_plan(p);
		}

		void ifft(int Nx, int Ny, int Nz, fftw_complex* re , fftw_complex* fre){
			fftw_plan p;
			p = fftw_plan_dft_3d(Nx,Ny,Nz, fre, re, FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(p);
			fftw_destroy_plan(p);

			/// normalization
		  int N=Nx*Ny*Nz;
			int k = 0;
			for(int i= 0; i < Nx; ++ i)
				for(int j = 0; j < Ny; ++j)
					for(int l = 0; l < Ny; ++l){
						re[k] *= 1./N;++k;
					}
		}
};





#endif ///FFT3D_H__

