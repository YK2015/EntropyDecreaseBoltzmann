#include "QCalculator.h"
#include <stdlib.h>
#include <random>
#include <complex>



int main(int argc, char** argv)
{
	typedef double  r_type;
	typedef std::complex<double>  c_type;
	
	int N0 = std::atoi(argv[1]);

	int N = N0*N0*N0;
	std::vector<c_type> fin(N), Qout(N), Qout2(N);
	for (auto & x : fin) {
		x = c_type{static_cast<double>(rand() % 10), static_cast<double>(rand()%10)};
	}

	
	std::cout << "# Fast evaluation...\n";
	QCalculator<c_type,r_type> qc(N0,0);
	{
	  std::cout << "# FE::genB\n";
	//  TimeWatch<r_type> tw;
	  qc.genB();
	}
	{
	  std::cout << "# FE::calQ\n";
	 // TimeWatch<r_type> tw;
	  qc.calQ(Qout2,fin);
	}
	std::cout << "# End fast evaluation.\n\n";




	return 0;
}

