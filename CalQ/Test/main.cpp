#include "calBBox.h"
#include <stdlib.h>
#include <random>
#include <complex>

#include "TimeWatch.h"


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

	
	std::cout << "# Direct evaluation...\n";
	directEvaluation<c_type,r_type> de(N0,0);
	{
	  std::cout << "# DE::genB\t";
	  TimeWatch<r_type> tw;
	  de.genB();
	}
	{ 
	  TimeWatch<r_type> twatch;
	  de.calQ(Qout,fin);
	}

	
	std::cout <<"# End direct evaluation.\n\n";

	std::cout << "# Fast evaluation...\n";
	QCalculator<c_type,r_type> qc(N0,0);
	{
	  std::cout << "# FE::genB\t";
	  TimeWatch<r_type> tw;
	  qc.genB();
	}
	{
	  TimeWatch<r_type> tw;
	  qc.calQ(Qout2,fin);
	}
	std::cout << "# End fast evaluation.\n\n";


//	std::cout << "i\texact\tqc\t|qc-exact|\n";
	r_type l2err = 0., tmperr;
	for(int i = 0; i < N; ++i){
		tmperr = std::abs(Qout[i]-Qout2[i]);
//		std::cout << i << "\t" << Qout[i] 
//			<< "\t" << Qout2[i] << "\t" << tmperr 
//			<< "\n";
		l2err += tmperr*tmperr;
	}
	std::cout << "L2err = " << l2err << "\n";



	return 0;
}

