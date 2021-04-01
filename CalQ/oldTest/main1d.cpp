#include "calBBox1D.h"
#include <stdlib.h>
#include <random>
#include <iomanip>

#include "TimeWatch.h"


int main(int argc, char** argv)
{
	typedef double  value_type;
	
	int N0 = std::atoi(argv[1]);

	int N = N0;
	std::vector<value_type> fin(N), Qout(N), Qout2(N);
	for (auto & x : fin) x = rand() % 10;

	
	directEvaluation<value_type> de(N0,0);
	de.genB();
	std::cout << "Entry of B\n\t\t";
	for(int i = 0; i < N; ++i)
		std::cout << i - N/2 << "\t\t";
	std::cout << "\n";
	for(int i = 0; i < N; ++i){
	  std::cout << i - N/2 << "\t";
	  for(int j = 0 ; j < N; ++j){
		std::cout << std::fixed<< std::setprecision(8)
		  << de.B[i][j] << "\t";
	  }
	  std::cout << "\n";
	}
/*
	{ 
	  TimeWatch<double> twatch;
	  de.calQ(Qout,fin);
	}
*/
	
/*
	std::cout << "Fast evaluation...\n";
	QCalculator<value_type> qc(N0,0);
	std::cout << "FE::genB\n";
	qc.genB();
	std::cout << "FE::calQ\n";
	{
	  TimeWatch<double> tw;
	  qc.calQ(Qout2,fin);
	}
	std::cout << "End fast evaluation.\n";


	//std::cout << "i\texact\tqc\t|qc-exact|\n";
	value_type l2err = 0., tmperr;
	for(int i = 0; i < N; ++i){
		tmperr = std::abs(Qout[i]-Qout2[i]);
	//	std::cout << i << "\t" << Qout[i] 
	//		<< "\t" << Qout2[i] << "\t" << tmperr 
	//		<< "\n";
		l2err += tmperr*tmperr;
	}
	std::cout << "L2err = " << l2err << "\n";
*/


	return 0;
}

