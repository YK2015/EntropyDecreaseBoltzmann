#include "calBBox.h"
#include <stdlib.h>
#include <random>

int main(int argc, char** argv)
{
	
	int N0 = std::atoi(argv[1]);

//	int N = N0;
	int N = N0*N0*N0;
	std::vector<double> fin(N), Qout(N), Qout2(N);
	for (auto & x : fin) x = rand() % 10;


	directEvaluation<double> de(N0,0);
	de.genB();
	de.calQ(Qout,fin);

	QCalculator<double> qc(N,0);
	qc.genB();
	qc.calQ(Qout2,fin);


	std::cout << "i\texact\tqc\t|qc-exact|\n";
	double l2err = 0., tmperr;
	for(int i = 0; i < N; ++i){
		tmperr = fabs(Qout[i]-Qout2[i]);
		std::cout << i << "\t" << Qout[i] 
			<< "\t" << Qout2[i] << "\t" << tmperr 
			<< "\n";
		l2err += tmperr*tmperr;
	}
	std::cout << "L2err = " << l2err << "\n";



	return 0;
}

