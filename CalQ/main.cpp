#include "calBBox.h"
#include <stdlib.h>
#include <random>

int main(int argc, char** argv)
{
	
	int N = std::atoi(argv[1]);
	QCalculator<double,3> qcalculator(N);
	qcalculator.genB(0);

	std::vector<double> fin(N*N*N), Qout(N*N*N);
	for (auto & x : fin) x = rand() % 10;

	qcalculator.calQ(Qout,fin);
	qcalculator.displayQ();

	return 0;
}

