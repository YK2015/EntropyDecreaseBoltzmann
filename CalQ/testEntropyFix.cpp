#include "EntropyFix.h"
#include "TimeWatch.h"
#include <stdlib.h>
#include <cmath>
//#include <random>
//#include <complex>



int main(int argc, char** argv)
{
	typedef double  r_type;
	typedef int	i_type;	
	i_type N = std::atoi(argv[1]);
	//i_type N = 6;
	/*
	for (auto & x : fin) {
		std::cin >> x;
	}
	*/

	std::vector<r_type> f1(N), f2(N);
	//std::cout << "f^{n+1} = ";
	for (i_type i = 0; i < N; i++) {
		f1[i] = (r_type)i;
		//std::cout << f1[i] << " ";
	}
	//std::cout << '\n';

	//std::cout << "f(t_{n+1}) = ";
	for (i_type i = 0; i < N/2; i++) {
		f2[i] = floor( ((r_type)N - 1.0) / 2.0);
		//std::cout << f2[i] << " ";
	}
	for (i_type i = N/2; i < N; i++) {
		f2[i] = ceil( ((r_type)N - 1.0) / 2.0);
		//std::cout << f2[i] << " ";
	}
	//std::cout << '\n';	
	
	EntropyFix<i_type,r_type> ef(N);
	r_type entf1;
	{
		TimeWatch<r_type> tw;
		entf1 = ef.EntropyCal(f1);
	}
	std::cout << "H(f^{n+1}) = " << entf1 << std::endl;
	
	r_type entf2;
	entf2 = ef.EntropyCal(f2);
	std::cout << "H(f(t_{n+1})) = " << entf2 << std::endl;

	ef.FixSol(f1, entf2);
	/*
	std::cout << "f^{n+1} is revised to be: ";
	for (auto & x : f1) {
		std::cout << x << ' ';
	}
	std::cout << '\n';
	*/

	std::cout << "H(hatf^{n+1}) = " << ef.EntropyCal(f1) << std::endl;
	/*
	{
	  std::cout << "# FE::calQ\n";
	  TimeWatch<r_type> tw;
	  qc.calQ(Qout2,fin);
	}
	*/
	return 0;
}

