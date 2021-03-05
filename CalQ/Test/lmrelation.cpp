#include <iostream>
#include <vector>
#include <stdlib.h>
using namespace std;

int main(int argc, char ** argv)
{
	int n = atoi(argv[1]);
	int N = 2*n + 1;
	vector<double> Q(N);
	for(int i = 0; i < N; ++i) Q[i] = i-n;
	vector<double> f(2*N);
	for(int i = 0; i < n+1; ++i) f[i] = i;
	for(int i = n+1; i < N+n+1; ++i) f[i] = i-N;
	for(int i = N+n+1; i < 2*N; ++i) f[i] = i-4*n-2;


	cout << "Q:\n";
	for(int i = 0; i < Q.size(); ++i)
		cout << Q[i] << "\t";

	cout << "\nf:\n";
	for(int i = 0; i < f.size(); ++i)
		cout << f[i] << "\t";
	cout << "\n\n";


	cout << "k = l + m \n";
	for(int l = n+1; l < 3*n + 2; ++l){
		for(int k = n+1; k < 3*n + 2; ++k){
			int m = k - l + N; 	
			cout << Q[k-n-1] << " += " << f[l] << " + " << f[m] << "\t";
		}
		cout << "\n";
	}

	cout << "\nStorage of B \n";
	for(int l = n+1; l < 3*n + 2; ++l){
		for(int k = n+1; k < 3*n + 2; ++k){
			int m = k - l + N; 	
			cout << "B[" << f[l]+n << ", " << f[m]+n << "]\t";
		}
		cout << "\n";
	}

	cout << "\nStorage of B \n";
	for(int l = 0; l < N; ++l){
		for(int k = n-l; k < n-l+N; ++k){
			int m; 	
			if(k  < 0) m = k+N;
			else m = k%N;	
			cout << "B[" << l << ", " << m << "]\t";
		}
		cout << "\n";
	}




	
	return 0;
}
