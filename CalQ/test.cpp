/**
* Copyright (C) 2021 All rights reserved.
* @file test.cpp
* @brief  
* @author Yang Kuang (matkuan@nus.edu.sg)
* @date 2021-01-26
* @return 
*/

#include <iostream>
#include <vector>
using namespace std;

auto idx2to1(u_int i, u_int  j)
{
	return (i+1)*i/2 + j;
}

int main()
{
	int n = 5;
	vector<double> A (n*(n+1)/2);
	for(u_int i = 0; i < A.size(); ++i)
		A[i] = i+1;

	int k = 0;
	cout << "\t";
	for(int i = 0; i < n; ++i)
		cout << i << "\t";
	cout << "\n";
	for(int i = 0; i < n; ++i){
		cout << i << "\t";
		for(int j = 0; j < i+1; ++j)
			cout << A[k++] << "\t";
		cout << "\n";
	}

	u_int i,j;
	cin >> i >> j;
	auto idx = idx2to1(i,j);
	cout << "A["<<i<<","<<j << "] = " << A[idx] << "\n";
	idx = idx2to1(j,i);
	cout << "A["<<j<<","<<i << "] = " << A[idx] << "\n";
	return 0;
}
