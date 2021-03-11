/**
* Copyright (C) 2021 All rights reserved.
* @file TimeWatch.h
* @brief  
* @author Yang Kuang (matkuan@nus.edu.sg)
* @date 2021-01-06
* @return 
*/
#include <sys/time.h>
#include <chrono>
#include <iostream>
#include <iomanip>

#ifndef __TIMEWATCH_H__
#define __TIMEWATCH_H__

template <typename timeType>
class TimeWatch
{
	std::chrono::steady_clock::time_point stime;	/// start time
	timeType ctime;	

	public:
	TimeWatch(){stime=std::chrono::steady_clock::now();ctime=clock();};
	~TimeWatch(){
		auto etime = std::chrono::steady_clock::now();
		auto d = std::chrono::duration_cast<std::chrono::microseconds>(etime - stime);
		ctime=(timeType)(clock()-ctime)/CLOCKS_PER_SEC;
		timeType wtime = timeType(d.count()/1000000.);
		std::cout <<std::fixed << std::setprecision(2) 
							<< "#\tWTime: " << wtime << "s\t"
							<< "CTime: " << ctime << "s\n"
							<< "#\tSpeedup: " << ctime/wtime << "\n";
	};
	
};



#endif ///TIMEWATCH_H__

