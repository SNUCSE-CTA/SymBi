/**
 * @date 2020.03.25
 * @author Seunghwan Min
 */

#include <ctime>
#include <iostream>
#include <chrono>

using namespace std; 

class Timer{
private:
	chrono::high_resolution_clock::time_point s, e;
	double time;
public:
	Timer() : time(0.0) {}
	~Timer() {}

	void start()
	{
		s= chrono::high_resolution_clock::now();
	}

	void end()
	{
		e = chrono::high_resolution_clock::now();
		time += chrono::duration<double, milli> (e - s).count();
	}

	friend ostream& operator << (ostream& out, const Timer& t)
	{
		out << t.time;
		return out;
	}
};

