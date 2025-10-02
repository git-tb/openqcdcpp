#include <iostream>
#include <chrono>
#include <cmath>
#include <iomanip>

#include "integrate.h"

extern "C" {
    void gauss1_(double (*F)(double*), double* A, double* B,
                 int* LAMBDA, double* RESULT, double* EPS);
}

double testfunc(double *x)	{
	return std::sin(*x)*std::sin(*x);
}

int main(int argc, char** argv)	{

	const int NRUN = 1e4;
	double total(0.0);

	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);


	double result(0.0);
	double eps(0.0);
	double a(0.0), b(M_PI);
	int NDIV(1);
	double EPSABS(1e-10);
	double EPSREL(1e-10);


	const int WIDTH = 20;
	const int PREC = 15;
	std::cout << std::setprecision(PREC) << std::endl;

	start = std::chrono::high_resolution_clock::now();
	for(int i = 0; i < NRUN; i++)	{
		integrate([](double x){return testfunc(&x);}, a, b, NDIV, EPSABS, EPSREL);
	}
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	std::cout 	<< std::setw(15) << (double)duration.count()/(double)NRUN
				<< std::setw(3) << "ns"
				<< std::setw(WIDTH)<< integrate([](double x){return testfunc(&x);}, a, b, NDIV, EPSABS, EPSREL) 
				<< std::endl;


	start = std::chrono::high_resolution_clock::now();
	for(int i = 0; i < NRUN; i++)	{
		gauss1_(testfunc, &a, &b, &NDIV, &result, &eps);
	}
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	std::cout 	<< std::setw(15) << (double)duration.count()/(double)NRUN 
				<< std::setw(3) << "ns"
				<< std::setw(WIDTH) << result
				<< std::endl;

	std::cout 	<< std::setw(15+3+WIDTH) << M_PI/2.0 << std::endl;
	return 0;
}