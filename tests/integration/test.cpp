#include <iostream>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "integrate.h"

extern "C" {
    void gauss1_(double (*F)(double*), double* A, double* B,
                 int* LAMBDA, double* RESULT, double* EPS);
}

// const double A = 5;
// const double B = 10;
// double testfunc(double *x)	{
// 	return std::exp((*x)*A)*std::cos((*x)*B);
// }
// const double trueresult = (std::exp(A)*B*std::sin(B) + std::exp(A)*A*std::cos(B) - A)/(A*A+B*B);


double testfunc(double *x)	{
	return 1./std::sqrt(*x);
}
const double trueresult = 2;

/**
 * @brief This program tests 2 integration routines: 
 * 
 * gauss1(), a standard gaussian quadrature rule implemented in Fortran in openQCDrad
 * 	and
 * gls_integrate(), integration in the GSL library
 * 
 * The goal is to increase number of subdivisions step by step and monitor
 * the accuracy and computational cost.
 */
int main(int argc, char** argv)	{

	const int NRUN = 1e4;	///< number of integration runs over which the runtime is averaged
	
	/// timing variables
	double total(0.0);
	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	/// integration 
	double result(0.0);
	double eps(0.0);
	double x0(0.0), x1(1);
	int NDIV(1);				///< number of subintervals, modifed in the loop
	double EPSABS(1e-25);
	double EPSREL(1e-30);

	/// output formatting
	const int WIDTH = 23;
	const int PREC = 20;
	// std::cout << std::setprecision(PREC) << std::endl;

	/// output files
	std::ofstream fileout_cpp("integrate_gsl.dat");		///< C++ integration routine
	std::ofstream fileout_for("integrate_gauss1.dat");	///< openQCDrad integration routine

	fileout_cpp << std::setprecision(PREC);
	fileout_for << std::setprecision(PREC);

	fileout_cpp 	<< std::setw(WIDTH) << "runtime(avg)"
					<< std::setw(WIDTH) << ";"
					<< std::setw(WIDTH) << "integrationresult"
					<< std::endl;
	fileout_for 	<< std::setw(WIDTH) << "runtime(avg)"
					<< std::setw(WIDTH) << ";"
					<< std::setw(WIDTH) << "integrationresult"
					<< std::endl;
	
	for(int ndivpow = 0; ndivpow <= 12; ndivpow ++)	{
		NDIV = (int)std::round(std::pow(2,ndivpow));
		std::cout << "running... NDIV=" << NDIV << std::endl;

		/// cpp routine
		start = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < NRUN; i++)	{
			// integrate([](double x){return testfunc(&x);}, x0, x1, NDIV, EPSABS, EPSREL);
			integrate_STATICWORKSPACE([](double x){return testfunc(&x);}, x0, x1, NDIV, EPSABS, EPSREL);
		}
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		// std::cout 	<< std::setw(15) << (double)duration.count()/(double)NRUN
		// 			<< std::setw(3) << "ns"
		// 			<< std::setw(WIDTH)<< integrate([](double x){return testfunc(&x);}, a, b, NDIV, EPSABS, EPSREL) 
		// 			<< std::endl;

		fileout_cpp << std::setw(WIDTH) << (double)duration.count()/(double)NRUN
					<< std::setw(WIDTH) << ";"
					<< std::setw(WIDTH) << integrate([](double x){return testfunc(&x);}, x0, x1, NDIV, EPSABS, EPSREL)
					<< std::endl;

		/// fortran routine
		start = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < NRUN; i++)	{
			gauss1_(testfunc, &x0, &x1, &NDIV, &result, &eps);
		}
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		// std::cout 	<< std::setw(15) << (double)duration.count()/(double)NRUN 
		// 			<< std::setw(3) << "ns"
		// 			<< std::setw(WIDTH) << result
		// 			<< std::endl;

		fileout_for << std::setw(WIDTH) << (double)duration.count()/(double)NRUN
					<< std::setw(WIDTH) << ";"
					<< std::setw(WIDTH) << result
					<< std::endl;
	}
	// std::cout 	<< std::setw(15+3+WIDTH) << M_PI/2.0 << std::endl;

	fileout_cpp.close();
	fileout_for.close();

	std::cout << std::setprecision(PREC);
	std::cout << std::setw(20) << "True result:" << std::setw(WIDTH) << trueresult << std::endl;

	return 0;
}