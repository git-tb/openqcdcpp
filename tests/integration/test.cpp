#include <iostream>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <vector>

#include "integrate.h"

extern "C" {
    void gauss1_(double (*F)(double*), double* A, double* B,
                 int* LAMBDA, double* RESULT, double* EPS);
}


double testfunc(double *x)	{
	return 1./std::sqrt(*x);
}
const double trueresult = 2;

/**
 * @brief This program tests different integration routines: 
 * 
 * gauss1(), a standard (12-point?) gaussian quadrature rule implemented in Fortran in openQCDrad
 * 	and
 * gls_integrate(), integration in the GSL library
 * 
 * The goal is to increase number of subdivisions step by step and monitor
 * the accuracy and computational cost.
 */
int main(int argc, char** argv)	{

	const int NRUN = 1e4;	///< number of integration runs over which the runtime is averaged
	
	/// timing variables
	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	/// integration 
	double result(0.0);
	double err(0.0);
	double x0(0.0), x1(1.0);
	int NDIV(1);				///< number of subintervals, modifed in the loop
	double EPSABS(1e-30);
	double EPSREL(1e-30);

	/// output formatting
	const int PREC = 20;

	/// output file
	std::ofstream fileout("testdata.dat");

	fileout << std::setprecision(PREC);

	fileout 	<< "runtime;"
				<< "err;"
				<< "algorithm"
				<< std::endl;
	
	gsl_set_error_handler_off();
	std::vector<int> NDIV_list = {
		1,2,3,4,5,6,7,8,9,10,
		20,40,60,80,100,
		200,400,600,800,1000,
		2000,4000
	};
	for(int NDIV: NDIV_list)	{
		std::cout << "running... NDIV=" << NDIV << std::endl;

		/// gsl adaptive
		start = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < NRUN; i++)	{
			int KEY 	= 1;
			int ITER	= NDIV;

			gsl_integration_workspace *WORKSPACE = gsl_integration_workspace_alloc(ITER); 
    		gsl_function F;
			F.function	= [](double x, void* params) { return testfunc(&x);};

			int STATUS;
			STATUS = gsl_integration_qag(&F, x0, x1, EPSABS, EPSREL, ITER, KEY, WORKSPACE, &result, &err);
			// if (STATUS)	std::cerr << gsl_strerror(STATUS) << std::endl;

			gsl_integration_workspace_free(WORKSPACE);
		}
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		fileout	<< (double)duration.count()/(double)NRUN << ";"
				<< trueresult - result << ";"
				<< "gsl_qag"
				<< std::endl;

		/// gsl adaptive, reuse workspace
		start = std::chrono::high_resolution_clock::now();
		int KEY 	= 1;
		int ITER	= NDIV;		
		gsl_integration_workspace *WORKSPACE = gsl_integration_workspace_alloc(ITER); 
		gsl_function F;
		F.function	= [](double x, void* params) { return testfunc(&x);};	
		int STATUS;
		for(int i = 0; i < NRUN; i++)	{
			STATUS = gsl_integration_qag(&F, x0, x1, EPSABS, EPSREL, ITER, KEY, WORKSPACE, &result, &err);
		}
		gsl_integration_workspace_free(WORKSPACE);
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		fileout	<< (double)duration.count()/(double)NRUN << ";"
				<< trueresult - result << ";"
				<< "gsl_qag_ReuseWorkspace"
				<< std::endl;

		/// gsl non-adaptive
		start = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < NRUN; i++)	{
			ulong ITER	= NDIV;

			gsl_integration_workspace *WORKSPACE = gsl_integration_workspace_alloc(ITER); 
    		gsl_function F;
			F.function	= [](double x, void* params) { return testfunc(&x);};

			int STATUS;
			STATUS = gsl_integration_qng(&F, x0, x1, EPSABS, EPSREL, &result, &err, &ITER);

			gsl_integration_workspace_free(WORKSPACE);
		}
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		fileout	<< (double)duration.count()/(double)NRUN << ";"
				<< trueresult - result << ";"
				<< "gsl_qng"
				<< std::endl;

		/// openQCDrad
		start = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < NRUN; i++)	{
			gauss1_(testfunc, &x0, &x1, &NDIV, &result, &err);
		}
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		fileout	<< (double)duration.count()/(double)NRUN << ";"
				<< trueresult - result << ";"
				<< "gauss1_(openQCDrad)"
				<< std::endl;
				
		/// integrate
		start = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < NRUN; i++)	{
			result = integrate([](double x){return testfunc(&x);}, x0, x1, NDIV, EPSABS, EPSREL);
		}
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		fileout	<< (double)duration.count()/(double)NRUN << ";"
				<< trueresult - result << ";"
				<< "integrate(openQCD++)"
				<< std::endl;
	}
	fileout.close();

	return 0;
}