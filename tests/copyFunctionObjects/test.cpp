#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <functional>
#include <chrono>
#include <cmath>
#include <iomanip>

struct intargs {
    std::function<double(double)>* func;

	intargs() = delete; /// no accidental initialization
    intargs(std::function<double(double)>* func_) :
                func(func_) {}
};

double integrate(std::function<double(double)>& func, double a, double b, int ITER, double EPSABS, double EPSREL)	{
	/// GSL boilerplate
    int STATUS;
    int KEY(1);	///< 61 point Gauss-Kronrod rule https://www.gnu.org/software/gsl/doc/html/integration.html

    gsl_set_error_handler_off();
    gsl_integration_workspace *WORKSPACE = gsl_integration_workspace_alloc(ITER); 
    gsl_function F;
    double RESULT, ERR;

	///	pack test function into GSL object
	intargs* myintargs = new intargs(&func);
	F.params = myintargs;
	F.function = [](double x, void* params){
		intargs* args = (struct intargs *)params;
		return (*args->func)(x);
	};
	
    STATUS = gsl_integration_qag(&F, a, b, EPSABS, EPSREL, ITER, KEY, WORKSPACE, &RESULT, &ERR);
    // if (STATUS)	std::cerr << gsl_strerror(STATUS) << std::endl;
	
	delete myintargs;
	gsl_integration_workspace_free(WORKSPACE);

	return RESULT;
}

int main()	{
	int ITER		= 1000;
	double EPSABS	= 1e-3;
	double EPSREL	= 1e-3;
	double x0		= 0;
	double x1		= 5;
	uint NRUN		= (uint)1e6;
	auto start 		= std::chrono::high_resolution_clock::now();
	auto stop 		= std::chrono::high_resolution_clock::now();
	auto duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	double time1 	= (double)duration.count()/(double)NRUN;
	double time2 	= (double)duration.count()/(double)NRUN;

	double a		= 1.;
	std::function<double(double)> testfunc = [a](double x){ return std::sin(a*x);};
	
	start 		= std::chrono::high_resolution_clock::now();
	for(uint i = 0; i < NRUN; i++)	{
		std::function<double(double)> loctestfunc = [a](double x){ return std::sin(a*x);};
		integrate(loctestfunc,x0,x1,ITER,EPSABS,EPSREL);
	}
	stop 		= std::chrono::high_resolution_clock::now();
	duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	time1 		= (double)duration.count()/(double)NRUN;
	
	start 		= std::chrono::high_resolution_clock::now();
	for(uint i = 0; i < NRUN; i++)	{
		integrate(testfunc,x0,x1,ITER,EPSABS,EPSREL);
	}
	stop 		= std::chrono::high_resolution_clock::now();
	duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	time2 		= (double)duration.count()/(double)NRUN;

	uint WIDTH	= 20;
	uint PREC	= 10;
	std::cout	<< std::setprecision(10);
	std::cout	<< std::setw(WIDTH)	<< time1 << " ns"	<< std::endl
				<< std::setw(WIDTH)	<< time2 << " ns"	<< std::endl;

	return 0;
}