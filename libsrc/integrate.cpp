#include "integrate.h"

double integrate(const std::function<double(double)>& func, double a, double b, int ITER, double EPSABS, double EPSREL)	{
	/// GSL boilerplate
    int STATUS;
    int KEY(1);	///< 61 point Gauss-Kronrod rule https://www.gnu.org/software/gsl/doc/html/integration.html

    gsl_set_error_handler_off();
    gsl_integration_workspace *WORKSPACE = gsl_integration_workspace_alloc(ITER); 
    gsl_function F;
    double RESULT, ERR;

	///	pack test function into GSL object
	intargs myintargs(func);
	F.params = &myintargs;
	F.function = [](double x, void* params){
		intargs* args = (struct intargs *)params;
		return args->func(x);
	};

    STATUS = gsl_integration_qag(&F, a, b, EPSABS, EPSREL, ITER, KEY, WORKSPACE, &RESULT, &ERR);
    if (STATUS)	std::cerr << gsl_strerror(STATUS) << std::endl;
	
	gsl_integration_workspace_free(WORKSPACE);

	return RESULT;
}

double integrate_STATICWORKSPACE(const std::function<double(double)>& func, double a, double b, int ITER, double EPSABS, double EPSREL)	{
	/// GSL boilerplate
    int STATUS;
    int KEY(1);	///< 61 point Gauss-Kronrod rule https://www.gnu.org/software/gsl/doc/html/integration.html


	if((not STATICWORKSPACE) || (STATICWORKSPACE->limit != ITER))	{
		gsl_set_error_handler_off();
		STATICWORKSPACE = gsl_integration_workspace_alloc(ITER); 
		std::cout << "workspace created!" << std::endl;
	}
    double RESULT, ERR;
	
	///	pack test function into GSL object
	gsl_function F;
	intargs myintargs(func);
	F.params = &myintargs;
	F.function = [](double x, void* params){
		intargs* args = (struct intargs *)params;
		return args->func(x);
	};

    STATUS = gsl_integration_qag(&F, a, b, EPSABS, EPSREL, ITER, KEY, STATICWORKSPACE, &RESULT, &ERR);
    if (STATUS)	std::cerr << gsl_strerror(STATUS) << std::endl;
	
	return RESULT;
}