#include <iostream>
#include <functional>	// std::function<()>
#include <cmath>		
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

/// @brief Integration with GSL works through gsl_function objects that can be defined e.g. with lambdas ( i.e. func = [](double, void*){ return something; })
///	but cannot be given access to external variables via the capture bracket [] of a usual lambda. Therefore we must pack all additional variables needed
///	inside the gsl_function into 1 object that is passed as a poitner void* to the function (besides the integral variable).
struct intargs {
    std::function<double(double)>& func;
	double& func1;

	intargs() = delete; /// NO ACCIDENTAL INITIALIZATION
    intargs(std::function<double(double)>& func_, double& func1_) :
                func(func_), func1(func1_) {}
};
/// @brief Integrate the plus distribution of a test function func from 0 to 1.
/// @param func
/// @param del Numerical integration from [0,1-del], approximation on [1-del,1].
/// @return <double> The evaluated integral.
double integratePlus(std::function<double(double)>& func, double del = 1e-10)	{

	/// GSL boilerplate
	int ITERATIONS = 100000;
    int STATUS;
    int KEY(6);	// 61 point Gauss-Kronrod rule https://www.gnu.org/software/gsl/doc/html/integration.html
    gsl_set_error_handler_off();
    gsl_integration_workspace *WORKSPACE = gsl_integration_workspace_alloc(ITERATIONS); 
    gsl_function F;  
    double EPSABS(1e-5), EPSREL(0.0);  
    double RESULT, ERR;

	///	pack test function into GSL object
	double func1 = func(1);
	intargs myintargs(func, func1);
	F.params = &myintargs;
	F.function = [](double x, void* params){
		intargs* args = (struct intargs *)params;
		return (args->func(x)-args->func1)/(1.0-x);
	};

    STATUS = gsl_integration_qag(&F, 0, 1-del, EPSABS, EPSREL, ITERATIONS, KEY, WORKSPACE, &RESULT, &ERR);
    if (STATUS)	std::cerr << gsl_strerror(STATUS) << std::endl;

	return RESULT + 2.0*(func(1-del/2.0)-func1);
}

int main(int argc, char** argv) {
	for(uint i = 0; i < (uint)argc; i++)	{
		std::cout << "arg" << i << ": " << argv[i] << std::endl;
	}
	if(argc != 3)	{
		std::cout << "[ERROR] usage: ./plusdist <double>A <double>del" << std::endl;
		return 1;
	}
	double A, del;
	try	{
		A = std::stod(argv[1]);
		del = std::stod(argv[2]);
	}	catch(const std::exception& e)	{
		std::cout << "[ERROR] usage: ./plusdist <double>A <double>del" << std::endl;
		return 1;
	}

	std::function<double(double)> testfunc = [A](double z)	{
		return 1/((z-(1+A))*(z-(1+A)));
	};

	std::cout << integratePlus(testfunc,del) << std::endl;

	return 0;
}