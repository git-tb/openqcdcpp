#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <iostream>
#include <functional>	/// std::function
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

/**
 * @brief Integration with GSL works through gsl_function objects that can be defined e.g.
 * with lambdas ( i.e. func = [](double, void*){ return something; }) but cannot be given
 * access to external variables via the capture bracket [] of a usual lambda. Therefore we
 * must pack all additional variables needed inside the gsl_function into 1 arg object that
 * is passed as a pointer void* to the function (besides the integration variable). 
 */

struct intargs {
    std::function<double(double)>* func;

	intargs() = delete; /// no accidental initialization
    intargs(std::function<double(double)>* func_) :
                func(func_) {}
};

/**
 * @brief light weight interface to gsl_integration that integrates a function func from a to b.
 * 
 * @param func monovariate function
 * @param a lower integration limit
 * @param b upper integration limit
 * @param ITER number of subdivisions
 * @param EPSABS absolute precision goal
 * @param EPSREL relative precision goal
 * 
 * @todo reduce overhead from copying the function object
 *			-> switch to function pointer?
 *			-> at the moment this is difficult because we have to pass
 *			parameters [x,Q2] to our functions which require a complete
 *			std::function<...> object. We could switch to a call signature
 *			given by 
 *				>> integrate(double(*func)(double,double,double), double par1, double par2,...)
 *			making this wrapper less useful.
 * @todo (YIELDS NO SPEEDUP) reduce overhead from recreating the gsl_integration_workspace every time integrate() is called
 */
double integrate(std::function<double(double)> func, double a, double b, int ITER, double EPSABS, double EPSREL);

#endif