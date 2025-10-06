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
 * is passed as a poitner void* to the function (besides the integration variable). 
 */
struct intargs {
    const std::function<double(double)>& func;

	intargs() = delete; /// no accidental initialization
    intargs(const std::function<double(double)>& func_) :
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
 * @todo reduce overhead from copying the function object -> switch to function pointer?
 * @todo reduce overhead from recreating the gsl_integration_workspace every time integrate() is called
 */
double integrate(const std::function<double(double)>& func, double a, double b, int ITER, double EPSABS, double EPSREL);


inline gsl_integration_workspace *STATICWORKSPACE = NULL;
/**
 * @brief Same integration routine, but this time with a workspace that is only recreated
 * if its too small
 */
double integrate_STATICWORKSPACE(const std::function<double(double)>& func, double a, double b, int ITER, double EPSABS, double EPSREL);

#endif