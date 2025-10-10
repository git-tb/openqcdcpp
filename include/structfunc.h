#ifndef STRUCTFUNC_H
#define STRUCTFUNC_H

#include "pdf.h"
#include "constants.h"
#include "integrate.h"
#include "coefffunc.h"

#include <cmath>		// use math functions from std::

/// @brief electromagnetic structure function F2 from photon interaction
/// @param x hadronic scaling variable x=Q2/(2*P*q) with proton momentum P
/// @param Q2 minus photon momentum squared
double F2(double x, double Q2);

/// @brief Computes all contributions to F2 that need to be numerically integrated over, i.e.
///	all contributions from coefficient functions that are not proportional to delta(1-z).
/// @param z partonic scaling variable
/// @param Q2 minus photon momentum squared
/// @param x hadronic/Bjorken scaling variable
double F2integrand(double z, double Q2, double x);

/// @brief transformation of F2integrand that samples closer to small z
double F2integrand_logtrafo1(double t, double Q2, double x);

/// @brief transformation of F2integrand that samples closer to large z
double F2integrand_logtrafo2(double t, double Q2, double x);

#endif