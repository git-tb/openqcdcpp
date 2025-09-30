#ifndef STRUCTFUNC_H
#define STRUCTFUNC_H

#include "pdf.h"
#include "constants.h"
#include "integrate.h"
#include "coefffunc.h"

#include <cmath>		// use math functions from std::

double F2(double x, double Q2);
double F2integrand(double z, double Q2, double x);
double F2integrand_logtrafo1(double t, double Q2, double x);
double F2integrand_logtrafo2(double t, double Q2, double x);

/// @brief electromagnetic structure function F2 from photon interaction
/// @param x hadronic scaling variable x=Q2/(2*P*q) with proton momentum P
/// @param Q2 minus photon momentum squared
double F2(double x, double Q2)	{
	double result(0.);

	/// lo
	result += xfiQi2sum(x, Q2);

	/// @todo running coupling, DONE
	const double muR2 = Q2;
	double a4pi = Pdf::get()->alphasQ2(muR2)/(4.*M_PI);	///< alphas/4Pi at some reference scale muR
	double runcorr(1.0);
	if(QCDORDER::F2ORDER >= 1)	{
		double logQ2muR2 = std::log(Q2/muR2);
		double b0 = QCD::beta0(), b1 = QCD::beta1();
		runcorr += - a4pi * b0 * logQ2muR2;

		if(QCDORDER::F2ORDER >= 2)	{
			runcorr += a4pi * a4pi * (
				+	b1 * logQ2muR2
				-	std::pow(b0 * logQ2muR2,2)
			);
		}
	}

	/// higher orders, non-local parts
	if(QCDORDER::F2ORDER >= 1)	{
		//// naive way, samples uniformly across integration domain
		// result += integrate(
		// 	[x,Q2](double z){return F2integrand(z,Q2,x);},
		// x, 1, PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);

		if(x >= PRECISION::XTHRESH)	{
			result += integrate(
				[x,Q2](double t){return F2integrand_logtrafo2(t,Q2,x);},
				std::log(PRECISION::DELTA), std::log(1-x), PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);
		} else {
			result += integrate(
				[x,Q2](double t){return F2integrand_logtrafo1(t,Q2,x);},
				std::log(x), std::log(PRECISION::XTHRESH), PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);
			result += integrate(
				[x,Q2](double t){return F2integrand_logtrafo2(t,Q2,x);},
				std::log(PRECISION::DELTA), std::log(1-PRECISION::XTHRESH), PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);
		}
	}

	/// higher orders, local parts
	if(QCDORDER::F2ORDER >= 1)	{
		result += a4pi * runcorr * c2q_ns_1_0_local() * xfiQi2sum(x, Q2);
		result += a4pi * runcorr * c2q_ns_1_0_localplus(x) * xfiQi2sum(x, Q2);
	}

	return result;
};


/// @brief Computes all contributions to F2 that need to be numerically integrated over, i.e.
///	all contributions from coefficient functions that are not proportional to delta(1-z).
/// @param z partonic scaling variable
/// @param Q2 minus photon momentum squared
/// @param x hadronic/Bjorken scaling variable
double F2integrand(double z, double Q2, double x)	{
	double result(0.);
	
	/// @todo running coupling, DONE
	///	@todo make the running coupling a parameter of this function such that we do not
	///	recalculate logarithms of the renormalization scale
	const double muR2 = Q2;
	double a4pi = Pdf::get()->alphasQ2(muR2)/(4.*M_PI);	///< alphas/4Pi at some reference scale muR
	double runcorr(1.0);
	if(QCDORDER::F2ORDER >= 1)	{
		double logQ2muR2 = std::log(Q2/muR2);
		double b0 = QCD::beta0(), b1 = QCD::beta1();
		runcorr += - a4pi * b0 * logQ2muR2;

		if(QCDORDER::F2ORDER >= 2)	{
			runcorr += a4pi * a4pi * (
				+	b1 * logQ2muR2
				-	std::pow(b0 * logQ2muR2,2)
			);
		}
	}

	
	/// nlo
	if(QCDORDER::F2ORDER >= 1)	{
		result += a4pi * runcorr * c2q_ns_1_0_reg(z) * xfiQi2sum(x/z, Q2);
		result += a4pi * runcorr * c2q_ns_1_0_plus(z) * ( xfiQi2sum(x/z, Q2) - xfiQi2sum(x, Q2) );
		result += a4pi * runcorr * c2g_1_0(z) * QCD::sumQi2() * Pdf::get()->xfxQ2(G, x/z, Q2);
	}

	return result;
}

/// @brief transformation of F2integrand that samples closer to small z
double F2integrand_logtrafo1(double t, double Q2, double x)	{
	double z = std::exp(t);
	return F2integrand(z,Q2,x)*z;
}

/// @brief transformation of F2integrand that samples closer to large z
double F2integrand_logtrafo2(double t, double Q2, double x)	{
	double z = 1.-std::exp(t);
	return F2integrand(z,Q2,x)*(1-z);
}

#endif