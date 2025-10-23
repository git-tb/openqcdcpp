#include "structfunc.h"
#include "pdf.h"
#include "constants.h"
#include "integrate.h"
#include "coefffunc.h"

#include <cmath>		// use math functions from std::

/// @brief electromagnetic structure function F2 from photon interaction
/// @param x hadronic scaling variable x=Q2/(2*P*q) with proton momentum P
/// @param Q2 minus photon momentum squared
double F2(double x, double Q2)	{
	double result_loc(0.);	///< local contributions (coefficient functions \propto \delta(1-z))
	double result_int(0.);	///< non-local contributions (-> integration required)

	/// lo
	result_loc += xfiQi2sum(x, Q2);

	/// @todo running coupling, DONE
	/**
	 * Here's how I understand it so far, based on (Schwartz QFT, eq. 26.102)
	 * Let a  = \alphas_s/4Pi
	 * and aR = \alpha_s^R(\mu)/4Pi	...the renormalized coupling
	 * 
	 * F = F(0) + a * F(1) + a**2 * F(2) + ...
	 * 
	 * For renormalization of amplitudes, replace a -> aR (1 - aR * beta0*ln(Q2/muR2) + aR**2 * (-beta1*ln(Q2/muR2 + beta0**2 * ln**2(Q2/muR2))) + ...)
	 * Organizing F in order of the renormalized coupling we get
	 * 
	 * F = F(0) + aR * F(1) + aR**2 ( F(2) - beta0*ln(Q2/muR2)  F(1) ) + aR**3 ( F(3) - beta0*ln(Q2/muR2) F(2) + (....) F(1))
	 * 
	 * We instead organize F in contributions F(i), each of which comes with a different truncation of the factor from the running coupling.
	 */
	const double muR2 = Q2;
	double a4pi = Pdf::alphas(muR2)/(4.*M_PI);			///< alphas/4Pi at some reference scale muR
	double runcorr_ci_1(1.0);							///< correction to powers of alphasR(mu) associated to N^1LO coefficient function ci_1(...)
	double runcorr_ci_2(1.0);							///< correction to powers of alphasR(mu) associated to N^2LO coefficient function ci_2(...)
	if(QCDORDER::F2ORDER >= 2)	{
		double	logQ2muR2 	= 	std::log(Q2/muR2);
		double 	b0 			= 	QCD::beta0();
		
		runcorr_ci_1 		+= 	- a4pi * b0 * logQ2muR2;

		if(QCDORDER::F2ORDER >= 3)	{
			double b1 		= 	QCD::beta1();

			runcorr_ci_2 	+=	- 2*a4pi * b0 * logQ2muR2;;
			runcorr_ci_1 	+=	a4pi * a4pi * (
								- b1 * logQ2muR2
								+ std::pow(b0 * logQ2muR2,2)
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
			result_int += integrate(
				[x,Q2](double t){return F2integrand_logtrafo2(t,Q2,x);},
				std::log(PRECISION::DELTA), std::log(1-x), PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);
		} else {
			result_int += integrate(
				[x,Q2](double t){return F2integrand_logtrafo1(t,Q2,x);},
				std::log(x), std::log(PRECISION::XTHRESH), PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);
			result_int += integrate(
				[x,Q2](double t){return F2integrand_logtrafo2(t,Q2,x);},
				std::log(PRECISION::DELTA), std::log(1-PRECISION::XTHRESH), PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);
		}
	}

	/// higher orders, local parts
	if(QCDORDER::F2ORDER >= 1)	{
		result_loc += a4pi * runcorr_ci_1 * c2q_ns_1_0_local() * xfiQi2sum(x, Q2);
		result_loc += a4pi * runcorr_ci_1 * c2q_ns_1_0_localplus(x) * xfiQi2sum(x, Q2);
		if(QCDORDER::F2ORDER >= 2)	{
			switch (APPROX::LEVEL)
			{
			case APPROX::APPR1:
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_0_local_approx() * xfiQi2sum(x, Q2);
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_1_local_approx() * xfiQi2sum(x, Q2) *  QCD::NF;

				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_0_localplus_approx(x) * xfiQi2sum(x, Q2);
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_1_localplus_approx(x) * xfiQi2sum(x, Q2) *  QCD::NF;

				result_loc += a4pi * a4pi * runcorr_ci_2 * c2g_2_0_local_approx() * QCD::sumQi2() * Pdf::xf(G, x, Q2);
				break;
			case APPROX::APPR2:
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_0_local_approx2() * xfiQi2sum(x, Q2);
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_1_local_approx2() * xfiQi2sum(x, Q2) *  QCD::NF;

				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_0_localplus_approx2(x) * xfiQi2sum(x, Q2);
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_1_localplus_approx2(x) * xfiQi2sum(x, Q2) *  QCD::NF;

				result_loc += a4pi * a4pi * runcorr_ci_2 * c2g_2_0_local_approx2() * QCD::sumQi2() * Pdf::xf(G, x, Q2);
				break;
			case APPROX::EXACT:
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_01_local_exact() * xfiQi2sum(x, Q2);
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_01_localplus_exact(x) * xfiQi2sum(x, Q2);

				result_loc += a4pi * a4pi * runcorr_ci_2 * c2g_2_0_local_exact() * QCD::sumQi2() * Pdf::xf(G, x, Q2);
				break;
			default:
				std::cout << "ERROR in determining approximation level of coefficient functions" << std::endl;
				abort();
				break;
			}
		}
	}

	return result_loc + result_int;
};


/// @brief Computes all contributions to F2 that need to be numerically integrated over, i.e.
///	all contributions from coefficient functions that are not proportional to delta(1-z).
/// @param z partonic scaling variable
/// @param Q2 minus photon momentum squared
/// @param x hadronic/Bjorken scaling variable
double F2integrand(double z, double Q2, double x)	{
	double result_ns_reg(0.0);
	double result_ns_plus(0.0);
	double result_ps(0.0);
	double result_g(0.0);
	
	///	@todo make the running coupling a parameter of this function such that we do not
	///	recalculate logarithms of the renormalization scale
	/// @todo running coupling, DONE
	const double muR2 = Q2;
	double a4pi = Pdf::alphas(muR2)/(4.*M_PI);	///< alphas/4Pi at some reference scale muR
	double runcorr_ci_1(1.0);							///< correction to powers of alphasR(mu) associated to N^1LO coefficient function ci_1(...)
	double runcorr_ci_2(1.0);							///< correction to powers of alphasR(mu) associated to N^2LO coefficient function ci_2(...)
	if(QCDORDER::F2ORDER >= 2)	{
		double	logQ2muR2 	= 	std::log(Q2/muR2);
		double 	b0 			= 	QCD::beta0();
		
		runcorr_ci_1 		+= 	- a4pi * b0 * logQ2muR2;

		if(QCDORDER::F2ORDER >= 3)	{
			double b1 		= 	QCD::beta1();

			runcorr_ci_2 	+=	- 2*a4pi * b0 * logQ2muR2;;
			runcorr_ci_1 	+=	a4pi * a4pi * (
								- b1 * logQ2muR2
								+ std::pow(b0 * logQ2muR2,2)
			);
		}
	}

	/// reduce number of calls to PDF sampling
	double xfiQi2sumxzQ2		= xfiQi2sum(x/z, Q2);
	double xfiQi2sumxQ2			= xfiQi2sum(x,Q2);
	double xfGxzQ2				= Pdf::xf(G, x/z, Q2);
	double xfiSingletSumxzQ2	= xfiSingletSum(x/z, Q2);

	/// nlo
	if(QCDORDER::F2ORDER >= 1)	{
		result_ns_reg 	+= a4pi * runcorr_ci_1 * c2q_ns_1_0_reg(z) * xfiQi2sumxzQ2;
		result_ns_plus 	+= a4pi * runcorr_ci_1 * c2q_ns_1_0_plus(z) * ( xfiQi2sumxzQ2 - xfiQi2sumxQ2 );
		result_g 		+= a4pi * runcorr_ci_1 * c2g_1_0(z) * QCD::sumQi2() * xfGxzQ2;
		/// nnlo
		if(QCDORDER::F2ORDER >= 2)	{
			switch (APPROX::LEVEL)
			{
			case APPROX::APPR1:
				result_ns_reg 	+= a4pi * a4pi * runcorr_ci_2 * ( c2q_ns_2_0_reg_approx(z) + QCD::NF * c2q_ns_2_1_reg_approx(z) )* xfiQi2sumxzQ2;
				result_ns_plus 	+= a4pi * a4pi * runcorr_ci_2 * ( c2q_ns_2_0_plus_approx(z) + QCD::NF * c2q_ns_2_1_plus_approx(z) ) * ( xfiQi2sumxzQ2 - xfiQi2sumxQ2 ) ;
				result_g 		+= a4pi * a4pi * runcorr_ci_2 * c2g_2_0_reg_approx(z) * QCD::sumQi2() * xfGxzQ2;
				result_ps 		+= a4pi * a4pi * runcorr_ci_2 * c2q_ps_2_0_reg_approx(z) * QCD::sumQi2() * xfiSingletSumxzQ2;
				break;
			case APPROX::APPR2:
				result_ns_reg 	+= a4pi * a4pi * runcorr_ci_2 * ( c2q_ns_2_0_reg_approx2(z) + QCD::NF * c2q_ns_2_1_reg_approx2(z) )* xfiQi2sumxzQ2;
				result_ns_plus 	+= a4pi * a4pi * runcorr_ci_2 * ( c2q_ns_2_0_plus_approx2(z) + QCD::NF * c2q_ns_2_1_plus_approx2(z) ) * ( xfiQi2sumxzQ2 - xfiQi2sumxQ2 ) ;
				result_g 		+= a4pi * a4pi * runcorr_ci_2 * c2g_2_0_reg_approx2(z) * QCD::sumQi2() * xfGxzQ2;
				result_ps 		+= a4pi * a4pi * runcorr_ci_2 * c2q_ps_2_0_reg_approx2(z) * QCD::sumQi2() * xfiSingletSumxzQ2;
				break;
			case APPROX::EXACT:
				result_ns_reg 	+= a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_01_reg_exact(z) * xfiQi2sumxzQ2;
				result_ns_plus 	+= a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_01_plus_exact(z) * ( xfiQi2sumxzQ2 - xfiQi2sumxQ2 ) ;
				result_g 		+= a4pi * a4pi * runcorr_ci_2 * c2g_2_0_reg_exact(z) * QCD::sumQi2() * xfGxzQ2;
				result_ps 		+= a4pi * a4pi * runcorr_ci_2 * c2q_ps_2_0_reg_exact(z) * QCD::sumQi2() * xfiSingletSumxzQ2;
				break;
			default:
				std::cout << "ERROR in determining approximation level of coefficient functions" << std::endl;
				abort();
				break;
			}
		}
	}

	return result_ns_reg + result_ns_plus + result_ps + result_g;
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