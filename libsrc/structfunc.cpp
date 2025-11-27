#include "structfunc.h"
#include "pdf.h"
#include "constants.h"
#include "integrate.h"
#include "coefffunc.h"
#include "hqcoeff.h"

#include <cmath>		// use math functions from std::

/// @brief electromagnetic structure function F2 from photon interaction
/// @param x hadronic scaling variable x=Q2/(2*P*q) with proton momentum P
/// @param Q2 minus photon momentum squared
///	@param orderAlps highest included power of alpha_s
double F2(double x, double Q2, int orderAlps, int nf)	{
	double result_loc(0.);	///< local contributions (coefficient functions \propto \delta(1-z))
	double result_int(0.);	///< non-local contributions (-> integration required)

	/// lo
	result_loc += Pdf::xfiQi2sum(x, Q2);

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
	if(orderAlps >= 2)	{
		double	logQ2muR2 	= 	std::log(Q2/muR2);
		double 	b0 			= 	QCD::beta0(nf);
		
		runcorr_ci_1 		+= 	- a4pi * b0 * logQ2muR2;

		if(orderAlps >= 3)	{
			double b1 		= 	QCD::beta1(nf);

			runcorr_ci_2 	+=	- 2*a4pi * b0 * logQ2muR2;;
			runcorr_ci_1 	+=	a4pi * a4pi * (
								- b1 * logQ2muR2
								+ std::pow(b0 * logQ2muR2,2)
			);
		}
	}

	/// higher orders, non-local parts
	if(orderAlps >= 1)	{
		//// naive way, samples uniformly across integration domain
		// result_int += integrate(
		// 	[x,Q2](double z){return F2integrand(z,Q2,x);},
		// 	x, 1, PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);

		if(x >= PRECISION::XTHRESH)	{
			result_int += integrate(
				[x,Q2,orderAlps,nf](double t){return F2integrand_logtrafo2(t,Q2,x,orderAlps, nf);},
				std::log(PRECISION::DELTA), std::log(1-x), PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);
		} else {
			result_int += integrate(
				[x,Q2,orderAlps,nf](double t){return F2integrand_logtrafo1(t,Q2,x,orderAlps, nf);},
				std::log(x), std::log(PRECISION::XTHRESH), PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);
			result_int += integrate(
				[x,Q2,orderAlps,nf](double t){return F2integrand_logtrafo2(t,Q2,x,orderAlps, nf);},
				std::log(PRECISION::DELTA), std::log(1-PRECISION::XTHRESH), PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL);
		}
	}

	/// higher orders, local parts
	if(orderAlps >= 1)	{
		result_loc += a4pi * runcorr_ci_1 * c2q_ns_1_0_local() * Pdf::xfiQi2sum(x, Q2);
		result_loc += a4pi * runcorr_ci_1 * c2q_ns_1_0_localplus(x) * Pdf::xfiQi2sum(x, Q2);
		if(orderAlps >= 2)	{
			switch (APPROX::LEVEL)
			{
			case APPROX::APPR1:
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_0_local_approx() * Pdf::xfiQi2sum(x, Q2);
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_1_local_approx() * Pdf::xfiQi2sum(x, Q2) *  nf;

				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_0_localplus_approx(x) * Pdf::xfiQi2sum(x, Q2);
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_1_localplus_approx(x) * Pdf::xfiQi2sum(x, Q2) *  nf;

				result_loc += a4pi * a4pi * runcorr_ci_2 * c2g_2_0_local_approx() * QCD::sumQi2(nf) * Pdf::xf(G, x, Q2);
				break;
			case APPROX::APPR2:
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_0_local_approx2() * Pdf::xfiQi2sum(x, Q2);
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_1_local_approx2() * Pdf::xfiQi2sum(x, Q2) *  nf;

				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_0_localplus_approx2(x) * Pdf::xfiQi2sum(x, Q2);
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_1_localplus_approx2(x) * Pdf::xfiQi2sum(x, Q2) *  nf;

				result_loc += a4pi * a4pi * runcorr_ci_2 * c2g_2_0_local_approx2() * QCD::sumQi2(nf) * Pdf::xf(G, x, Q2);
				break;
			case APPROX::EXACT:
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_01_local_exact(nf) * Pdf::xfiQi2sum(x, Q2);
				result_loc += a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_01_localplus_exact(x,nf) * Pdf::xfiQi2sum(x, Q2);

				result_loc += a4pi * a4pi * runcorr_ci_2 * c2g_2_0_local_exact() * QCD::sumQi2(nf) * Pdf::xf(G, x, Q2);
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



double F2integrand(double z, double Q2, double x, int orderAlps, int nf)	{
	double result_ns_reg(0.0);
	double result_ns_plus(0.0);
	double result_ps(0.0);
	double result_g(0.0);
	
	///	@todo make the running coupling a parameter of this function such that we do not
	///	recalculate logarithms of the renormalization scale
	/// @todo DONE running coupling
	const double muR2 = Q2;
	double a4pi = Pdf::alphas(muR2)/(4.*M_PI);	///< alphas/4Pi at some reference scale muR
	double runcorr_ci_1(1.0);					///< correction to powers of alphasR(mu) associated to N^1LO coefficient function ci_1(...)
	double runcorr_ci_2(1.0);					///< correction to powers of alphasR(mu) associated to N^2LO coefficient function ci_2(...)
	if(orderAlps >= 2)	{
		double	logQ2muR2 	= 	std::log(Q2/muR2);
		double 	b0 			= 	QCD::beta0(nf);
		
		runcorr_ci_1 		+= 	- a4pi * b0 * logQ2muR2;

		if(orderAlps >= 3)	{
			double b1 		= 	QCD::beta1(nf);

			runcorr_ci_2 	+=	- 2*a4pi * b0 * logQ2muR2;;
			runcorr_ci_1 	+=	a4pi * a4pi * (
								- b1 * logQ2muR2
								+ std::pow(b0 * logQ2muR2,2)
			);
		}
	}

	/// reduce number of calls to PDF sampling
	double xfiQi2sumxzQ2		= Pdf::xfiQi2sum(x/z, Q2);
	double xfiQi2sumxQ2			= Pdf::xfiQi2sum(x,Q2);
	double xfGxzQ2				= Pdf::xf(G, x/z, Q2);
	double xfiSingletSumxzQ2	= Pdf::xfiSingletSum(x/z, Q2);

	/// nlo
	if(orderAlps >= 1)	{
		result_ns_reg 	+= a4pi * runcorr_ci_1 * c2q_ns_1_0_reg(z) * xfiQi2sumxzQ2;
		result_ns_plus 	+= a4pi * runcorr_ci_1 * c2q_ns_1_0_plus(z) * ( xfiQi2sumxzQ2 - xfiQi2sumxQ2 );
		result_g 		+= a4pi * runcorr_ci_1 * c2g_1_0(z) * QCD::sumQi2(nf) * xfGxzQ2;
		/// nnlo
		if(orderAlps >= 2)	{
			switch (APPROX::LEVEL)
			{
			case APPROX::APPR1:
				result_ns_reg 	+= a4pi * a4pi * runcorr_ci_2 * ( c2q_ns_2_0_reg_approx(z) + nf * c2q_ns_2_1_reg_approx(z) )* xfiQi2sumxzQ2;
				result_ns_plus 	+= a4pi * a4pi * runcorr_ci_2 * ( c2q_ns_2_0_plus_approx(z) + nf * c2q_ns_2_1_plus_approx(z) ) * ( xfiQi2sumxzQ2 - xfiQi2sumxQ2 ) ;
				result_g 		+= a4pi * a4pi * runcorr_ci_2 * c2g_2_0_reg_approx(z) * QCD::sumQi2(nf) * xfGxzQ2;
				result_ps 		+= a4pi * a4pi * runcorr_ci_2 * c2q_ps_2_0_reg_approx(z) * QCD::sumQi2(nf) * xfiSingletSumxzQ2;
				break;
			case APPROX::APPR2:
				result_ns_reg 	+= a4pi * a4pi * runcorr_ci_2 * ( c2q_ns_2_0_reg_approx2(z) + nf * c2q_ns_2_1_reg_approx2(z) )* xfiQi2sumxzQ2;
				result_ns_plus 	+= a4pi * a4pi * runcorr_ci_2 * ( c2q_ns_2_0_plus_approx2(z) + nf * c2q_ns_2_1_plus_approx2(z) ) * ( xfiQi2sumxzQ2 - xfiQi2sumxQ2 ) ;
				result_g 		+= a4pi * a4pi * runcorr_ci_2 * c2g_2_0_reg_approx2(z) * QCD::sumQi2(nf) * xfGxzQ2;
				result_ps 		+= a4pi * a4pi * runcorr_ci_2 * c2q_ps_2_0_reg_approx2(z) * QCD::sumQi2(nf) * xfiSingletSumxzQ2;
				break;
			case APPROX::EXACT:
				result_ns_reg 	+= a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_01_reg_exact(z, nf) * xfiQi2sumxzQ2;
				result_ns_plus 	+= a4pi * a4pi * runcorr_ci_2 * c2q_ns_2_01_plus_exact(z, nf) * ( xfiQi2sumxzQ2 - xfiQi2sumxQ2 ) ;
				result_g 		+= a4pi * a4pi * runcorr_ci_2 * c2g_2_0_reg_exact(z) * QCD::sumQi2(nf) * xfGxzQ2;
				result_ps 		+= a4pi * a4pi * runcorr_ci_2 * c2q_ps_2_0_reg_exact(z) * QCD::sumQi2(nf) * xfiSingletSumxzQ2;
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
double F2integrand_logtrafo1(double t, double Q2, double x, int orderAlps, int nf)	{
	double z = std::exp(t);
	return F2integrand(z,Q2,x,orderAlps,nf)*z;
}

/// @brief transformation of F2integrand that samples closer to large z
double F2integrand_logtrafo2(double t, double Q2, double x, int orderAlps, int nf)	{
	double z = 1.-std::exp(t);
	return F2integrand(z,Q2,x,orderAlps,nf)*(1-z);
}



/// @todo fix all issues with factorization scale, basically in all functions
double F2heavy(double x, double Q2, int orderAlps, int nlight)	{
	double result(0.0);
	double a = 1./(1.+4.*std::pow(QCD::QMASSES[nlight],2)/Q2);
	double muR2 = Q2 + 4*std::pow(QCD::QMASSES[nlight],2);

	if(x >= a) return 0.;

	result += integrate([x,Q2,nlight,orderAlps, muR2](double logz){ return F2heavyintegrand(std::exp(logz), Q2, x, orderAlps, nlight, muR2);},
		std::log(x), std::log(a), PRECISION::ITER, PRECISION::EPSABS, PRECISION::EPSREL
	);

	return result;
}

double F2heavyintegrand(double z, double Q2, double x, int orderAlps, int nlight)	{
	return F2heavyintegrand(z,Q2,x,orderAlps, nlight,Q2);
}

/// @brief implementation of [Laenen, Riemersma, Smith, v. Neerven; Nucl. Phys. B, 392 (1), 162, 1993]
/// equations (6.2) @ LO, (6.4) & (6.5) @ NLO
double F2heavyintegrand(double z, double Q2, double x, int orderAlps, int nlight, double muR2)	{
	double result(0.0);

	/// @todo move multiplicative prefactors out of this function
	/// @todo call PDFs only once
	///	@todo make the running coupling a parameter of this function such that we do not
	///	recalculate logarithms of the renormalization scale
	double alps = Pdf::alphas(muR2);
	double runcorr_ci_1(1.0);							///< correction to powers of alphasR(mu) associated to O(alphas) coefficient function ci_1(...)
	double runcorr_ci_2(1.0);							///< correction to powers of alphasR(mu) associated to O(alphas^2) coefficient function ci_2(...)
	// if(QCDORDER::F2ORDER >= 2)	{
	// 	double	logQ2muR2 	= 	std::log(Q2/muR2);
	// 	double 	b0 			= 	QCD::beta0();
		
	// 	runcorr_ci_1 		+= 	- alps / (4.0*M_PI) * b0 * logQ2muR2;

	// 	if(QCDORDER::F2ORDER >= 3)	{
	// 		double b1 		= 	QCD::beta1();

	// 		runcorr_ci_2 	+=	- 2.0*alps/(4.0*M_PI) * b0 * logQ2muR2;;
	// 		runcorr_ci_1 	+=	alps/(4.0*M_PI) * alps/(4.0*M_PI) * (
	// 							- b1 * logQ2muR2
	// 							+ std::pow(b0 * logQ2muR2,2)
	// 		);
	// 	}
	// }

	double s	= Q2*(1.0/z-1.0);
	double eta	= s/(4.0*std::pow(QCD::QMASSES[nlight],2))-1.0;
	double chi	= Q2/std::pow(QCD::QMASSES[nlight],2);

	// LO is O(alphas)	
	if(orderAlps >= 1)	{
		result += Q2 * alps * runcorr_ci_1 / (4 * M_PI*M_PI) * std::pow(QCD::QCHARGES[nlight]/QCD::QMASSES[nlight],2) * (
			Pdf::xf(G,x/z,muR2) * ch2g_1_0(eta,chi)
		);

		if(orderAlps >= 2)	{
			double L	= std::log(muR2/std::pow(QCD::QMASSES[nlight],2));

			/// gluon
			result	+=  Q2 * std::pow(alps,2) / (M_PI * std::pow(QCD::QMASSES[nlight], 2)) * (
				std::pow(QCD::QCHARGES[nlight],2) * Pdf::xf(G,x/z,muR2) * (
					ch2_g_2_0(eta,chi) + L * ch2_g_2_1(eta,chi)
				)
			);

			/// quark
			result	+=  Q2 * std::pow(alps,2) / (M_PI * std::pow(QCD::QMASSES[nlight], 2)) * (
				std::pow(QCD::QCHARGES[nlight],2) * Pdf::xfiSingletSum(x/z,muR2) * (
					ch2_q_2_0_Hcoupl(eta,chi) + L * ch2_q_2_1_Hcoupl(eta,chi)
				)
			);

			result	+=  Q2 * std::pow(alps,2) / (M_PI * std::pow(QCD::QMASSES[nlight], 2)) * (
				Pdf::xfiQi2sum(x/z,muR2) * (
					ch2_q_2_0_Lcoupl(eta,chi) + L * ch2_q_2_1_Lcoupl(eta,chi)
				)
			);

			if(orderAlps >= 3)	{
				result += Q2 * std::pow(alps,3) * 4 / std::pow(QCD::QMASSES[nlight], 2) * (
					std::pow(QCD::QCHARGES[nlight],2) * Pdf::xf(G,x/z,muR2) * (
						ch2_g_3_0(eta,chi) + L * ch2_g_3_1(eta,chi) + L * L * ch2_g_3_2(eta,chi)
					)
				);
			}
		}
	}

	return result;
}