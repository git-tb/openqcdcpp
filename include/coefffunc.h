#ifndef COEFFFUNC_H
#define COEFFFUNC_H

#include <cmath>		// use math functions from std::

#include "constants.h"

/**
 * @brief The naming of the coefficient functions follows the scheme
 * 		cnf_ps/ns_j_k_desc
 * with	c - "coefficient function"
 * 		n - for structure function fn
 * 		f - quark (q) or gluon (g)
 * 		ps/ns - pure singlet (ps) or non-singlet (ns) in case of f=q
 * 		j - associated QCD order (alpha_s/(4*Pi))^j
 * 		k - associated power of Nf
 * 		desc - further description (i.e. regular, plus distribution, ...)
 */



///
///
/// O(alphas)

/// @brief from [Zijlstra, van Neerven; Nucl. Phys. B. 383 (3), 525, 1992] eq. B.2
double c2q_ns_1_0_reg(double z)	{
	double result(0.);
	result += 	QCD::CF * (								
					- 2. * (1.+z) * std::log(1.-z)			
					- 2. * (1.+z*z)/(1.-z) * std::log(z)	
					+ 6. + 4.*z ); 
	return result;
}

/// @brief from [Zijlstra, van Neerven; Nucl. Phys. B. 383 (3), 525, 1992] eq. B.2
double c2q_ns_1_0_plus(double z)	{
	double result(0.);
	result += QCD::CF * ( 4. * std::log(1.-z)/(1.-z) - 3./(1.-z) );
	return result;
}

/// @brief from [Zijlstra, van Neerven; Nucl. Phys. B. 383 (3), 525, 1992] eq. B.2
double c2q_ns_1_0_local()	{
	double result(0.);
	result += - QCD::CF * ( 2.*M_PI*M_PI/3. + 9.);
	return result;
}

/// @brief from [Zijlstra, van Neerven; Nucl. Phys. B. 383 (3), 525, 1992] eq. B.2
///	contains the local contribution from integrating the plus distributions over (x,1)
///	instead of (0,1)
double c2q_ns_1_0_localplus(double x)	{
	double result(0.);
	double lx = std::log(1-x);
	result += QCD::CF * ( 2. * lx*lx - 3. * lx );
	return result;
}

/// @brief from [Zijlstra, van Neerven; Nucl. Phys. B. 383 (3), 525, 1992] eq. B.6
double c2g_1_0(double z)	{
	double result(0.);
	result += QCD::TR * 4 * (
					+ (1. - 2.*z + 2.*z*z) * std::log((1.-z)/z)
					- 1 + 8.*z*(1.-z)
				);
	return result;
}










///
///
/// O(alphas^2)

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_0_local_approx()	{	
	return -338.046;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_1_local_approx()	{	
	return 46.8405;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_0_localplus_approx(double z)	{
	double result(0.0);
	double L1 = std::log(1-z);
	result +=	14.222 * std::pow(L1,4) / 4.0
				- 61.3333 * std::pow(L1,3) / 3.0
				- 31.105 * std::pow(L1,2) / 2.0
				+ 188.64 * L1;
	return result;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_1_localplus_approx(double z)	{
	double result(0.0);
	double L1 = std::log(1-z);
	result +=	1.77778 * std::pow(L1, 3) / 3.0
				- 8.5926 * std::pow(L1,2) / 2.0
				+ 6.3489 * L1;
	return result;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_0_plus_approx(double z)	{
	double result(0.0);
	double L1 = std::log(1-z);
	result +=	(
				14.2222 * std::pow(L1,3)
				- 61.3333 * std::pow(L1,2)
				- 31.105 * L1
				+ 188.64 ) / (1.0 - z);
	return result;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_1_plus_approx(double z)	{
	double result(0.0);
	double L1 = std::log(1-z);
	result +=	(
				1.77778 * std::pow(L1, 2)
				- 8.5926 * L1
				+ 6.3489 ) / (1.0 - z);
	return result;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_0_reg_approx(double z)	{
	double result(0.0);
	double L0 = std::log(z);
	double L1 = std::log(1-z);
	result +=	- 17.19 * std::pow(L1,3)
				+ 71.08 * std::pow(L1,2)
				- 660.7 * L1
				+ L1 * L0 * ( -174.8 * L1 + 95.09 * L0 )
				- 2.835 * std::pow(L0, 3)
				- 17.08 * std::pow(L0, 2)
				+ 5.986 * L0
				- 1008 * z
				- 69.59;
	return result;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_1_reg_approx(double z)	{
	double result(0.0);
	double L0 = std::log(z);
	double L1 = std::log(1-z);
	result +=	- 1.707 * std::pow(L1, 2)
				+ 22.95 * L1
				+ L1 * L0 * ( 3.036 * L0 + 17.97 )
				+ 2.244 * std::pow(L0,2)
				+ 5.770 * L0
				- 37.91 * z
				- 5.691;
	return result;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 588 (1), 345, 2000] eq. 3.6
double c2g_2_0_reg_approx(double z)	{
	double result(0.0);
	double L0 = std::log(z);
	double L1 = std::log(1-z);
	result +=	( 6.445 + 209.4 * (1-z)) * std::pow(L1,3)
				- 24.00 * std::pow(L1,2)
				+ (1494.0/z - 1483) * L1
				+ L1 * L0 * (-871.8 * L1 - 724.1 * L0)
				+ 5.319 * std::pow(L0,3)
				- 59.48 * std::pow(L0,2)
				- 284.8 * L0
				+ 11.90 / z 
				+ 392.4;
	return result;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 588 (1), 345, 2000] eq. 3.6
double c2g_2_0_local_approx()	{
	return -0.28;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 588 (1), 345, 2000] eq. 3.6
double c2q_ps_2_0_reg_approx(double z)	{
	double result(0.0);
	double L0 = std::log(z);
	double L1 = std::log(1-z);
	result +=	- 0.101 * (1 - z) * std::pow(L1,3) 
				- (24.75 - 13.80 * z) * L1 * std::pow(L0,2)
				+ 30.23 * L1 * L0
				+ 4.310 * std::pow(L0,3)
				- 2.086 * std::pow(L0,2)
				+ 39.78 * L0
				+ 5.290 * (1.0/z - 1);
	return result;
}

#endif