#include "coefffunc.h"
#include <cmath>		// use math functions from std::

#include "constants.h"
#include "chaplin.h"

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
	result += - QCD::CF * ( 4.*MATH::ZETA2 + 9.);
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

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.8
double c2q_ns_2_0_local_approx2()	{	
	return -338.513;
}

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.8
double c2q_ns_2_1_local_approx2()	{	
	return 46.8531;
}

/// @brief eq. (B.5) in [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005]
/// This expression was obtained after some rewriting and reordering with Mathematica.
double c2q_ns_2_01_local_exact()	{
	double result(0.0);
	result += 	std::pow(
QCD::CF,2)*(331./8. + 69*MATH::ZETA2 + 6*std::pow(MATH::ZETA2,2) - \
78*MATH::ZETA3) + 
QCD::CA*
QCD::CF*(-5465./72. + -251./3.*MATH::ZETA2 + \
71./5.*std::pow(MATH::ZETA2,2) + 140./3.*MATH::ZETA3) + 
QCD::CF*
QCD::NF*(457./36. + 38./3.*MATH::ZETA2 + 4./3.*MATH::ZETA3);
	return result;
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

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.8
double c2q_ns_2_0_localplus_approx2(double z)	{
	double result(0.0);
	double L1 = std::log(1-z);
	result +=	128./9. * std::pow(L1,4) / 4.
				- 184./3. * std::pow(L1,3) / 3.
				- 31.1052 * std::pow(L1,2) / 2.
				+ 188.64 * L1;
	return result;
}

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.8
double c2q_ns_2_1_localplus_approx2(double z)	{
	double result(0.0);
	double L1 = std::log(1-z);
	result +=	16./9. * std::pow(L1, 3) / 3.
				- 232./27. * std::pow(L1,2) / 2.
				+ 6.34888 * L1;
	return result;
}

/// @brief eq. (B.5) in [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005]
/// This expression was obtained after some rewriting and reordering with Mathematica.
double c2q_ns_2_01_localplus_exact(double x)	{
	double result(0.0);
	double L1x	=	std::log(1-x);
	result 		+= 	QCD::CA*
QCD::CF*((-3155./54. + 44./3.*MATH::ZETA2 + 40*MATH::ZETA3)*L1x + \
(367./18. - 4*MATH::ZETA2)*std::pow(L1x,2) + -22./9.*std::pow(L1x,3)) \
+ 
QCD::CF*
QCD::NF*((247./27. + -8./3.*MATH::ZETA2)*L1x + -29./9.*std::pow(L1x,2) \
+ 4./9.*std::pow(L1x,3)) + std::pow(
QCD::CF,2)*((51./2. + 36*MATH::ZETA2 - 8*MATH::ZETA3)*L1x + (-27./2. - \
16*MATH::ZETA2)*std::pow(L1x,2) - 6*std::pow(L1x,3) + \
2*std::pow(L1x,4));
	return result;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_0_plus_approx(double z)	{
	double result(0.0);
	double L1 = std::log(1-z);
	result +=	( 14.2222 * std::pow(L1,3)
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

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.8
double c2q_ns_2_0_plus_approx2(double z)	{
	double result(0.0);
	double L1 = std::log(1-z);
	result +=	( 128./9. * std::pow(L1,3)
				- 184./3. * std::pow(L1,2)
				- 31.1052 * L1
				+ 188.641 ) / (1.0 - z);
	return result;
}

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.8
double c2q_ns_2_1_plus_approx2(double z)	{
	double result(0.0);
	double L1 = std::log(1-z);
	result +=	( 16./9. * std::pow(L1, 2)
				- 232./27. * L1
				+ 6.34888 ) / (1.0 - z);
	return result;
}

/// @brief eq. (B.5) in [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005]
/// This expression was obtained after some rewriting and reordering with Mathematica.
double c2q_ns_2_01_plus_exact(double z)	{
	double result(0.0);
	double L1z	= std::log(1-z);
	result += (
QCD::CA*
QCD::CF*(3155./54. + -44./3.*MATH::ZETA2 - 40*MATH::ZETA3 + (-367./9. \
+ 8*MATH::ZETA2)*L1z + 22./3.*std::pow(L1z,2)) + 
QCD::CF*
QCD::NF*(-247./27. + 8./3.*MATH::ZETA2 + 58./9.*L1z + \
-4./3.*std::pow(L1z,2)) + std::pow(
QCD::CF,2)*(-51./2. - 36*MATH::ZETA2 + 8*MATH::ZETA3 + (27 + \
32*MATH::ZETA2)*L1z + 18*std::pow(L1z,2) - 8*std::pow(L1z,3)))/(-1 + \
z);
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

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.8
double c2q_ns_2_0_reg_approx2(double z)	{
	double result(0.0);
	double L0 = std::log(z);
	double L1 = std::log(1-z);
	result +=	- 17.74 * std::pow(L1,3)
				+ 72.24 * std::pow(L1,2)
				- 628.8 * L1
				- 181.0
				- 806.7 * z
				+ 0.719 * z * std::pow(L0,4)
				+ L1 * L0 * ( 37.75 * L0 - 147.1 * L1 )
				- 28.384 * L0
				- 20.70 * std::pow(L0, 2)
				- 80./27. * std::pow(L0, 3);
	return result;
}

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.8
double c2q_ns_2_1_reg_approx2(double z)	{
	double result(0.0);
	double L0 = std::log(z);
	double L1 = std::log(1-z);
	result +=	- 1.500 * std::pow(L1, 2)
				+ 24.87 * L1
				- 7.8109
				- 17.82 * z
				- 12.97 * std::pow(z,2)
				- 0.185 * z * std::pow(L0,3)
				+ 8.113 * L0 * L1
				+ 16./3. * L0
				+ 20./9. * std::pow(L0,2);
	return result;
}

/// @brief eq. (B.5) in [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005]
/// This expression was obtained after some rewriting and reordering with Mathematica.
double c2q_ns_2_01_reg_exact(double z)	{
	double result(0.0);
	double L1z		= std::log(1-z);
	double HPL1_m1z	= HPL1(-1,z);
	double HPL1_0z	= HPL1(0,z);
	double HPL1_1z	= HPL1(1,z);
	double HPL2_m10z	= HPL2(-1,0,z);
	double HPL2_00z	= HPL2(0,0,z);
	double HPL2_01z	= HPL2(0,1,z);
	double HPL2_10z	= HPL2(1,0,z);
	double HPL2_11z	= HPL2(1,1,z);
	double HPL3_m1m10z	= HPL3(-1,-1,0,z);
	double HPL3_m100z	= HPL3(-1,0,0,z);
	double HPL3_m101z	= HPL3(-1,0,1,z);
	double HPL3_0m10z	= HPL3(0,-1,0,z);
	double HPL3_000z	= HPL3(0,0,0,z);
	double HPL3_001z	= HPL3(0,0,1,z);
	double HPL3_010z	= HPL3(0,1,0,z);
	double HPL3_011z	= HPL3(0,1,1,z);
	double HPL3_100z	= HPL3(1,0,0,z);
	double HPL3_101z	= HPL3(1,0,1,z);
	double HPL3_110z	= HPL3(1,1,0,z);
	double HPL3_111z	= HPL3(1,1,1,z);
	result += QCD::CA*
QCD::CF*(3709./135. + -8./5./z + 17626./135.*z + -72./5.*std::pow(z,2) \
+ (12 - 56*z)*MATH::ZETA3 + (164./5.*HPL1_0z)/(-1 + std::pow(z,2)) + \
(-8./5.*HPL1_0z)/(z*(-1 + std::pow(z,2))) + (-118./5.*z*HPL1_0z)/(-1 + \
std::pow(z,2)) + (799./15.*std::pow(z,2)*HPL1_0z)/(-1 + std::pow(z,2)) \
+ (1693./15.*std::pow(z,3)*HPL1_0z)/(-1 + std::pow(z,2)) + \
(-72./5.*std::pow(z,4)*HPL1_0z)/(-1 + std::pow(z,2)) + 56./9.*HPL1_1z \
+ 668./9.*z*HPL1_1z + MATH::ZETA2*(-44./3. + -104./3.*z + \
72./5.*std::pow(z,3) - (8*z*HPL1_0z)/(-1 + std::pow(z,2)) - \
(8*std::pow(z,3)*HPL1_0z)/(-1 + std::pow(z,2)) - 8*HPL1_1z - \
32*z*HPL1_1z) - 36*HPL2_m10z + (-8./5.*HPL2_m10z)/std::pow(z,2) - \
20*z*HPL2_m10z + 72./5.*std::pow(z,3)*HPL2_m10z + 55./3.*HPL2_00z + \
115./3.*z*HPL2_00z + -72./5.*std::pow(z,3)*HPL2_00z + 44./3.*HPL2_01z \
+ 44./3.*z*HPL2_01z + 22./3.*HPL2_10z + 22./3.*z*HPL2_10z + \
22./3.*HPL2_11z + 22./3.*z*HPL2_11z - 8*HPL3_m1m10z + \
56*z*HPL3_m1m10z + 16*HPL3_m100z - 40*z*HPL3_m100z + 8*HPL3_m101z - \
8*z*HPL3_m101z + 16*HPL3_0m10z + 12*z*HPL3_000z + 8*z*HPL3_001z + \
(-28*MATH::ZETA3 + MATH::ZETA2*(20*HPL1_m1z + 24*z*HPL1_m1z + \
36*std::pow(z,2)*HPL1_m1z) + 32*HPL3_m1m10z - 40*HPL3_m100z - \
16*HPL3_m101z - 24*HPL3_0m10z + 12*HPL3_000z + 8*HPL3_001z)/(1 + z) + \
4*HPL3_100z + 28*z*HPL3_100z + 4*HPL3_101z + 4*z*HPL3_101z - \
4*HPL3_110z - 4*z*HPL3_110z + (36*MATH::ZETA3 + 367./9.*HPL1_1z + \
110./3.*HPL2_00z + 88./3.*HPL2_01z + 44./3.*HPL2_10z + \
44./3.*HPL2_11z + 24*HPL3_0m10z + 12*HPL3_000z + 8*HPL3_001z + \
16*HPL3_100z + 8*HPL3_101z - 8*HPL3_110z + MATH::ZETA2*(-44./3. - \
24*HPL1_1z - 8*L1z) + 367./9.*L1z + -22./3.*std::pow(L1z,2))/(-1 + z)) \
+ 
QCD::CF*
QCD::NF*(-158./27. + -488./27.*z + (8./3. + 8./3.*z)*MATH::ZETA2 - \
(4*HPL1_0z)/(-1 + std::pow(z,2)) + \
(-26./3.*std::pow(z,2)*HPL1_0z)/(-1 + std::pow(z,2)) + \
(-38./3.*std::pow(z,3)*HPL1_0z)/(-1 + std::pow(z,2)) + \
-32./9.*HPL1_1z + -68./9.*z*HPL1_1z + -10./3.*HPL2_00z + \
-10./3.*z*HPL2_00z + -8./3.*HPL2_01z + -8./3.*z*HPL2_01z + \
-4./3.*HPL2_10z + -4./3.*z*HPL2_10z + -4./3.*HPL2_11z + \
-4./3.*z*HPL2_11z + (8./3.*MATH::ZETA2 + -58./9.*HPL1_1z + \
-20./3.*HPL2_00z + -16./3.*HPL2_01z + -8./3.*HPL2_10z + \
-8./3.*HPL2_11z + -58./9.*L1z + 4./3.*std::pow(L1z,2))/(-1 + z)) + \
std::pow(
QCD::CF,2)*(-124./5. + 16./5./z + -461./5.*z + 144./5.*std::pow(z,2) + \
(-64 + 72*z)*MATH::ZETA3 + (-93./5.*HPL1_0z)/(-1 + std::pow(z,2)) + \
(16./5.*HPL1_0z)/(z*(-1 + std::pow(z,2))) + (101./5.*z*HPL1_0z)/(-1 + \
std::pow(z,2)) + (-276./5.*std::pow(z,2)*HPL1_0z)/(-1 + std::pow(z,2)) \
+ (-502./5.*std::pow(z,3)*HPL1_0z)/(-1 + std::pow(z,2)) + \
(144./5.*std::pow(z,4)*HPL1_0z)/(-1 + std::pow(z,2)) + 16*HPL1_1z - \
68*z*HPL1_1z + MATH::ZETA2*(-32 - 8*z + -144./5.*std::pow(z,3) - \
(24*HPL1_0z)/(-1 + std::pow(z,2)) - (8*z*HPL1_0z)/(-1 + std::pow(z,2)) \
- (40*std::pow(z,2)*HPL1_0z)/(-1 + std::pow(z,2)) - \
(24*std::pow(z,3)*HPL1_0z)/(-1 + std::pow(z,2)) - 16*HPL1_1z + \
32*z*HPL1_1z) + 72*HPL2_m10z + (16./5.*HPL2_m10z)/std::pow(z,2) + \
40*z*HPL2_m10z + -144./5.*std::pow(z,3)*HPL2_m10z + 24*HPL2_00z - \
4*z*HPL2_00z + 144./5.*std::pow(z,3)*HPL2_00z + 32*HPL2_01z + \
48*z*HPL2_01z + 32*HPL2_10z + 32*z*HPL2_10z + 28*HPL2_11z + \
36*z*HPL2_11z + 16*HPL3_m1m10z - 112*z*HPL3_m1m10z - 32*HPL3_m100z + \
80*z*HPL3_m100z - 16*HPL3_m101z + 16*z*HPL3_m101z - 32*HPL3_0m10z + \
30*HPL3_000z + 6*z*HPL3_000z + (56*MATH::ZETA3 + \
MATH::ZETA2*(-40*HPL1_m1z - 48*z*HPL1_m1z - 72*std::pow(z,2)*HPL1_m1z) \
- 64*HPL3_m1m10z + 80*HPL3_m100z + 32*HPL3_m101z + 48*HPL3_0m10z - \
24*HPL3_000z - 16*HPL3_001z)/(1 + z) + 40*HPL3_001z + 24*z*HPL3_001z + \
28*HPL3_010z + 28*z*HPL3_010z + 32*HPL3_011z + 32*z*HPL3_011z + \
20*HPL3_100z - 28*z*HPL3_100z + 24*HPL3_101z + 24*z*HPL3_101z + \
32*HPL3_110z + 32*z*HPL3_110z + 24*HPL3_111z + 24*z*HPL3_111z + \
(-72*MATH::ZETA3 - 27*HPL1_1z + 6*HPL2_00z + 24*HPL2_01z + \
36*HPL2_10z + 36*HPL2_11z - 48*HPL3_0m10z + 16*HPL3_000z + \
48*HPL3_001z + 48*HPL3_010z + 56*HPL3_011z + 24*HPL3_100z + \
48*HPL3_101z + 64*HPL3_110z + 48*HPL3_111z + MATH::ZETA2*(12 - \
16*HPL1_1z - 32*L1z) - 27*L1z - 18*std::pow(L1z,2) + \
8*std::pow(L1z,3))/(-1 + z));
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

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.10
double c2g_2_0_reg_approx2(double z)	{
	double result(0.0);
	double L0 = std::log(z);
	double L1 = std::log(1-z);
	result +=	( 58./9. * std::pow(L1,3)
				- 24 * std::pow(L1,2)
				- 34.88 * L1
				+ 30.586
				- ( 25.08 + 760.3 * z + 29.65 * std::pow(L1,3) ) * (1-z)
				+ 1204 * z * std::pow(L0,2)
				+ L0 * L1 * ( 293.8 + 711.2 * z + 1043 * L0 )
				+ 115.6 * L0 
				- 7.109 * std::pow(L0,2)
				+ 70./9. * std::pow(L0,3)
				+ 11.9033 * (1. - z)/z);
	return result;
}

/// @brief eq. (B.6) in [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005]
/// This expression was obtained after some rewriting and reordering with Mathematica.
double c2g_2_0_reg_exact(double z)	{
	double result(0.0);
	double HPL1_m1z	= HPL1(-1,z);
	double HPL1_0z	= HPL1(0,z);
	double HPL1_1z	= HPL1(1,z);
	double HPL2_m10z	= HPL2(-1,0,z);
	double HPL2_00z	= HPL2(0,0,z);
	double HPL2_01z	= HPL2(0,1,z);
	double HPL2_10z	= HPL2(1,0,z);
	double HPL2_11z	= HPL2(1,1,z);
	double HPL3_m1m10z	= HPL3(-1,-1,0,z);
	double HPL3_m100z	= HPL3(-1,0,0,z);
	double HPL3_m101z	= HPL3(-1,0,1,z);
	double HPL3_0m10z	= HPL3(0,-1,0,z);
	double HPL3_000z	= HPL3(0,0,0,z);
	double HPL3_001z	= HPL3(0,0,1,z);
	double HPL3_010z	= HPL3(0,1,0,z);
	double HPL3_011z	= HPL3(0,1,1,z);
	double HPL3_100z	= HPL3(1,0,0,z);
	double HPL3_101z	= HPL3(1,0,1,z);
	double HPL3_110z	= HPL3(1,1,0,z);
	double HPL3_111z	= HPL3(1,1,1,z);
	result += QCD::CF*(-647./15. + 8./15./z + 239./5.*z + -36./5.*std::pow(z,2) + \
(32 + 72*std::pow(z,2))*MATH::ZETA3 + -236./15.*HPL1_0z + \
(-8./15.*HPL1_0z)/z + 113./5.*z*HPL1_0z + \
-216./5.*std::pow(z,2)*HPL1_0z - 14*HPL1_1z + 40*z*HPL1_1z - \
24*std::pow(z,2)*HPL1_1z + MATH::ZETA2*(16 + -104./3.*z + \
72*std::pow(z,2) + 96./5.*std::pow(z,3) - 16*HPL1_m1z - 32*z*HPL1_m1z \
- 16*std::pow(z,2)*HPL1_m1z + 16*HPL1_0z - 32*z*HPL1_0z + \
48*std::pow(z,2)*HPL1_0z + 8*HPL1_1z - 16*z*HPL1_1z + \
32*std::pow(z,2)*HPL1_1z) + 48*HPL2_m10z + \
(8./15.*HPL2_m10z)/std::pow(z,2) + 64./3.*z*HPL2_m10z + \
96./5.*std::pow(z,3)*HPL2_m10z - 3*HPL2_00z + 44./3.*z*HPL2_00z - \
72*std::pow(z,2)*HPL2_00z + -96./5.*std::pow(z,3)*HPL2_00z - \
16*HPL2_01z + 56*z*HPL2_01z - 72*std::pow(z,2)*HPL2_01z - 26*HPL2_10z \
+ 80*z*HPL2_10z - 72*std::pow(z,2)*HPL2_10z - 26*HPL2_11z + \
80*z*HPL2_11z - 72*std::pow(z,2)*HPL2_11z - 32*HPL3_m1m10z - \
64*z*HPL3_m1m10z - 32*std::pow(z,2)*HPL3_m1m10z + 16*HPL3_m100z + \
32*z*HPL3_m100z + 16*std::pow(z,2)*HPL3_m100z + 32*HPL3_0m10z + \
32*std::pow(z,2)*HPL3_0m10z - 10*HPL3_000z + 20*z*HPL3_000z - \
40*std::pow(z,2)*HPL3_000z - 16*HPL3_001z + 32*z*HPL3_001z - \
48*std::pow(z,2)*HPL3_001z - 12*HPL3_010z + 24*z*HPL3_010z - \
32*std::pow(z,2)*HPL3_010z - 16*HPL3_011z + 32*z*HPL3_011z - \
40*std::pow(z,2)*HPL3_011z - 4*HPL3_100z + 8*z*HPL3_100z - \
24*std::pow(z,2)*HPL3_100z - 24*HPL3_101z + 48*z*HPL3_101z - \
48*std::pow(z,2)*HPL3_101z - 16*HPL3_110z + 32*z*HPL3_110z - \
32*std::pow(z,2)*HPL3_110z - 20*HPL3_111z + 40*z*HPL3_111z - \
40*std::pow(z,2)*HPL3_111z) + 
QCD::CA*(239./9. + 344./27./z + 1072./9.*z + -4493./27.*std::pow(z,2) \
+ (4 - 48*z + 24*std::pow(z,2))*MATH::ZETA3 + 58*HPL1_0z + \
584./3.*z*HPL1_0z + -2090./9.*std::pow(z,2)*HPL1_0z + 62./3.*HPL1_1z + \
(-104./9.*HPL1_1z)/z + 454./3.*z*HPL1_1z + \
-1570./9.*std::pow(z,2)*HPL1_1z + MATH::ZETA2*(8 + -16./3./z - 144*z + \
148*std::pow(z,2) - 4*HPL1_m1z - 8*z*HPL1_m1z - \
16*std::pow(z,2)*HPL1_m1z - 8*HPL1_0z - 64*z*HPL1_0z + \
16*std::pow(z,2)*HPL1_0z + 8*HPL1_1z - 16*z*HPL1_1z + \
8*std::pow(z,2)*HPL1_1z) - 24*HPL2_m10z + (-16./3.*HPL2_m10z)/z + \
80./3.*std::pow(z,2)*HPL2_m10z - 2*HPL2_00z + 176*z*HPL2_00z + \
-388./3.*std::pow(z,2)*HPL2_00z - 8*HPL2_01z + 144*z*HPL2_01z - \
148*std::pow(z,2)*HPL2_01z - 4*HPL2_10z + (16./3.*HPL2_10z)/z + \
80*z*HPL2_10z + -268./3.*std::pow(z,2)*HPL2_10z - 4*HPL2_11z + \
(16./3.*HPL2_11z)/z + 72*z*HPL2_11z + -244./3.*std::pow(z,2)*HPL2_11z \
+ 8*HPL3_m1m10z + 16*z*HPL3_m1m10z + 8*HPL3_m100z + 16*z*HPL3_m100z + \
24*std::pow(z,2)*HPL3_m100z + 8*HPL3_m101z + 16*z*HPL3_m101z + \
16*std::pow(z,2)*HPL3_m101z + 16*std::pow(z,2)*HPL3_0m10z + \
20*HPL3_000z + 56*z*HPL3_000z + 8*HPL3_001z + 64*z*HPL3_001z - \
16*std::pow(z,2)*HPL3_001z + 48*z*HPL3_010z - \
16*std::pow(z,2)*HPL3_010z + 48*z*HPL3_011z - \
16*std::pow(z,2)*HPL3_011z - 12*HPL3_100z + 24*z*HPL3_100z - \
16*std::pow(z,2)*HPL3_100z - 4*HPL3_101z + 8*z*HPL3_101z - \
8*std::pow(z,2)*HPL3_101z - 12*HPL3_110z + 24*z*HPL3_110z - \
24*std::pow(z,2)*HPL3_110z - 4*HPL3_111z + 8*z*HPL3_111z - \
8*std::pow(z,2)*HPL3_111z);
	return result;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 588 (1), 345, 2000] eq. 3.6
double c2g_2_0_local_approx()	{
	return -0.28;
}

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.10
double c2g_2_0_local_approx2()	{
	return 0;
}

/// @brief eq. (B.7) in [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005]
/// (just put here to stick to the structure "approx <--> exact")
double c2g_2_0_local_exact()	{
	return 0;
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

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.9
double c2q_ps_2_0_reg_approx2(double z)	{
	double result(0.0);
	double L0 = std::log(z);
	double L1 = std::log(1-z);
	result +=	( 8./3. * std::pow(L1,2) - 32./3. * L1 + 9.8937 ) * (1-z)
				+ (9.57 - 13.41 * z + 0.08 * std::pow(L1,3)) * std::pow(1-z,2)
				+ 5.667 * z * std::pow(L0,3)
				- std::pow(L0,2) * L1 * (20.26 - 33.93 * z)
				+ 43.36 * (1-z) * L0
				- 1.053 * std::pow(L0,2)
				+ 40./9. * std::pow(L0,3)
				+ 5.2903 * std::pow(1-z,2)/z;
	return result;
}

/// @brief eq. (B.7) in [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005]
/// This expression was obtained after some rewriting and reordering with Mathematica.
double c2q_ps_2_0_reg_exact(double z)	{
	double result(0.0);
	double HPL1_0z	= HPL1(0,z);
	double HPL1_1z	= HPL1(1,z);
	double HPL2_m10z	= HPL2(-1,0,z);
	double HPL2_00z	= HPL2(0,0,z);
	double HPL2_01z	= HPL2(0,1,z);
	double HPL2_10z	= HPL2(1,0,z);
	double HPL2_11z	= HPL2(1,1,z);
	double HPL3_000z	= HPL3(0,0,0,z);
	double HPL3_001z	= HPL3(0,0,1,z);
	double HPL3_010z	= HPL3(0,1,0,z);
	double HPL3_011z	= HPL3(0,1,1,z);
	result += QCD::CF*(158./9. + 344./27./z + -422./9.*z + 448./27.*std::pow(z,2) + \
(-8 - 8*z)*MATH::ZETA3 + 56*HPL1_0z + -88./3.*z*HPL1_0z + \
-128./9.*std::pow(z,2)*HPL1_0z + MATH::ZETA2*(-16./3./z - 16*z + \
16*std::pow(z,2) - 16*HPL1_0z - 16*z*HPL1_0z) + 104./3.*HPL1_1z + \
(-104./9.*HPL1_1z)/z + -80./3.*z*HPL1_1z + \
32./9.*std::pow(z,2)*HPL1_1z - 16*HPL2_m10z + (-16./3.*HPL2_m10z)/z - \
16*z*HPL2_m10z + -16./3.*std::pow(z,2)*HPL2_m10z - 2*HPL2_00z + \
30*z*HPL2_00z + -64./3.*std::pow(z,2)*HPL2_00z - \
16*std::pow(z,2)*HPL2_01z + 4*HPL2_10z + (16./3.*HPL2_10z)/z - \
4*z*HPL2_10z + -16./3.*std::pow(z,2)*HPL2_10z + 4*HPL2_11z + \
(16./3.*HPL2_11z)/z - 4*z*HPL2_11z + -16./3.*std::pow(z,2)*HPL2_11z + \
20*HPL3_000z + 20*z*HPL3_000z + 16*HPL3_001z + 16*z*HPL3_001z + \
8*HPL3_010z + 8*z*HPL3_010z + 8*HPL3_011z + 8*z*HPL3_011z);
	return result;
}