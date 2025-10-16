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
	double HPL1_1 	= HPL1(1,x);
	double HPL2_10 	= HPL2(1,0,x);
	double HPL2_11 	= HPL2(1,1,x);
	double HPL3_100 	= HPL3(1,0,0,x);
	double HPL3_101 	= HPL3(1,0,1,x);
	double HPL3_110 	= HPL3(1,1,0,x);
	double HPL3_111 	= HPL3(1,1,1,x);
	double HPL4_10m10 	= HPL4(1,0,-1,0,x);
	double HPL4_1000 	= HPL4(1,0,0,0,x);
	double HPL4_1001 	= HPL4(1,0,0,1,x);
	double HPL4_1010 	= HPL4(1,0,1,0,x);
	double HPL4_1011 	= HPL4(1,0,1,1,x);
	double HPL4_1100 	= HPL4(1,1,0,0,x);
	double HPL4_1101 	= HPL4(1,1,0,1,x);
	double HPL4_1110 	= HPL4(1,1,1,0,x);
	double HPL4_1111 	= HPL4(1,1,1,1,x);
	result += QCD::CF*
QCD::NF*(-247./27.*HPL1_1 + 16./3.*MATH::ZETA2*HPL1_1 + \
-38./3.*HPL2_10 + -58./9.*HPL2_11 + -20./3.*HPL3_100 + \
-16./3.*HPL3_101 + -8./3.*HPL3_110 + -8./3.*HPL3_111) + 
QCD::CA*
QCD::CF*(3155./54.*HPL1_1 - 4*MATH::ZETA3*HPL1_1 + 239./3.*HPL2_10 + \
MATH::ZETA2*(-88./3.*HPL1_1 - 8*HPL2_10 - 24*HPL2_11) + \
367./9.*HPL2_11 + 110./3.*HPL3_100 + 88./3.*HPL3_101 + \
44./3.*HPL3_110 + 44./3.*HPL3_111 + 24*HPL4_10m10 + 12*HPL4_1000 + \
8*HPL4_1001 + 16*HPL4_1100 + 8*HPL4_1101 - 8*HPL4_1110) + std::pow(
QCD::CF,2)*(-51./2.*HPL1_1 - 64*MATH::ZETA3*HPL1_1 - 61*HPL2_10 + \
MATH::ZETA2*(-24*HPL1_1 - 48*HPL2_10 - 16*HPL2_11) - 27*HPL2_11 + \
6*HPL3_100 + 24*HPL3_101 + 36*HPL3_110 + 36*HPL3_111 - 48*HPL4_10m10 + \
16*HPL4_1000 + 48*HPL4_1001 + 48*HPL4_1010 + 56*HPL4_1011 + \
24*HPL4_1100 + 48*HPL4_1101 + 64*HPL4_1110 + 48*HPL4_1111);
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
	double HPL1_0	= HPL1(0,z);
	double HPL1_1	= HPL1(1,z);
	double HPL2_00	= HPL2(0,0,z);
	double HPL2_01	= HPL2(0,1,z);
	double HPL2_10	= HPL2(1,0,z);
	double HPL2_11	= HPL2(1,1,z);
	double HPL3_0m10	= HPL3(0,-1,0,z);
	double HPL3_000	= HPL3(0,0,0,z);
	double HPL3_001	= HPL3(0,0,1,z);
	double HPL3_010	= HPL3(0,1,0,z);
	double HPL3_011	= HPL3(0,1,1,z);
	double HPL3_100	= HPL3(1,0,0,z);
	double HPL3_101	= HPL3(1,0,1,z);
	double HPL3_110	= HPL3(1,1,0,z);
	double HPL3_111	= HPL3(1,1,1,z);
	result += (
QCD::CF*
QCD::NF*(247./27. + -16./3.*MATH::ZETA2 + 38./3.*HPL1_0 + \
58./9.*HPL1_1 + 20./3.*HPL2_00 + 16./3.*HPL2_01 + 8./3.*HPL2_10 + \
8./3.*HPL2_11) + 
QCD::CA*
QCD::CF*(-3155./54. + 4*MATH::ZETA3 + -239./3.*HPL1_0 + \
-367./9.*HPL1_1 + MATH::ZETA2*(88./3. + 8*HPL1_0 + 24*HPL1_1) + \
-110./3.*HPL2_00 + -88./3.*HPL2_01 + -44./3.*HPL2_10 + \
-44./3.*HPL2_11 - 24*HPL3_0m10 - 12*HPL3_000 - 8*HPL3_001 - \
16*HPL3_100 - 8*HPL3_101 + 8*HPL3_110) + std::pow(
QCD::CF,2)*(51./2. + 64*MATH::ZETA3 + 61*HPL1_0 + 27*HPL1_1 + \
MATH::ZETA2*(24 + 48*HPL1_0 + 16*HPL1_1) - 6*HPL2_00 - 24*HPL2_01 - \
36*HPL2_10 - 36*HPL2_11 + 48*HPL3_0m10 - 16*HPL3_000 - 48*HPL3_001 - \
48*HPL3_010 - 56*HPL3_011 - 24*HPL3_100 - 48*HPL3_101 - 64*HPL3_110 - \
48*HPL3_111))/(1 - z);
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
	double HPL1_m1	= HPL1(-1,z);
	double HPL1_0	= HPL1(0,z);
	double HPL1_1	= HPL1(1,z);
	double HPL2_m10	= HPL2(-1,0,z);
	double HPL2_00	= HPL2(0,0,z);
	double HPL2_01	= HPL2(0,1,z);
	double HPL2_10	= HPL2(1,0,z);
	double HPL2_11	= HPL2(1,1,z);
	double HPL3_m1m10	= HPL3(-1,-1,0,z);
	double HPL3_m100	= HPL3(-1,0,0,z);
	double HPL3_m101	= HPL3(-1,0,1,z);
	double HPL3_0m10	= HPL3(0,-1,0,z);
	double HPL3_000	= HPL3(0,0,0,z);
	double HPL3_001	= HPL3(0,0,1,z);
	double HPL3_010	= HPL3(0,1,0,z);
	double HPL3_011	= HPL3(0,1,1,z);
	double HPL3_100	= HPL3(1,0,0,z);
	double HPL3_101	= HPL3(1,0,1,z);
	double HPL3_110	= HPL3(1,1,0,z);
	double HPL3_111	= HPL3(1,1,1,z);
	result += QCD::CF*
QCD::NF*(-158./27. + -488./27.*z + (8./3. + 8./3.*z)*MATH::ZETA2 + \
-26./3.*HPL1_0 + -38./3.*z*HPL1_0 + -32./9.*HPL1_1 + -68./9.*z*HPL1_1 \
+ -10./3.*HPL2_00 + -10./3.*z*HPL2_00 + -8./3.*HPL2_01 + \
-8./3.*z*HPL2_01 + -4./3.*HPL2_10 + -4./3.*z*HPL2_10 + -4./3.*HPL2_11 \
+ -4./3.*z*HPL2_11) + 
QCD::CA*
QCD::CF*(3709./135. + -8./5./z + 17626./135.*z + -72./5.*std::pow(z,2) \
+ (12 - 56*z - 28/(1 + z))*MATH::ZETA3 + 583./15.*HPL1_0 + \
(8./5.*HPL1_0)/z + 1693./15.*z*HPL1_0 + -72./5.*std::pow(z,2)*HPL1_0 + \
(8*HPL1_0)/(1 + z) + 56./9.*HPL1_1 + 668./9.*z*HPL1_1 + \
MATH::ZETA2*(-44./3. + -104./3.*z + 72./5.*std::pow(z,3) - 12*HPL1_m1 \
+ 36*z*HPL1_m1 + (32*HPL1_m1)/(1 + z) - 8*z*HPL1_0 - (8*HPL1_0)/(1 + \
z) - 8*HPL1_1 - 32*z*HPL1_1) - 36*HPL2_m10 + \
(-8./5.*HPL2_m10)/std::pow(z,2) - 20*z*HPL2_m10 + \
72./5.*std::pow(z,3)*HPL2_m10 + 55./3.*HPL2_00 + 115./3.*z*HPL2_00 + \
-72./5.*std::pow(z,3)*HPL2_00 + 44./3.*HPL2_01 + 44./3.*z*HPL2_01 + \
22./3.*HPL2_10 + 22./3.*z*HPL2_10 + 22./3.*HPL2_11 + 22./3.*z*HPL2_11 \
- 8*HPL3_m1m10 + 56*z*HPL3_m1m10 + (32*HPL3_m1m10)/(1 + z) + \
16*HPL3_m100 - 40*z*HPL3_m100 - (40*HPL3_m100)/(1 + z) + 8*HPL3_m101 - \
8*z*HPL3_m101 - (16*HPL3_m101)/(1 + z) + 16*HPL3_0m10 - \
(24*HPL3_0m10)/(1 + z) + 12*z*HPL3_000 + (12*HPL3_000)/(1 + z) + \
8*z*HPL3_001 + (8*HPL3_001)/(1 + z) + 4*HPL3_100 + 28*z*HPL3_100 + \
4*HPL3_101 + 4*z*HPL3_101 - 4*HPL3_110 - 4*z*HPL3_110) + std::pow(
QCD::CF,2)*(-124./5. + 16./5./z + -461./5.*z + 144./5.*std::pow(z,2) + \
(-64 + 72*z + 56/(1 + z))*MATH::ZETA3 + -132./5.*HPL1_0 + \
(-16./5.*HPL1_0)/z + -502./5.*z*HPL1_0 + 144./5.*std::pow(z,2)*HPL1_0 \
- (16*HPL1_0)/(1 + z) + 16*HPL1_1 - 68*z*HPL1_1 + MATH::ZETA2*(-32 - \
8*z + -144./5.*std::pow(z,3) + 24*HPL1_m1 - 72*z*HPL1_m1 - \
(64*HPL1_m1)/(1 + z) - 40*HPL1_0 - 24*z*HPL1_0 + (16*HPL1_0)/(1 + z) - \
16*HPL1_1 + 32*z*HPL1_1) + 72*HPL2_m10 + \
(16./5.*HPL2_m10)/std::pow(z,2) + 40*z*HPL2_m10 + \
-144./5.*std::pow(z,3)*HPL2_m10 + 24*HPL2_00 - 4*z*HPL2_00 + \
144./5.*std::pow(z,3)*HPL2_00 + 32*HPL2_01 + 48*z*HPL2_01 + \
32*HPL2_10 + 32*z*HPL2_10 + 28*HPL2_11 + 36*z*HPL2_11 + 16*HPL3_m1m10 \
- 112*z*HPL3_m1m10 - (64*HPL3_m1m10)/(1 + z) - 32*HPL3_m100 + \
80*z*HPL3_m100 + (80*HPL3_m100)/(1 + z) - 16*HPL3_m101 + \
16*z*HPL3_m101 + (32*HPL3_m101)/(1 + z) - 32*HPL3_0m10 + \
(48*HPL3_0m10)/(1 + z) + 30*HPL3_000 + 6*z*HPL3_000 - \
(24*HPL3_000)/(1 + z) + 40*HPL3_001 + 24*z*HPL3_001 - \
(16*HPL3_001)/(1 + z) + 28*HPL3_010 + 28*z*HPL3_010 + 32*HPL3_011 + \
32*z*HPL3_011 + 20*HPL3_100 - 28*z*HPL3_100 + 24*HPL3_101 + \
24*z*HPL3_101 + 32*HPL3_110 + 32*z*HPL3_110 + 24*HPL3_111 + \
24*z*HPL3_111);
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
	double HPL1_m1	= HPL1(-1,z);
	double HPL1_0	= HPL1(0,z);
	double HPL1_1	= HPL1(1,z);
	double HPL2_m10	= HPL2(-1,0,z);
	double HPL2_00	= HPL2(0,0,z);
	double HPL2_01	= HPL2(0,1,z);
	double HPL2_10	= HPL2(1,0,z);
	double HPL2_11	= HPL2(1,1,z);
	double HPL3_m1m10	= HPL3(-1,-1,0,z);
	double HPL3_m100	= HPL3(-1,0,0,z);
	double HPL3_m101	= HPL3(-1,0,1,z);
	double HPL3_0m10	= HPL3(0,-1,0,z);
	double HPL3_000	= HPL3(0,0,0,z);
	double HPL3_001	= HPL3(0,0,1,z);
	double HPL3_010	= HPL3(0,1,0,z);
	double HPL3_011	= HPL3(0,1,1,z);
	double HPL3_100	= HPL3(1,0,0,z);
	double HPL3_101	= HPL3(1,0,1,z);
	double HPL3_110	= HPL3(1,1,0,z);
	double HPL3_111	= HPL3(1,1,1,z);
	result += QCD::CF*(-647./15. + 8./15./z + 239./5.*z + -36./5.*std::pow(z,2) + \
(32 + 72*std::pow(z,2))*MATH::ZETA3 + -236./15.*HPL1_0 + \
(-8./15.*HPL1_0)/z + 113./5.*z*HPL1_0 + -216./5.*std::pow(z,2)*HPL1_0 \
- 14*HPL1_1 + 40*z*HPL1_1 - 24*std::pow(z,2)*HPL1_1 + MATH::ZETA2*(16 \
+ -104./3.*z + 72*std::pow(z,2) + 96./5.*std::pow(z,3) - 16*HPL1_m1 - \
32*z*HPL1_m1 - 16*std::pow(z,2)*HPL1_m1 + 16*HPL1_0 - 32*z*HPL1_0 + \
48*std::pow(z,2)*HPL1_0 + 8*HPL1_1 - 16*z*HPL1_1 + \
32*std::pow(z,2)*HPL1_1) + 48*HPL2_m10 + \
(8./15.*HPL2_m10)/std::pow(z,2) + 64./3.*z*HPL2_m10 + \
96./5.*std::pow(z,3)*HPL2_m10 - 3*HPL2_00 + 44./3.*z*HPL2_00 - \
72*std::pow(z,2)*HPL2_00 + -96./5.*std::pow(z,3)*HPL2_00 - 16*HPL2_01 \
+ 56*z*HPL2_01 - 72*std::pow(z,2)*HPL2_01 - 26*HPL2_10 + 80*z*HPL2_10 \
- 72*std::pow(z,2)*HPL2_10 - 26*HPL2_11 + 80*z*HPL2_11 - \
72*std::pow(z,2)*HPL2_11 - 32*HPL3_m1m10 - 64*z*HPL3_m1m10 - \
32*std::pow(z,2)*HPL3_m1m10 + 16*HPL3_m100 + 32*z*HPL3_m100 + \
16*std::pow(z,2)*HPL3_m100 + 32*HPL3_0m10 + \
32*std::pow(z,2)*HPL3_0m10 - 10*HPL3_000 + 20*z*HPL3_000 - \
40*std::pow(z,2)*HPL3_000 - 16*HPL3_001 + 32*z*HPL3_001 - \
48*std::pow(z,2)*HPL3_001 - 12*HPL3_010 + 24*z*HPL3_010 - \
32*std::pow(z,2)*HPL3_010 - 16*HPL3_011 + 32*z*HPL3_011 - \
40*std::pow(z,2)*HPL3_011 - 4*HPL3_100 + 8*z*HPL3_100 - \
24*std::pow(z,2)*HPL3_100 - 24*HPL3_101 + 48*z*HPL3_101 - \
48*std::pow(z,2)*HPL3_101 - 16*HPL3_110 + 32*z*HPL3_110 - \
32*std::pow(z,2)*HPL3_110 - 20*HPL3_111 + 40*z*HPL3_111 - \
40*std::pow(z,2)*HPL3_111) + 
QCD::CA*(239./9. + 344./27./z + 1072./9.*z + -4493./27.*std::pow(z,2) \
+ (4 - 48*z + 24*std::pow(z,2))*MATH::ZETA3 + 58*HPL1_0 + \
584./3.*z*HPL1_0 + -2090./9.*std::pow(z,2)*HPL1_0 + 62./3.*HPL1_1 + \
(-104./9.*HPL1_1)/z + 454./3.*z*HPL1_1 + \
-1570./9.*std::pow(z,2)*HPL1_1 + MATH::ZETA2*(8 + -16./3./z - 144*z + \
148*std::pow(z,2) - 4*HPL1_m1 - 8*z*HPL1_m1 - \
16*std::pow(z,2)*HPL1_m1 - 8*HPL1_0 - 64*z*HPL1_0 + \
16*std::pow(z,2)*HPL1_0 + 8*HPL1_1 - 16*z*HPL1_1 + \
8*std::pow(z,2)*HPL1_1) - 24*HPL2_m10 + (-16./3.*HPL2_m10)/z + \
80./3.*std::pow(z,2)*HPL2_m10 - 2*HPL2_00 + 176*z*HPL2_00 + \
-388./3.*std::pow(z,2)*HPL2_00 - 8*HPL2_01 + 144*z*HPL2_01 - \
148*std::pow(z,2)*HPL2_01 - 4*HPL2_10 + (16./3.*HPL2_10)/z + \
80*z*HPL2_10 + -268./3.*std::pow(z,2)*HPL2_10 - 4*HPL2_11 + \
(16./3.*HPL2_11)/z + 72*z*HPL2_11 + -244./3.*std::pow(z,2)*HPL2_11 + \
8*HPL3_m1m10 + 16*z*HPL3_m1m10 + 8*HPL3_m100 + 16*z*HPL3_m100 + \
24*std::pow(z,2)*HPL3_m100 + 8*HPL3_m101 + 16*z*HPL3_m101 + \
16*std::pow(z,2)*HPL3_m101 + 16*std::pow(z,2)*HPL3_0m10 + 20*HPL3_000 \
+ 56*z*HPL3_000 + 8*HPL3_001 + 64*z*HPL3_001 - \
16*std::pow(z,2)*HPL3_001 + 48*z*HPL3_010 - 16*std::pow(z,2)*HPL3_010 \
+ 48*z*HPL3_011 - 16*std::pow(z,2)*HPL3_011 - 12*HPL3_100 + \
24*z*HPL3_100 - 16*std::pow(z,2)*HPL3_100 - 4*HPL3_101 + 8*z*HPL3_101 \
- 8*std::pow(z,2)*HPL3_101 - 12*HPL3_110 + 24*z*HPL3_110 - \
24*std::pow(z,2)*HPL3_110 - 4*HPL3_111 + 8*z*HPL3_111 - \
8*std::pow(z,2)*HPL3_111);
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
	double HPL1_0	= HPL1(0,z);
	double HPL1_1	= HPL1(1,z);
	double HPL2_m10	= HPL2(-1,0,z);
	double HPL2_00	= HPL2(0,0,z);
	double HPL2_01	= HPL2(0,1,z);
	double HPL2_10	= HPL2(1,0,z);
	double HPL2_11	= HPL2(1,1,z);
	double HPL3_000	= HPL3(0,0,0,z);
	double HPL3_001	= HPL3(0,0,1,z);
	double HPL3_010	= HPL3(0,1,0,z);
	double HPL3_011	= HPL3(0,1,1,z);
	result += QCD::CF*(158./9. + 344./27./z + -422./9.*z + 448./27.*std::pow(z,2) + \
(-8 - 8*z)*MATH::ZETA3 + 56*HPL1_0 + -88./3.*z*HPL1_0 + \
-128./9.*std::pow(z,2)*HPL1_0 + MATH::ZETA2*(-16./3./z - 16*z + \
16*std::pow(z,2) - 16*HPL1_0 - 16*z*HPL1_0) + 104./3.*HPL1_1 + \
(-104./9.*HPL1_1)/z + -80./3.*z*HPL1_1 + 32./9.*std::pow(z,2)*HPL1_1 - \
16*HPL2_m10 + (-16./3.*HPL2_m10)/z - 16*z*HPL2_m10 + \
-16./3.*std::pow(z,2)*HPL2_m10 - 2*HPL2_00 + 30*z*HPL2_00 + \
-64./3.*std::pow(z,2)*HPL2_00 - 16*std::pow(z,2)*HPL2_01 + 4*HPL2_10 + \
(16./3.*HPL2_10)/z - 4*z*HPL2_10 + -16./3.*std::pow(z,2)*HPL2_10 + \
4*HPL2_11 + (16./3.*HPL2_11)/z - 4*z*HPL2_11 + \
-16./3.*std::pow(z,2)*HPL2_11 + 20*HPL3_000 + 20*z*HPL3_000 + \
16*HPL3_001 + 16*z*HPL3_001 + 8*HPL3_010 + 8*z*HPL3_010 + 8*HPL3_011 + \
8*z*HPL3_011);
	return result;
}