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

/// @brief eq. (B.5) in [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005]
/// This expression was obtained after some rewriting and reordering with Mathematica.
///
/*
double c2q_ns_2_0_exact(double z)	{
	double result(0.0);
	result += (std::pow(QCD::CF,2)*(-93 - 209*z))/4. + (QCD::CA*QCD::CF*(139 + 3159*z))/36. - (4*std::pow(QCD::CF,2)*(7 + 13*z)*MATH::ZETA2)/5. - 4*std::pow(QCD::CF,2)*(1 - 19*z)*MATH::ZETA3 + 
   (std::pow(QCD::CF,2)*(331./8. + 69*MATH::ZETA2 + 6*std::pow(MATH::ZETA2,2) - 78*MATH::ZETA3) - QCD::CA*QCD::CF*(5465./72. + (251*MATH::ZETA2)/3. - (71*std::pow(MATH::ZETA2,2))/5. - (140*MATH::ZETA3)/3.) + 
      QCD::CF*QCD::NF*(457./36. + (38*MATH::ZETA2)/3. + (4*MATH::ZETA3)/3.))*DiracDelta(-1 + z) + (8*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(-1 + 2*z)*(1 - HPL1(0,z)))/(5.*(-2 + z)) - 
   (std::pow(QCD::CF,2)*(43 + 63*z)*HPL1(0,z))/2. + (QCD::CA*QCD::CF*(71 + 323*z)*HPL1(0,z))/6. + (std::pow(QCD::CF,2)*(59 - 109*z)*HPL1(1,z))/2. - (17*QCD::CA*QCD::CF*(5 - 19*z)*HPL1(1,z))/6. + 
   (296*QCD::CF*(-0.5*QCD::CA + QCD::CF)*HPL2(-1,0,z))/5. + (8*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(1 - 1/z)*(-1 - 2*z)*HPL2(-1,0,z))/(5.*(-2 - z)) + (136*QCD::CF*(-0.5*QCD::CA + QCD::CF)*z*HPL2(-1,0,z))/5. + 
   (72*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(1 - z)*(1 + 2*z + 2*std::pow(z,2))*HPL2(-1,0,z))/5. - 
   (72*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(1 - 2*z + 2*std::pow(z,2))*(-1 - HPL1(0,z) + (1 + z)*(MATH::ZETA2 - HPL2(0,0,z))))/5. - (4*QCD::CA*QCD::CF*(9 + 16*z)*(MATH::ZETA2 - HPL2(0,0,z)))/5. + 
   (std::pow(QCD::CF,2)*(33 + 37*z)*HPL2(0,0,z))/5. + 2*std::pow(QCD::CF,2)*(5 + 9*z)*(2*HPL2(0,1,z) + HPL2(1,1,z)) + 
   QCD::NF*((QCD::CF*(-23 - 243*z))/18. - (QCD::CF*(7 + 19*z)*HPL1(0,z))/3. - (QCD::CF*(1 + 13*z)*HPL1(1,z))/3. + 
      (QCD::CF*(-247 + 144*MATH::ZETA2 - 342*HPL1(0,z) - 174*HPL1(1,z) - 180*HPL2(0,0,z) - 144*HPL2(0,1,z) - 72*HPL2(1,0,z) - 72*HPL2(1,1,z)))/54. - 
      (QCD::CF*z*(247 - 144*MATH::ZETA2 + 342*HPL1(0,z) + 174*HPL1(1,z) + 180*HPL2(0,0,z) + 144*HPL2(0,1,z) + 72*HPL2(1,0,z) + 72*HPL2(1,1,z)))/54.) - 
   8*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(1 + 5*z)*(MATH::ZETA2*HPL1(-1,z) + 2*HPL3(-1,-1,0,z) - HPL3(-1,0,0,z)) + 16*std::pow(QCD::CF,2)*HPL3(0,-1,0,z) - 8*QCD::CA*QCD::CF*(5*z*MATH::ZETA3 + HPL3(0,-1,0,z)) + 
   4*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(-1 + z + 2/(1 + z))*(7*MATH::ZETA3 - 8*MATH::ZETA2*HPL1(-1,z) - 2*HPL1(0,z) + 2*MATH::ZETA2*HPL1(0,z) - 8*HPL3(-1,-1,0,z) + 10*HPL3(-1,0,0,z) + 4*HPL3(-1,0,1,z) + 
      6*HPL3(0,-1,0,z) - 3*HPL3(0,0,0,z) - 2*HPL3(0,0,1,z)) + 2*std::pow(QCD::CF,2)*(1 + z)*
    (-4*MATH::ZETA2*HPL1(0,z) + 7*HPL2(1,0,z) + 5*HPL3(0,0,0,z) + 4*HPL3(0,0,1,z) + 2*HPL3(0,1,0,z) + 2*HPL3(0,1,1,z)) + 
   8*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(1 - 5*z)*(-(MATH::ZETA2*HPL1(1,z)) + HPL3(1,0,0,z)) + 
   (QCD::CA*QCD::CF*(3155 - 1584*MATH::ZETA2 - 216*MATH::ZETA3 + 4302*HPL1(0,z) - 432*MATH::ZETA2*HPL1(0,z) + 2202*HPL1(1,z) - 1296*MATH::ZETA2*HPL1(1,z) + 1980*HPL2(0,0,z) + 1584*HPL2(0,1,z) + 
        792*HPL2(1,0,z) + 792*HPL2(1,1,z) + 1296*HPL3(0,-1,0,z) + 648*HPL3(0,0,0,z) + 432*HPL3(0,0,1,z) + 864*HPL3(1,0,0,z) + 432*HPL3(1,0,1,z) - 432*HPL3(1,1,0,z)))/108.\
    - (QCD::CA*QCD::CF*z*(-3155 + 1584*MATH::ZETA2 + 216*MATH::ZETA3 - 4302*HPL1(0,z) + 432*MATH::ZETA2*HPL1(0,z) - 2202*HPL1(1,z) + 1296*MATH::ZETA2*HPL1(1,z) - 1980*HPL2(0,0,z) - 1584*HPL2(0,1,z) - 
        792*HPL2(1,0,z) - 792*HPL2(1,1,z) - 1296*HPL3(0,-1,0,z) - 648*HPL3(0,0,0,z) - 432*HPL3(0,0,1,z) - 864*HPL3(1,0,0,z) - 432*HPL3(1,0,1,z) + 432*HPL3(1,1,0,z)))/108.\
    + ((QCD::CF*QCD::NF*(247 - 144*MATH::ZETA2 + 342*HPL1(0,z) + 174*HPL1(1,z) + 180*HPL2(0,0,z) + 144*HPL2(0,1,z) + 72*HPL2(1,0,z) + 72*HPL2(1,1,z)))/27. + 
      (QCD::CA*QCD::CF*(-3155 + 1584*MATH::ZETA2 + 216*MATH::ZETA3 - 4302*HPL1(0,z) + 432*MATH::ZETA2*HPL1(0,z) - 2202*HPL1(1,z) + 1296*MATH::ZETA2*HPL1(1,z) - 1980*HPL2(0,0,z) - 1584*HPL2(0,1,z) - 
           792*HPL2(1,0,z) - 792*HPL2(1,1,z) - 1296*HPL3(0,-1,0,z) - 648*HPL3(0,0,0,z) - 432*HPL3(0,0,1,z) - 864*HPL3(1,0,0,z) - 432*HPL3(1,0,1,z) + 432*HPL3(1,1,0,z)))/
       54. + (std::pow(QCD::CF,2)*(51 + 48*MATH::ZETA2 + 128*MATH::ZETA3 + 122*HPL1(0,z) + 96*MATH::ZETA2*HPL1(0,z) + 54*HPL1(1,z) + 32*MATH::ZETA2*HPL1(1,z) - 12*HPL2(0,0,z) - 48*HPL2(0,1,z) - 
           72*HPL2(1,0,z) - 72*HPL2(1,1,z) + 96*HPL3(0,-1,0,z) - 32*HPL3(0,0,0,z) - 96*HPL3(0,0,1,z) - 96*HPL3(0,1,0,z) - 112*HPL3(0,1,1,z) - 48*HPL3(1,0,0,z) - 
           96*HPL3(1,0,1,z) - 128*HPL3(1,1,0,z) - 96*HPL3(1,1,1,z)))/2.)/(1 - z) - 
   (std::pow(QCD::CF,2)*z*(51 + 48*MATH::ZETA2 + 128*MATH::ZETA3 + 122*HPL1(0,z) + 96*MATH::ZETA2*HPL1(0,z) + 54*HPL1(1,z) + 32*MATH::ZETA2*HPL1(1,z) - 12*HPL2(0,0,z) - 48*HPL2(0,1,z) - 72*HPL2(1,0,z) - 
        72*HPL2(1,1,z) + 96*HPL3(0,-1,0,z) - 32*HPL3(0,0,0,z) - 96*HPL3(0,0,1,z) - 96*HPL3(0,1,0,z) - 112*HPL3(0,1,1,z) - 48*HPL3(1,0,0,z) - 96*HPL3(1,0,1,z) - 
        128*HPL3(1,1,0,z) - 96*HPL3(1,1,1,z)))/4. + (std::pow(QCD::CF,2)*(-51 - 48*MATH::ZETA2 - 128*MATH::ZETA3 - 122*HPL1(0,z) - 96*MATH::ZETA2*HPL1(0,z) - 54*HPL1(1,z) - 32*MATH::ZETA2*HPL1(1,z) + 
        12*HPL2(0,0,z) + 48*HPL2(0,1,z) + 72*HPL2(1,0,z) + 72*HPL2(1,1,z) - 96*HPL3(0,-1,0,z) + 32*HPL3(0,0,0,z) + 96*HPL3(0,0,1,z) + 96*HPL3(0,1,0,z) + 
        112*HPL3(0,1,1,z) + 48*HPL3(1,0,0,z) + 96*HPL3(1,0,1,z) + 128*HPL3(1,1,0,z) + 96*HPL3(1,1,1,z)))/4.
	return result;
}
*/

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_0_local_approx()	{	
	return -338.046;
}

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_1_local_approx()	{	
	return 46.8405;
}

/// @brief from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005] eq. 4.8
double c2q_ns_2_0_local_approx1()	{	
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
	result += (std::pow(QCD::CF,2)*(331./8. + 69*MATH::ZETA2 + 6*std::pow(MATH::ZETA2,2) - 78*MATH::ZETA3) - QCD::CA*QCD::CF*(5465./72. + (251*MATH::ZETA2)/3. - (71*std::pow(MATH::ZETA2,2))/5. - (140*MATH::ZETA3)/3.) + 
      QCD::CF*QCD::NF*(457./36. + (38*MATH::ZETA2)/3. + (4*MATH::ZETA3)/3.));
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
	result += (QCD::CF*QCD::NF*(-494 + 288*MATH::ZETA2 - 684*HPL2(1,0,x) - 348*HPL2(1,1,x) - 360*HPL3(1,0,0,x) - 288*HPL3(1,0,1,x) - 144*HPL3(1,1,0,x) - 144*HPL3(1,1,1,x)))/54. + 
   (QCD::CA*QCD::CF*(3155 - 1584*MATH::ZETA2 - 216*MATH::ZETA3 + 4302*HPL2(1,0,x) - 432*MATH::ZETA2*HPL2(1,0,x) - 6*(-367 + 216*MATH::ZETA2)*HPL2(1,1,x) + 1980*HPL3(1,0,0,x) + 1584*HPL3(1,0,1,x) + 792*HPL3(1,1,0,x) + 792*HPL3(1,1,1,x) + 
        1296*HPL4(1,0,-1,0,x) + 648*HPL4(1,0,0,0,x) + 432*HPL4(1,0,0,1,x) + 864*HPL4(1,1,0,0,x) + 432*HPL4(1,1,0,1,x) - 432*HPL4(1,1,1,0,x)))/54. + 
   (std::pow(QCD::CF,2)*(-27*(51 + 48*MATH::ZETA2 + 128*MATH::ZETA3) - 3294*HPL2(1,0,x) - 2592*MATH::ZETA2*HPL2(1,0,x) - 54*(27 + 16*MATH::ZETA2)*HPL2(1,1,x) + 324*HPL3(1,0,0,x) + 1296*HPL3(1,0,1,x) + 1944*HPL3(1,1,0,x) + 
        1944*HPL3(1,1,1,x) - 2592*HPL4(1,0,-1,0,x) + 864*HPL4(1,0,0,0,x) + 2592*HPL4(1,0,0,1,x) + 2592*HPL4(1,0,1,0,x) + 3024*HPL4(1,0,1,1,x) + 1296*HPL4(1,1,0,0,x) + 2592*HPL4(1,1,0,1,x) + 
        3456*HPL4(1,1,1,0,x) + 2592*HPL4(1,1,1,1,x)))/54.;
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
				+ 188.64 ) / (1.0 - z);
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
	// double HPL0z	=	HPL1(0,z);
	// double HPL1z	=	HPL1(1,z);
	// double HPLm1z	=	HPL1(-1,z);
	// double HPLm10z	=	HPL2(-1,0,z);
	// double HPL00z	=	HPL2(0,0,z);
	// double HPL01z	=	HPL2(0,1,z);
	// double HPL11z	=	HPL2(1,1,z);
	// double HPL10z	=	HPL2(1,0,z);
	// double HPLm1m10z	=	HPL3(-1,-1,0,z);
	// double HPLm100z	=	HPL3(-1,0,0,z);
	// double HPL0m10z	=	HPL3(0,-1,0,z);
	// double HPLm101z	=	HPL3(-1,0,1,z);
	// double HPL000z	=	HPL3(0,0,0,z);
	// double HPL001z	=	HPL3(0,0,1,z);
	// double HPL010z	=	HPL3(0,1,0,z);
	// double HPL011z	=	HPL3(0,1,1,z);
	// double HPL100z	=	HPL3(1,0,0,z);
	// double HPL101z	=	HPL3(1,0,1,z);
	// double HPL110z	=	HPL3(1,1,0,z);
	// double HPL111z	=	HPL3(1,1,1,z);
	double result(0.0);
	result += ((QCD::CF*QCD::NF*(247 - 144*MATH::ZETA2 + 342*HPL1(0,z) + 174*HPL1(1,z) + 180*HPL2(0,0,z) + 144*HPL2(0,1,z) + 72*HPL2(1,0,z) + 72*HPL2(1,1,z)))/27. + 
      (QCD::CA*QCD::CF*(-3155 + 1584*MATH::ZETA2 + 216*MATH::ZETA3 - 4302*HPL1(0,z) + 432*MATH::ZETA2*HPL1(0,z) - 2202*HPL1(1,z) + 1296*MATH::ZETA2*HPL1(1,z) - 1980*HPL2(0,0,z) - 1584*HPL2(0,1,z) - 
           792*HPL2(1,0,z) - 792*HPL2(1,1,z) - 1296*HPL3(0,-1,0,z) - 648*HPL3(0,0,0,z) - 432*HPL3(0,0,1,z) - 864*HPL3(1,0,0,z) - 432*HPL3(1,0,1,z) + 432*HPL3(1,1,0,z)))/
       54. + (std::pow(QCD::CF,2)*(51 + 48*MATH::ZETA2 + 128*MATH::ZETA3 + 122*HPL1(0,z) + 96*MATH::ZETA2*HPL1(0,z) + 54*HPL1(1,z) + 32*MATH::ZETA2*HPL1(1,z) - 12*HPL2(0,0,z) - 48*HPL2(0,1,z) - 
           72*HPL2(1,0,z) - 72*HPL2(1,1,z) + 96*HPL3(0,-1,0,z) - 32*HPL3(0,0,0,z) - 96*HPL3(0,0,1,z) - 96*HPL3(0,1,0,z) - 112*HPL3(0,1,1,z) - 48*HPL3(1,0,0,z) - 
           96*HPL3(1,0,1,z) - 128*HPL3(1,1,0,z) - 96*HPL3(1,1,1,z)))/2.)/(1 - z);
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
				+ 28.384 * L0
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
				+ 2.244 * std::pow(L0,2)
				+ 16./3. * L0
				+ 20./9 * std::pow(L0,2);
	return result;
}

/// @brief eq. (B.5) in [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005]
/// This expression was obtained after some rewriting and reordering with Mathematica.
double c2q_ns_2_01_reg_exact(double z)	{
	// double HPL0z	=	HPL1(0,z);
	// double HPL1z	=	HPL1(1,z);
	// double HPLm1z	=	HPL1(-1,z);
	// double HPLm10z	=	HPL2(-1,0,z);
	// double HPL00z	=	HPL2(0,0,z);
	// double HPL01z	=	HPL2(0,1,z);
	// double HPL11z	=	HPL2(1,1,z);
	// double HPL10z	=	HPL2(1,0,z);
	// double HPLm1m10z	=	HPL3(-1,-1,0,z);
	// double HPLm100z	=	HPL3(-1,0,0,z);
	// double HPL0m10z	=	HPL3(0,-1,0,z);
	// double HPLm101z	=	HPL3(-1,0,1,z);
	// double HPL000z	=	HPL3(0,0,0,z);
	// double HPL001z	=	HPL3(0,0,1,z);
	// double HPL010z	=	HPL3(0,1,0,z);
	// double HPL011z	=	HPL3(0,1,1,z);
	// double HPL100z	=	HPL3(1,0,0,z);
	// double HPL101z	=	HPL3(1,0,1,z);
	// double HPL110z	=	HPL3(1,1,0,z);
	// double HPL111z	=	HPL3(1,1,1,z);
	double result(0.0);
	result += (std::pow(QCD::CF,2)*(-93 - 209*z))/4. + (QCD::CA*QCD::CF*(139 + 3159*z))/36. - (4*std::pow(QCD::CF,2)*(7 + 13*z)*MATH::ZETA2)/5. - 4*std::pow(QCD::CF,2)*(1 - 19*z)*MATH::ZETA3 + 
   (8*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(-1 + 2*z)*(1 - HPL1(0,z)))/(5.*(-2 + z)) - 
   (std::pow(QCD::CF,2)*(43 + 63*z)*HPL1(0,z))/2. + (QCD::CA*QCD::CF*(71 + 323*z)*HPL1(0,z))/6. + (std::pow(QCD::CF,2)*(59 - 109*z)*HPL1(1,z))/2. - (17*QCD::CA*QCD::CF*(5 - 19*z)*HPL1(1,z))/6. + 
   (296*QCD::CF*(-0.5*QCD::CA + QCD::CF)*HPL2(-1,0,z))/5. + (8*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(1 - 1/z)*(-1 - 2*z)*HPL2(-1,0,z))/(5.*(-2 - z)) + (136*QCD::CF*(-0.5*QCD::CA + QCD::CF)*z*HPL2(-1,0,z))/5. + 
   (72*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(1 - z)*(1 + 2*z + 2*std::pow(z,2))*HPL2(-1,0,z))/5. - 
   (72*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(1 - 2*z + 2*std::pow(z,2))*(-1 - HPL1(0,z) + (1 + z)*(MATH::ZETA2 - HPL2(0,0,z))))/5. - (4*QCD::CA*QCD::CF*(9 + 16*z)*(MATH::ZETA2 - HPL2(0,0,z)))/5. + 
   (std::pow(QCD::CF,2)*(33 + 37*z)*HPL2(0,0,z))/5. + 2*std::pow(QCD::CF,2)*(5 + 9*z)*(2*HPL2(0,1,z) + HPL2(1,1,z)) + 
   QCD::NF*((QCD::CF*(-23 - 243*z))/18. - (QCD::CF*(7 + 19*z)*HPL1(0,z))/3. - (QCD::CF*(1 + 13*z)*HPL1(1,z))/3. + 
      (QCD::CF*(-247 + 144*MATH::ZETA2 - 342*HPL1(0,z) - 174*HPL1(1,z) - 180*HPL2(0,0,z) - 144*HPL2(0,1,z) - 72*HPL2(1,0,z) - 72*HPL2(1,1,z)))/54. - 
      (QCD::CF*z*(247 - 144*MATH::ZETA2 + 342*HPL1(0,z) + 174*HPL1(1,z) + 180*HPL2(0,0,z) + 144*HPL2(0,1,z) + 72*HPL2(1,0,z) + 72*HPL2(1,1,z)))/54.) - 
   8*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(1 + 5*z)*(MATH::ZETA2*HPL1(-1,z) + 2*HPL3(-1,-1,0,z) - HPL3(-1,0,0,z)) + 16*std::pow(QCD::CF,2)*HPL3(0,-1,0,z) - 8*QCD::CA*QCD::CF*(5*z*MATH::ZETA3 + HPL3(0,-1,0,z)) + 
   4*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(-1 + z + 2/(1 + z))*(7*MATH::ZETA3 - 8*MATH::ZETA2*HPL1(-1,z) - 2*HPL1(0,z) + 2*MATH::ZETA2*HPL1(0,z) - 8*HPL3(-1,-1,0,z) + 10*HPL3(-1,0,0,z) + 4*HPL3(-1,0,1,z) + 
      6*HPL3(0,-1,0,z) - 3*HPL3(0,0,0,z) - 2*HPL3(0,0,1,z)) + 2*std::pow(QCD::CF,2)*(1 + z)*
    (-4*MATH::ZETA2*HPL1(0,z) + 7*HPL2(1,0,z) + 5*HPL3(0,0,0,z) + 4*HPL3(0,0,1,z) + 2*HPL3(0,1,0,z) + 2*HPL3(0,1,1,z)) + 
   8*QCD::CF*(-0.5*QCD::CA + QCD::CF)*(1 - 5*z)*(-(MATH::ZETA2*HPL1(1,z)) + HPL3(1,0,0,z)) + 
   (QCD::CA*QCD::CF*(3155 - 1584*MATH::ZETA2 - 216*MATH::ZETA3 + 4302*HPL1(0,z) - 432*MATH::ZETA2*HPL1(0,z) + 2202*HPL1(1,z) - 1296*MATH::ZETA2*HPL1(1,z) + 1980*HPL2(0,0,z) + 1584*HPL2(0,1,z) + 
        792*HPL2(1,0,z) + 792*HPL2(1,1,z) + 1296*HPL3(0,-1,0,z) + 648*HPL3(0,0,0,z) + 432*HPL3(0,0,1,z) + 864*HPL3(1,0,0,z) + 432*HPL3(1,0,1,z) - 432*HPL3(1,1,0,z)))/108.\
    - (QCD::CA*QCD::CF*z*(-3155 + 1584*MATH::ZETA2 + 216*MATH::ZETA3 - 4302*HPL1(0,z) + 432*MATH::ZETA2*HPL1(0,z) - 2202*HPL1(1,z) + 1296*MATH::ZETA2*HPL1(1,z) - 1980*HPL2(0,0,z) - 1584*HPL2(0,1,z) - 
        792*HPL2(1,0,z) - 792*HPL2(1,1,z) - 1296*HPL3(0,-1,0,z) - 648*HPL3(0,0,0,z) - 432*HPL3(0,0,1,z) - 864*HPL3(1,0,0,z) - 432*HPL3(1,0,1,z) + 432*HPL3(1,1,0,z)))/108.\
    - (std::pow(QCD::CF,2)*z*(51 + 48*MATH::ZETA2 + 128*MATH::ZETA3 + 122*HPL1(0,z) + 96*MATH::ZETA2*HPL1(0,z) + 54*HPL1(1,z) + 32*MATH::ZETA2*HPL1(1,z) - 12*HPL2(0,0,z) - 48*HPL2(0,1,z) - 72*HPL2(1,0,z) - 
        72*HPL2(1,1,z) + 96*HPL3(0,-1,0,z) - 32*HPL3(0,0,0,z) - 96*HPL3(0,0,1,z) - 96*HPL3(0,1,0,z) - 112*HPL3(0,1,1,z) - 48*HPL3(1,0,0,z) - 96*HPL3(1,0,1,z) - 
        128*HPL3(1,1,0,z) - 96*HPL3(1,1,1,z)))/4. + (std::pow(QCD::CF,2)*(-51 - 48*MATH::ZETA2 - 128*MATH::ZETA3 - 122*HPL1(0,z) - 96*MATH::ZETA2*HPL1(0,z) - 54*HPL1(1,z) - 32*MATH::ZETA2*HPL1(1,z) + 
        12*HPL2(0,0,z) + 48*HPL2(0,1,z) + 72*HPL2(1,0,z) + 72*HPL2(1,1,z) - 96*HPL3(0,-1,0,z) + 32*HPL3(0,0,0,z) + 96*HPL3(0,0,1,z) + 96*HPL3(0,1,0,z) + 
        112*HPL3(0,1,1,z) + 48*HPL3(1,0,0,z) + 96*HPL3(1,0,1,z) + 128*HPL3(1,1,0,z) + 96*HPL3(1,1,1,z)))/4.;
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
	result += (QCD::CF*(-117 + 121*z))/3. - (4*QCD::CF*(111 - 176*z)*MATH::ZETA2)/15. - 4*QCD::CF*(1 - 18*z)*MATH::ZETA3 + (4*QCD::CF*(-1 + 2*z)*(1 - HPL1(0,z)))/(15.*(-2 + z)) + (QCD::CF*(16 - 61*z)*HPL1(0,z))/3. - 2*QCD::CF*(1 - 8*z)*HPL1(1,z) + 
   (4*QCD::CF*(217 + ((1 - 1/z)*(-1 - 2*z))/(-2 - z) + 117*z)*HPL2(-1,0,z))/15. - (48*QCD::CF*(1 - z)*(1 + 2*z + 2*std::pow(z,2))*HPL2(-1,0,z))/5. + (QCD::CF*(639 - 1004*z)*HPL2(0,0,z))/15. + 4*QCD::CF*(5 - 4*z)*HPL2(0,1,z) + 
   2*QCD::CF*(5 + 4*z)*(HPL2(1,0,z) + HPL2(1,1,z)) - 8*QCD::CF*(1 + 2*z)*(MATH::ZETA2*HPL1(-1,z) + 2*HPL3(-1,-1,0,z) - HPL3(-1,0,0,z)) + 
   8*QCD::CF*(1 + 2*z + 2*std::pow(z,2))*(-(MATH::ZETA2*HPL1(-1,z)) - 2*HPL3(-1,-1,0,z) + HPL3(-1,0,0,z) + 2*HPL3(0,-1,0,z)) + 
   2*QCD::CF*(1 - 2*z)*(-4*MATH::ZETA2*HPL1(0,z) - 4*MATH::ZETA2*HPL1(1,z) + 8*HPL3(0,-1,0,z) + 5*HPL3(0,0,0,z) + 4*HPL3(0,0,1,z) + 2*HPL3(0,1,0,z) + 2*HPL3(0,1,1,z) + 4*HPL3(1,0,0,z)) + 
   (2*QCD::CF*(1 - 2*z + 2*std::pow(z,2))*(-9 + 90*MATH::ZETA3 - 54*HPL1(0,z) + 60*MATH::ZETA2*HPL1(0,z) - 30*HPL1(1,z) + 40*MATH::ZETA2*HPL1(1,z) + 6*(19 + 4*z)*(MATH::ZETA2 - HPL2(0,0,z)) - 90*HPL2(0,1,z) - 90*HPL2(1,0,z) - 
        90*HPL2(1,1,z) - 50*HPL3(0,0,0,z) - 60*HPL3(0,0,1,z) - 40*HPL3(0,1,0,z) - 50*HPL3(0,1,1,z) - 30*HPL3(1,0,0,z) - 60*HPL3(1,0,1,z) - 40*HPL3(1,1,0,z) - 50*HPL3(1,1,1,z)))/5. + 
   QCD::CA*((7*(105 - 46*z))/6. - (2*(107 - 10*z)*MATH::ZETA2)/3. + ((1567 - 338*z)*HPL1(0,z))/9. + ((289 - 52*z)*HPL1(1,z))/3. + (8*((-1 - 2*z)/(-2 - z) - 3*(4 + 3*z))*HPL2(-1,0,z))/3. + 
      (4*(47 + 35*z)*HPL2(0,0,z))/3. + 2*(33 - 2*z)*HPL2(0,1,z) + 2*(23 - 6*z)*HPL2(1,0,z) + 6*(7 - 2*z)*HPL2(1,1,z) + 
      (4*(-1 + 2*z)*(43 - 18*MATH::ZETA2 - 39*HPL1(1,z) + 18*HPL2(1,0,z) + 18*HPL2(1,1,z)))/(27.*(-2 + z)) - 4*(1 + 2*z)*(-(MATH::ZETA2*HPL1(-1,z)) - 2*HPL3(-1,-1,0,z) + HPL3(-1,0,0,z) + 2*HPL3(0,-1,0,z)) + 
      (4*(1 + 2*z + 2*std::pow(z,2))*(-6*MATH::ZETA2*HPL1(-1,z) + 10*HPL2(-1,0,z) + 9*HPL3(-1,0,0,z) + 6*HPL3(-1,0,1,z) + 6*HPL3(0,-1,0,z)))/3. + 4*(5 + 14*z)*HPL3(0,0,0,z) - 
      8*(1 + 3*z)*(MATH::ZETA3 + 2*MATH::ZETA2*HPL1(0,z) - 2*HPL3(0,0,1,z)) + 8*(1 + 4*z)*(HPL3(0,1,0,z) + HPL3(0,1,1,z)) - 4*(1 - 2*z)*(-(MATH::ZETA2*HPL1(1,z)) + HPL3(1,0,0,z)) - 
      ((1 - 2*z + 2*std::pow(z,2))*(4493 - 3996*MATH::ZETA2 - 648*MATH::ZETA3 + 6270*HPL1(0,z) - 432*MATH::ZETA2*HPL1(0,z) + 4710*HPL1(1,z) - 216*MATH::ZETA2*HPL1(1,z) + 3492*HPL2(0,0,z) + 3996*HPL2(0,1,z) + 2412*HPL2(1,0,z) + 
           2196*HPL2(1,1,z) + 432*HPL3(0,0,1,z) + 432*HPL3(0,1,0,z) + 432*HPL3(0,1,1,z) + 432*HPL3(1,0,0,z) + 216*HPL3(1,0,1,z) + 648*HPL3(1,1,0,z) + 216*HPL3(1,1,1,z)))/54.);
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
	result += 12*QCD::CF*(1 - z)*(HPL2(1,0,z) + HPL2(1,1,z)) + QCD::CF*((22*(3 - 5*z))/3. - (8*(5 - z)*MATH::ZETA2)/3. + (8*(71 - 49*z)*HPL1(0,z))/9. + (4*(16 - 13*z)*HPL1(1,z))/3. + 
      (8*(-1 + (-1 - 2*z)/(-2 - z) - 2*z - 2*std::pow(z,2))*HPL2(-1,0,z))/3. + 8*(1 - 2*z)*HPL2(0,1,z) + 
      (8*(1 - 2*z + 2*std::pow(z,2))*(28 + 27*MATH::ZETA2 - 24*HPL1(0,z) + 6*HPL1(1,z) - 36*HPL2(0,0,z) - 27*HPL2(0,1,z) - 9*HPL2(1,0,z) - 9*HPL2(1,1,z)))/27. + 
      (4*(-1 + 2*z)*(43 - 18*MATH::ZETA2 - 39*HPL1(1,z) + 18*HPL2(1,0,z) + 18*HPL2(1,1,z)))/(27.*(-2 + z)) -
      (2*(1 + z)*(12*MATH::ZETA3 + 24*MATH::ZETA2*HPL1(0,z) + 12*HPL2(-1,0,z) - 13*HPL2(0,0,z) - 30*HPL3(0,0,0,z) - 24*HPL3(0,0,1,z) - 12*HPL3(0,1,0,z) - 12*HPL3(0,1,1,z)))/3.);
	return result;
}