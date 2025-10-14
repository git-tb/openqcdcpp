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
	result += - QCD::CF * ( 2.*MATH::ZETA2 + 9.);
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

/// @brief from [Zijlstra, van Neerven; Nucl. Phys. B. 383 (3), 525, 1992] eq. B.2
/// and [Zijlstra, van Neerven; Nucl. Phys. B. 272 (1), 127, 1991] eq. 9,10
// double c2q_ns_2_0_local()	{
// 	double result(0.0);
// 	result	+=	QCD::CF * QCD::CF * ( 
// 				6. * std::pow(MATH::ZETA2,2)
// 				- 78. * MATH::ZETA3
// 				+ 69. * MATH::ZETA2
// 				+ 331./8. );
// 	result	+=	QCD::CA * QCD::CF * (
// 				71./5. * std::pow(MATH::ZETA2,2)
// 				+ 140./3. * MATH::ZETA3
// 				- 251./3. * MATH::ZETA2
// 				- 5465./72. );	
// 	return result;
// }

/// @brief from [van Neerven, Vogt; Nucl. Phys. B. 568 (1-2), 263, 2000] eq. 3.2
double c2q_ns_2_1_local_approx()	{	
	return 46.8405;
}

/// @brief from [Zijlstra, van Neerven; Nucl. Phys. B. 383 (3), 525, 1992] eq. B.2
/// and [Zijlstra, van Neerven; Nucl. Phys. B. 272 (1), 127, 1991] eq. 9,10
// double c2q_ns_2_1_local()	{
// 	double result(0.0);
// 	result	+=	QCD::CF *	(
// 				4./3. * MATH::ZETA3
// 				+ 38./3. * MATH::ZETA2
// 				+ 457./36. );
// 	return result;
// }

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


/**
 * @brief The exact 2nd order non-singlet coefficient function from [Vogt, Moch; Nucl. Phys. B, 724 (1-2), 3, 2005]
 * pasted here for reference. This expression was obtained after some rewriting and reordering with Mathematica.
 */
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

double c2q_ns_2_01_local_exact(double z)	{
	double result(0.0);
	result += (std::pow(QCD::CF,2)*(331./8. + 69*MATH::ZETA2 + 6*std::pow(MATH::ZETA2,2) - 78*MATH::ZETA3) - QCD::CA*QCD::CF*(5465./72. + (251*MATH::ZETA2)/3. - (71*std::pow(MATH::ZETA2,2))/5. - (140*MATH::ZETA3)/3.) + 
      QCD::CF*QCD::NF*(457./36. + (38*MATH::ZETA2)/3. + (4*MATH::ZETA3)/3.));
	return result;
}

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