#ifndef HQCOEFF_H
#define HQCOEFF_H

#include <libInterpolate/Interpolate.hpp>

/// @brief returns index of gridpoint in xarr[Nx] to the left or equal to x
uint findIdx(const double& x, const double* xarr, const std::size_t& Nx);

/// @brief 1 dimensional quadratic interpolation (3pt Lagrange)
double myInterp1D_3pt(	const double& x, 
						const double x1, const double x2, const double x3,
						const double f1, const double f2, const double f3
					);
/**
 * @brief Find the interpolated value of a function f(x,y), where f is given
 * as a 2D grid farr[Nx][Ny] with respect to the 1D discretized coordinates xarr[Nx]
 * yarr[Ny];
 * 
 * NOTE: farr should be passed as &farr[0][0]
 */
double myInterp2D(const double& x, const double& y, const double* xarr, const std::size_t& Nx, const double* yarr, const std::size_t& Ny, const double* farr);

/**
 * @brief The naming of the heavy quark coefficient functions follows the scheme
 * 		chkf_ps/ns_j_l_desc
 * with	ch - "heavy quark coefficient function"
 * 		k - for structure function Fk
 * 		f - quark (q) or gluon (g)
 * 		ps/ns - pure singlet (ps) or non-singlet (ns) in case of f=q
 * 		j - associated QCD order alpha_s^j
 * 		l - associated power of ln(muR2/m2)
 * 		desc - further description (i.e. asymptotic, ...)
 */

///
///
/// O(alphas)
/// Following [Riemersma, Smith, van Neerven; Phys. Let. B, 347 (1-2), 143, 1995]

double chLg_1_0(double eta, double chi);
double chTg_1_0(double eta, double chi);
double ch2g_1_0(double eta, double chi);


///
///
/// O(alphas^2)
/// Following [Riemersma, Smith, van Neerven; Phys. Let. B, 347 (1-2), 143, 1995]
/// coefficient functions are given in terms of analytically known asymptotics ("asymp")
/// and a remainder which is small and stored in a tabulation


/// ================================
/// ====== GLUON CONTRIBUTION ======
/// ================================

///
/// O(ln(muR2/m2)^0)

double I_hqhelper(double chi);
double J_hqhelper(double chi);

/// Longitudinal part
double chL_g_2_0(double eta, double chi);
/// @todo Rename the G and E functions, it seems they dont depend on the label "g" and are used for the quark contribution as well
double chL_g_2_0_asympG(double eta, double chi);
double chL_g_2_0_A_asympE(double eta, double chi);
double chL_g_2_0_F_asympE(double eta, double chi);
double chL_g_2_0_A_interp(double eta, double chi);
double chL_g_2_0_F_interp(double eta, double chi);
////// tables and methods for interpolation
/// @todo maybe in all interpolations it would be smart to interpolate the log of the values, if they span several order of magnitude
const int Neta = 73;
const int Nchi = 49;
extern double chL_g_2_0_A_table[Nchi][Neta];
extern double chL_g_2_0_F_table[Nchi][Neta];
extern double ch_g_2_logetalist[Neta]; /// reuse those
extern double ch_g_2_logchilist[Nchi]; /// reuse those
extern bool chL_g_2_0_A_interper_initialized;
extern bool chL_g_2_0_F_interper_initialized;
extern _2D::BicubicInterpolator<double> chL_g_2_0_A_interper;
extern _2D::BicubicInterpolator<double> chL_g_2_0_F_interper;

/// Transversal part
double chT_g_2_0(double eta, double chi);
/// @todo Rename the G and E functions, it seems they dont depend on the label "g" and are used for the quark contribution as well
double chT_g_2_0_asympG(double eta, double chi);
double chT_g_2_0_A_asympE(double eta, double chi);
double chT_g_2_0_F_asympE(double eta, double chi);
double chT_g_2_0_A_interp(double eta, double chi);
double chT_g_2_0_F_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chT_g_2_0_A_table[49][73];
extern double chT_g_2_0_F_table[49][73];
extern bool chT_g_2_0_A_interper_initialized;
extern bool chT_g_2_0_F_interper_initialized;
extern _2D::BicubicInterpolator<double> chT_g_2_0_A_interper;
extern _2D::BicubicInterpolator<double> chT_g_2_0_F_interper;

/// 2 part
double ch2_g_2_0(double eta, double chi);


///
/// O(ln(muR2/m2)^1)

/// Longitudinal part
double chL_g_2_1(double eta, double chi);
/// @todo Rename the G and E functions, it seems they dont depend on the label "g" and are used for the quark contribution as well
double chL_g_2_1_asympG(double eta, double chi);
double chL_g_2_1_A_asympE(double eta, double chi);
double chL_g_2_1_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chL_g_2_1_table[49][73];
extern bool chL_g_2_1_interper_initialized;
extern _2D::BicubicInterpolator<double> chL_g_2_1_interper;

/// Transversal part
double chT_g_2_1(double eta, double chi);
/// @todo Rename the G and E functions, it seems they dont depend on the label "g" and are used for the quark contribution as well
double chT_g_2_1_asympG(double eta, double chi);
double chT_g_2_1_A_asympE(double eta, double chi);
double chT_g_2_1_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chT_g_2_1_table[49][73];
extern bool chT_g_2_1_interper_initialized;
extern _2D::BicubicInterpolator<double> chT_g_2_1_interper;

/// 2 part
double ch2_g_2_1(double eta, double chi);


/// ================================
/// ====== QUARK CONTRIBUTION ======
/// ================================

///
/// O(ln(muR2/m2)^0)

/// Hcoupl -- photon couples to heavy quark

/// Longitudinal part
double chL_q_2_0_Hcoupl(double eta, double chi);
double chL_q_2_0_Hcoupl_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chL_q_2_0_Hcoupl_table[49][73];
extern bool chL_q_2_0_Hcoupl_interper_initialized;
extern _2D::BicubicInterpolator<double> chL_q_2_0_Hcoupl_interper;

/// Transversal part
double chT_q_2_0_Hcoupl(double eta, double chi);
double chT_q_2_0_Hcoupl_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chT_q_2_0_Hcoupl_table[49][73];
extern bool chT_q_2_0_Hcoupl_interper_initialized;
extern _2D::BicubicInterpolator<double> chT_q_2_0_Hcoupl_interper;

/// 2 part
double ch2_q_2_0_Hcoupl(double eta, double chi);


/// Lcoupl -- photon couples to light quark

/// Longitudinal part
double chL_q_2_0_Lcoupl(double eta, double chi);
double chL_q_2_0_Lcoupl_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chL_q_2_0_Lcoupl_table[49][73];
extern bool chL_q_2_0_Lcoupl_interper_initialized;
extern _2D::BicubicInterpolator<double> chL_q_2_0_Lcoupl_interper;

/// Transversal part
double chT_q_2_0_Lcoupl(double eta, double chi);
double chT_q_2_0_Lcoupl_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chT_q_2_0_Lcoupl_table[49][73];
extern bool chT_q_2_0_Lcoupl_interper_initialized;
extern _2D::BicubicInterpolator<double> chT_q_2_0_Lcoupl_interper;

/// 2 part
double ch2_q_2_0_Lcoupl(double eta, double chi);

///
/// O(ln(muR2/m2)^1)

/// Hcoupl -- photon couples to heavy quark

/// Longitudinal part
double chL_q_2_1_Hcoupl(double eta, double chi);
double chL_q_2_1_Hcoupl_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chL_q_2_1_Hcoupl_table[49][73];
extern bool chL_q_2_1_Hcoupl_interper_initialized;
extern _2D::BicubicInterpolator<double> chL_q_2_1_Hcoupl_interper;

/// Transversal part
double chT_q_2_1_Hcoupl(double eta, double chi);
double chT_q_2_1_Hcoupl_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chT_q_2_1_Hcoupl_table[49][73];
extern bool chT_q_2_1_Hcoupl_interper_initialized;
extern _2D::BicubicInterpolator<double> chT_q_2_1_Hcoupl_interper;

/// 2 part
double ch2_q_2_1_Hcoupl(double eta, double chi);


/// Lcoupl -- photon couples to light quark

/// Longitudinal part
double chL_q_2_1_Lcoupl(double eta, double chi);
double chL_q_2_1_Lcoupl_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chL_q_2_1_Lcoupl_table[49][73];
extern bool chL_q_2_1_Lcoupl_interper_initialized;
extern _2D::BicubicInterpolator<double> chL_q_2_1_Lcoupl_interper;

/// Transversal part
double chT_q_2_1_Lcoupl(double eta, double chi);
double chT_q_2_1_Lcoupl_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chT_q_2_1_Lcoupl_table[49][73];
extern bool chT_q_2_1_Lcoupl_interper_initialized;
extern _2D::BicubicInterpolator<double> chT_q_2_1_Lcoupl_interper;

/// 2 part
double ch2_q_2_1_Lcoupl(double eta, double chi);

#endif