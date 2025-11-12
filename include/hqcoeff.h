#ifndef HQCOEFF_H
#define HQCOEFF_H

#include <libInterpolate/Interpolate.hpp>

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

/// following [Riemersma, Smith, van Neerven; Phys. Let. B, 347 (1-2), 143, 1995]
double chLg_1_0(double eta, double chi);
double chTg_1_0(double eta, double chi);
double ch2g_1_0(double eta, double chi);


///
///
/// O(alphas^2)
/// Following [Riemersma, Smith, van Neerven; Phys. Let. B, 347 (1-2), 143, 1995]
/// coefficient functions are given in terms of analytically known asymptotics ("asymp")
/// and a remainder which is small and stored in a tabulation


///
/// O(ln(muR2/m2)^0)

double I_hqhelper(double chi);
double J_hqhelper(double chi);

/// Longitudinal part
double chL_g_2_0(double eta, double chi);
double chL_g_2_0_asympG(double eta, double chi);
double chL_g_2_0_A_asympE(double eta, double chi);
double chL_g_2_0_F_asympE(double eta, double chi);
double chL_g_2_0_A_interp(double eta, double chi);
double chL_g_2_0_F_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chL_g_2_0_A_table[49][73];
extern double chL_g_2_0_F_table[49][73];
extern double ch_g_2_logetalist[73];
extern double ch_g_2_logchilist[49];
extern bool chL_g_2_0_A_interper_initialized;
extern bool chL_g_2_0_F_interper_initialized;
extern _2D::BicubicInterpolator<double> chL_g_2_0_A_interper;
extern _2D::BicubicInterpolator<double> chL_g_2_0_F_interper;

/// Transversal part
double chT_g_2_0(double eta, double chi);
double chT_g_2_0_asympG(double eta, double chi);
double chT_g_2_0_A_asympE(double eta, double chi);
double chT_g_2_0_F_asympE(double eta, double chi);
double chT_g_2_0_A_interp(double eta, double chi);
double chT_g_2_0_F_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chT_g_2_0_A_table[49][73];
extern double chT_g_2_0_F_table[49][73];
// extern double ch_g_2_logetalist[73]; /// reuse grid from above
// extern double ch_g_2_logchilist[49]; /// reuse grid from above
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
double chL_g_2_1_asympG(double eta, double chi);
double chL_g_2_1_A_asympE(double eta, double chi);
double chL_g_2_1_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chL_g_2_1_table[49][73];
/// extern double ch_g_2_logetalist[73]; /// reuse grid from above
/// extern double ch_g_2_logchilist[49]; /// reuse grid from above
extern bool chL_g_2_1_interper_initialized;
extern _2D::BicubicInterpolator<double> chL_g_2_1_interper;

/// Transversal part
double chT_g_2_1(double eta, double chi);
double chT_g_2_1_asympG(double eta, double chi);
double chT_g_2_1_A_asympE(double eta, double chi);
double chT_g_2_1_interp(double eta, double chi);
////// tables and methods for interpolation
extern double chT_g_2_1_table[49][73];
/// extern double ch_g_2_logetalist[73]; /// reuse grid from above
/// extern double ch_g_2_logchilist[49]; /// reuse grid from above
extern bool chT_g_2_1_interper_initialized;
extern _2D::BicubicInterpolator<double> chT_g_2_1_interper;

/// 2 part
double ch2_g_2_1(double eta, double chi);

#endif