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

double c2q_ns_1_0_reg(double z);
double c2q_ns_1_0_plus(double z);
double c2q_ns_1_0_local();
double c2q_ns_1_0_localplus(double x);
double c2g_1_0(double z);


///
///
/// O(alphas^2)

double c2q_ns_2_0_local_approx();
double c2q_ns_2_1_local_approx();
double c2q_ns_2_01_local_exact(double z);

double c2q_ns_2_0_localplus_approx(double z);
double c2q_ns_2_1_localplus_approx(double z);

double c2q_ns_2_0_plus_approx(double z);
double c2q_ns_2_1_plus_approx(double z);
double c2q_ns_2_01_plus_exact(double z);

double c2q_ns_2_0_reg_approx(double z);
double c2q_ns_2_1_reg_approx(double z);
double c2q_ns_2_01_reg_exact(double z);

double c2g_2_0_reg_approx(double z);
double c2g_2_0_local_approx();

double c2q_ps_2_0_reg_approx(double z);

#endif