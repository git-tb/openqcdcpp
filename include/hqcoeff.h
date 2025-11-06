#ifndef HQCOEFF_H
#define HQCOEFF_H

/**
 * @brief The naming of the heavy quark coefficient functions follows the scheme
 * 		chnf_ps/ns_j_k_desc
 * with	ch - "heavy quark coefficient function"
 * 		k - for structure function Fk
 * 		f - quark (q) or gluon (g)
 * 		ps/ns - pure singlet (ps) or non-singlet (ns) in case of f=q
 * 		j - associated QCD order (alpha_s/(4*Pi))^j
 * 		n - associated power of Nf
 * 		desc - further description (i.e. asymptotic, ...)
 */

///
///
/// O(alphas)

double chLg_1_0(double eta, double chi);
double chTg_1_0(double eta, double chi);
double ch2g_1_0(double eta, double chi);

#endif