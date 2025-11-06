#include "hqcoeff.h"
#include "constants.h"

#include <cmath>

/// @brief from [Riemersma, Smith, van Neerven; Phys. Let. B, 347 (1-2), 143, 1995] eq. (7)
double chLg_1_0(double eta, double chi)	{
	double result(0.0);
	double sqrteta		= std::sqrt(eta);
	double sqrt1peta	= std::sqrt(1+eta);
	result	+= M_PI/2.0 * QCD::TR * chi / std::pow(1+eta+chi/4,3) * (
		2*sqrt1peta*sqrteta
		- std::log((sqrt1peta+sqrteta)/(sqrt1peta-sqrteta))
	);
	return result;
}

/// @brief from [Riemersma, Smith, van Neerven; Phys. Let. B, 347 (1-2), 143, 1995] eq. (8)
double chTg_1_0(double eta, double chi) {
	double result(0.0);
	double sqrteta		= std::sqrt(eta);
	double sqrt1peta	= std::sqrt(1+eta);
	result	+= M_PI/2.0 * QCD::TR / std::pow(1+eta+chi/4,3) * (
		- 2.0* ( std::pow(1+eta-chi/4,2)+1+eta)*sqrteta/sqrt1peta
		+ (
			2.0 * std::pow(1+eta,2) + std::pow(chi,2)/8.0 + 1 + 2*eta
		)*(
			std::log((sqrt1peta+sqrteta)/(sqrt1peta-sqrteta))
		)
	);
	return result;
}

double ch2g_1_0(double eta, double chi) {
	return chLg_1_0(eta, chi) + chTg_1_0(eta,chi);
}