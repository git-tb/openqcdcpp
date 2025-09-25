#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "protectedobject.h"

namespace QCDORDER	{
	ProtectedObject<int>	NFLOOPS{2};
	ProtectedObject<int>	F2ORDER{2};
	#define kordf2 QCDORDER::F2ORDER
};

namespace QCD	{
	const double 	CF	= 4.0/3.0,
					CG	= 3.0,
					TR	= 0.5;
	ProtectedObject<int>	NF{3};

	/// @brief Computes the sum $\sum_{q_i,\bar q_i}Q_i^2$ of squared fractional quark charges.
	double sumQquark()	{
		assert(3 <= QCD::NF);
		assert(QCD::NF <= 5);

		double result(0.0);

		if(QCD::NF == 3)
			result += 4./9. + 1./9. + 1./9.;
		else if(QCD::NF == 4)
			result += 4./9. + 1./9. + 1./9. + 4./9.;
		else if(QCD::NF == 5)
			result += 4./9. + 1./9. + 1./9. + 4./9. + 1./9.;
		
		return result;
	}
}

namespace PRECISION	{
	ProtectedObject<double>	EPSABS{1e-10},
							EPSREL{1e-5};
	ProtectedObject<int>	ITER{1e5};
}

#endif