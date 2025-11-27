#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <assert.h>
#include "protectedobject.h"

namespace MATH	{
	extern const double ZETA2;
	extern const double ZETA3;
	extern const double ZETA4;
	extern const double LOG2;
}

namespace QCD	{
	extern const double 	CF;
	extern const double		CA;
	extern const double		TR;
	extern ProtectedObject<double[6]> QMASSES;
	extern ProtectedObject<double[6]> QCHARGES;

	/// @brief the sum $\sum_{q_i}Q_i^2$ of squared fractional quark charges over nf light flavors
	double sumQi2(int nf);

	/// @brief QCD beta function, O(alpha_s^2) taken from (Schwartz QFT, eq. 26.97)
	double beta0(int nf);

	/// @brief QCD beta function, O(alpha_s^3) taken from (Schwartz QFT, eq. 26.98)
	double beta1(int nf);
}

namespace PRECISION {

    /// Absolute integration error goal
    extern ProtectedObject<double> EPSABS;

    /// Relative integration error goal
    extern ProtectedObject<double> EPSREL;

    /// For x < XTHRESH, z-integrals over partonic structure functions
    /// are transformed logarithmically to sample more points at the
    /// boundaries at very large (z->1) and very small (z->x->0) values of z.
    /// For x > XTHRESH only the upper limit (z->1) is transformed and sampled
    /// in more detail.
    extern ProtectedObject<double> XTHRESH;

    /// z-integration is performed over (x,1-DELTA), as needed for the
    /// logarithmic variable transformation
    extern ProtectedObject<double> DELTA;

    /// Maximum number of subdivisions
    extern ProtectedObject<int> ITER;

}

namespace APPROX	{
	enum APPROXTYPE	{
		APPR1	=	1,
		APPR2	=	2,
		EXACT	=	3,
	};

	extern ProtectedObject<APPROX::APPROXTYPE> LEVEL;
}



#endif