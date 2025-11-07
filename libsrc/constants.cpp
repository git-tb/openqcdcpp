#include "constants.h"
#include "protectedobject.h"

#include <assert.h>
#include <cmath>

namespace QCDORDER	{
	ProtectedObject<int>	NFLOOPS{2};
	ProtectedObject<int>	F2ORDER{2};
	#define kordf2 QCDORDER::F2ORDER
};

namespace MATH	{
	const double ZETA2	=	M_PI*M_PI/6.;
	const double ZETA3	=	1.2020569031595942853997381615114499907649862923405;
	const double ZETA4	=	std::pow(M_PI,4)/90.;
	const double LOG2	=	0.69314718055994530941723212145817656807550013436026;
}

namespace QCD	{
	const double 	CF	= 4.0/3.0;
	const double	CA	= 3.0;
	const double	TR	= 0.5;

	/// @todo Get rid of NF as a global variable! it should be passed as a parameter to each relevant function call!
	ProtectedObject<int>	NF{3};
	ProtectedObject<double[6]>	QMASSES{{
		0,
		0,
		0,
		1.5,
		4.5,
		172
	}};

	ProtectedObject<double[6]>	QCHARGES{{
		-1./3.,
		 2./3.,
		-1./3.,
		 2./3.,
		-1./3.,
		 2./3.
	}};

	/// @brief computes the sum $\sum_{q_i}Q_i^2$ of squared fractional quark charges depending on the value of QCD::NF
	/// @todo formulate this in terms of QCHARGES
	double sumQi2()	{
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

	double beta0()	{
		assert(3 <= QCD::NF);
		assert(QCD::NF <= 5);
		
		return 11./3. * QCD::CA - 4./3. * QCD::TR * QCD::NF;
	}

	double beta1()	{
		assert(3 <= QCD::NF);
		assert(QCD::NF <= 5);
		
		return 34./3. * QCD::CA * QCD::CA - 20./3. * QCD::CA * QCD::TR * QCD::NF - 4. * QCD::CF * QCD::TR * QCD::NF;
	}
}

namespace PRECISION {

    /// Absolute integration error goal
    ProtectedObject<double> EPSABS{1e-10};

    /// Relative integration error goal
    ProtectedObject<double> EPSREL{1e-10};

    /// For x < XTHRESH, z-integrals over partonic structure functions
    /// are transformed logarithmically to sample more points at the
    /// boundaries at very large (z->1) and very small (z->x->0) values of z.
    /// For x > XTHRESH only the upper limit (z->1) is transformed and sampled
    /// in more detail.
    ProtectedObject<double> XTHRESH{0.1};

    /// z-integration is performed over (x,1-DELTA), as needed for the
    /// logarithmic variable transformation
    ProtectedObject<double> DELTA{1e-8};

    /// Maximum number of subdivisions
    ProtectedObject<int> ITER{1000000};

}

namespace APPROX	{
	ProtectedObject<APPROX::APPROXTYPE> LEVEL{APPROX::APPR1};
}
