#include <iostream>
#include <cmath>	// use math functions from std::

#include "structfunc.h"
#include "constants.h"

int main()	{
	Pdf::initialize("ABMP16_3_nnlo", 0);
	Pdf::printLHAPDFinfo();

	int NX = 5;
	double LOGXMIN(-5), LOGXMAX(-0.01);
	double DLOGX = (LOGXMAX-LOGXMIN)/(NX-1);

	int NQ2 = 5;
	double LOGQ2MIN(0.5), LOGQ2MAX(5);
	double DLOGQ2 = (LOGQ2MAX-LOGQ2MIN)/(NQ2-1);

	int NZ = 5;

	int WIDTH(20);
	int PREC(10);

	QCDORDER::F2ORDER.set(1);

	std::cout << std::scientific << std::setprecision(PREC);
	for(int iq2 = 0; iq2 < NQ2; iq2++)	{
		double logQ2	= LOGQ2MIN + iq2 * DLOGQ2;
		double Q2 		= std::pow(10,logQ2);
		std::cout	<< "Q2=" << Q2 << std::endl;
		std::cout	<< std::setw(WIDTH) << "x" << std::setw(WIDTH) << "z" << std::setw(WIDTH) << "F2integrand" << std::endl;
		for(int ix = 0; ix < NX; ix++)	{
			double logx	= LOGXMIN + ix * DLOGX;
			double x 	= std::pow(10,logx);

			double LOGZMIN(logx+1e-8), LOGZMAX(LOGXMAX+1e-7);
			double DLOGZ = (LOGZMAX-LOGZMIN)/(NZ-1);
			for(int iz = 0; iz < NZ; iz++)	{
				double logz	= LOGZMIN + iz * DLOGZ;
				double z = std::pow(10, logz);

				std::cout	<< std::setw(WIDTH) << x
							<< std::setw(WIDTH) << z
							<< std::setw(WIDTH) << F2integrand(z, Q2, x)
							<< std::setw(WIDTH) << xfiQi2sum(x,Q2)
							<< std::endl;
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
}