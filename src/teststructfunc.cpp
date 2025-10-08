#include <iostream>
#include <cmath>	// use math functions from std::

#include "structfunc.h"

int main()	{
	Pdf::initialize("ABMP16_3_nnlo", 0);
	Pdf::printLHAPDFinfo();

	PRECISION::EPSABS.set(10e-6);
	PRECISION::EPSREL.set(10e-6);
	QCDORDER::F2ORDER.set(1);

	double Q2_arr[9] = {
		3e2, 5e2, 1e3, 1.5e3, 2e3, 3e3, 5e3, 8e3, 1.5e4
	};
	double XMIN_arr[9] = {
		8e-3, 1.3e-2, 1.3e-2, 3.2e-2, 3.2e-2, 8e-2, 8e-2, 0.13e0, 0.32e0
	};
	double XMAX_arr[9] = {
		8e-2, 0.13e0, 0.25e0, 0.25e0, 0.25e0, 0.4e0, 0.4e0, 0.4e0, 0.4e0
	};

	/// output formatting
	const int PREC = 6;
	const int WIDTH = 15;

	std::cout	<< std::scientific << std::setprecision(PREC);
	std::cout	<< std::setw(WIDTH) << "Q2"
				<< std::setw(WIDTH) << "x"
				<< std::setw(WIDTH) << "F2@nlo"
				<< std::setw(WIDTH) << "F2@nnlo"
				<< std::endl;

	for(int i = 0; i < 9; i++)	{
		double Q2 = Q2_arr[i];

		int NX = 10;
		double DLNX = (std::log(XMAX_arr[i])-std::log(XMIN_arr[i]))/(NX-1);
		for(int j = 0; j < NX; j++)	{
			double x = XMIN_arr[i] * std::exp(DLNX * j);
			
			QCDORDER::F2ORDER.set(1);
			double F2nlo	= F2(x,Q2);
			QCDORDER::F2ORDER.set(2);
			double F2nnlo	= F2(x,Q2);

			std::cout	<< std::setw(WIDTH) << Q2
						<< std::setw(WIDTH) << x
						<< std::setw(WIDTH) << F2nlo
						<< std::setw(WIDTH) << F2nnlo
						<< std::endl;
		}
		std::cout << std::endl;
	}
}