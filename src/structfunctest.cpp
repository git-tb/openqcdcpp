#include <iostream>
#include <cmath>	// use math functions from std::

#include "structfunc.h"

int main()	{
	Pdf::initialize("ABMP16_3_nnlo", 0);
	Pdf::printLHAPDFinfo();

	// int NX = 15;
	// double LOGXMIN(-3), LOGXMAX(-0.01);
	// double DLOGX = (LOGXMAX-LOGXMIN)/(NX-1);

	// int NQ2 = 25;
	// double LOGQ2MIN(0.5), LOGQ2MAX(5);
	// double DLOGQ2 = (LOGQ2MAX-LOGQ2MIN)/(NQ2-1);

	// std::cout << "Q2\t\tx\t\tF2(x,Q2)\tF2-log10x" << std::endl;
	// std::cout << std::scientific << std::setprecision(3);
	// for(double logQ2 = LOGQ2MIN; logQ2 <= LOGQ2MAX; logQ2 += DLOGQ2)	{
	// 	double Q2 = std::pow(10,logQ2);
	// 	for(double logx = LOGXMIN; logx <= LOGXMAX; logx += DLOGX)	{
	// 		double x = std::pow(10,logx);
	// 		std::cout << Q2 << "\t" << x << "\t" << F2(x,Q2) << "\t" << F2(x,Q2) - logx << std::endl;
	// 	}
	// 	std::cout << std::endl;
	// }

	// double Q2 = 300;
	// double x = 0.008;
	// std::cout << Q2 << "\t" << x << "\t" << F2(x,Q2) << std::endl;

	double Q2_arr[9] = {
		3e2, 5e2, 1e3, 1.5e3, 2e3, 3e3, 5e3, 8e3, 1.5e4
	};
	double XMIN_arr[9] = {
		8e-3, 1.3e-2, 1.3e-2, 3.2e-2, 3.2e-2, 8e-2, 8e-2, 0.13e0, 0.32e0
	};
	double XMAX_arr[9] = {
		8e-2, 0.13e0, 0.25e0, 0.25e0, 0.25e0, 0.4e0, 0.4e0, 0.4e0, 0.4e0
	};

	std::cout << "Q2\t\tx\t\tF2(x,Q2)" << std::endl;
	std::cout << std::scientific << std::setprecision(3);

	for(int i = 0; i < 9; i++)	{
		double Q2 = Q2_arr[i];

		int NX = 10;
		double DLNX = (std::log(XMAX_arr[i])-std::log(XMIN_arr[i]))/(NX-1);
		for(int j = 0; j < NX; j++)	{
			double x = XMIN_arr[i] * std::exp(DLNX * j);
			std::cout << Q2 << "\t" << x << "\t" << F2(x,Q2) << std::endl;
		}
		std::cout << std::endl;
	}
}