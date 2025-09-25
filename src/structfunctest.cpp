#include <iostream>

#include "structfunc.h"

int main()	{
	std::cout << QCDORDER::NFLOOPS << std::endl;

	Pdf::initialize("ABMP16_3_nnlo", 0);
	Pdf::printLHAPDFinfo();

	int NX = 20;
	double XMIN(1e-3), XMAX(1-1e-3);
	for(int ix = 0; ix < NX; ix++)	{
		double x = XMIN + (double)ix/(double)(NX-1) * (XMAX - XMIN);
		std::cout << x << "\t" << F2(x,1) << std::endl;
	}

	return 0;
}