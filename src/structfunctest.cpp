#include <iostream>

#include "structfunc.h"

int main()	{
	std::cout << QCDORDER::NFLOOPS << std::endl;

	Pdf::initialize("ABMP16_3_nnlo", 0);
	Pdf::printLHAPDFinfo();

	return 0;
}