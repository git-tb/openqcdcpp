#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <ctime>

#include "LHAPDF/LHAPDF.h"
#include "structfunc.h"
#include "fortransymbols.h"

/// @brief quick redifiniton such that we can pass literal constants (e.g. 1.234) ass parameters
double f2qcd_(int nb, int nt, int ni, double xb, double q2) 	{
	return f2qcd_(&nb, &nt, &ni, &xb, &q2);
}

int main()	{
	Pdf::initialize("ABMP16_3_nnlo", 0);
	// Pdf::setSampling(SAMPLINGMETHOD::fromLHAPDF);
	Pdf::setSampling(SAMPLINGMETHOD::fromOPENQCDRAD);
	Pdf::printLHAPDFinfo();

	PRECISION::EPSABS.set(1e-5);
	PRECISION::EPSREL.set(1e-5);
	PRECISION::ITER.set(1000);
	forpreccontrol_.nf2qcd1	= 4;
	forpreccontrol_.nf2qcd2	= 4;
	
	QCDORDER::F2ORDER.set(1);
	foralpsrenorm_.kordf2_ 	= 1;
	qcdpar_.nfc 			= 3;
	qcdpar_.cf				= 4./3.;
	qcdpar_.qsum[0]			= 1./9.;
	qcdpar_.qsum[1]			= 5./9.;
	qcdpar_.qsum[2]			= 6./9.;
	qcdpar_.qsum[3]			= 10./9.;
	qcdpar_.qsum[4]			= 12./9.;
	qcdpar_.qsum[5]			= 16./9.;

	int NX = 5;
	double LOGXMIN(-5), LOGXMAX(-0.01);
	double DLOGX = (LOGXMAX-LOGXMIN)/(NX-1);

	int NQ2 = 5;
	double LOGQ2MIN(0.5), LOGQ2MAX(5);
	double DLOGQ2 = (LOGQ2MAX-LOGQ2MIN)/(NQ2-1);

	int NZ = 5;

	/// output formatting
	const int PREC = 7;
	const int WIDTH = 23;

	for(int iq2 = 0; iq2 < NQ2; iq2++)	{
		double logQ2	= LOGQ2MIN + iq2 * DLOGQ2;
		double Q2 		= std::pow(10,logQ2);
		std::cout	<< "Q2=" << Q2 << std::endl;
		std::cout	<< std::setw(WIDTH) << "x" 
					<< std::setw(WIDTH) << "z"
					<< std::setw(WIDTH) << "F2i@nlo"
					<< std::setw(WIDTH) << "F2i@nnlo"
					<< std::endl;
		for(int ix = 0; ix < NX; ix++)	{
			double logx	= LOGXMIN + ix * DLOGX;
			double x 	= std::pow(10,logx);

			double LOGZMIN(logx+1e-8), LOGZMAX(LOGXMAX+1e-7);
			double DLOGZ = (LOGZMAX-LOGZMIN)/(NZ-1);
			for(int iz = 0; iz < NZ; iz++)	{
				double logz	= LOGZMIN + iz * DLOGZ;
				double z = std::pow(10, logz);

				QCDORDER::F2ORDER.set(1);
				double F2i_nlo			= F2integrand(z, Q2, x);
				QCDORDER::F2ORDER.set(2);
				APPROX::LEVEL.set(APPROX::APPR1);
				double F2i_nnlo_apr1	= F2integrand(z, Q2, x);
				APPROX::LEVEL.set(APPROX::APPR2);
				double F2i_nnlo_apr2	= F2integrand(z, Q2, x);
				APPROX::LEVEL.set(APPROX::EXACT);
				double F2i_nnlo_ex		= F2integrand(z, Q2, x);

				std::cout	<< std::setw(WIDTH) << x
							<< std::setw(WIDTH) << z
							<< std::setw(WIDTH) << F2i_nlo
							<< std::setw(WIDTH) << F2i_nnlo_apr1
							<< std::endl
							<< std::setw(WIDTH) << " "
							<< std::setw(WIDTH) << " "
							<< std::setw(WIDTH) << " "
							<< std::setw(WIDTH) << F2i_nnlo_apr2
							<< std::endl
							<< std::setw(WIDTH) << " "
							<< std::setw(WIDTH) << " "
							<< std::setw(WIDTH) << " "
							<< std::setw(WIDTH) << F2i_nnlo_ex
							<< std::endl;
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
	
	Pdf::destroy();
}


// /**
//  * @brief test performance of openQCDrad F2 vs my own F2 @ NNLO
//  */
// int main(int argc, char** argv)	{
// 	/// output formatting
// 	const int PREC = 15;
// 	const int WIDTH = 22;

// 	///
// 	std::string pdfset = "ABMP16_3_nnlo";
// 	Pdf::initialize(pdfset, 0);

// 	/// initialization of openQCDrad grid variables
// 	initgridconst_();
// 	mypdffillgrid_witharg_(pdfset.c_str(), pdfset.length());

// 	const double 	xmin	= 5e-6,
// 					xmax	= 0.9999;
// 	const double 	q2min	= 1.0,
// 					q2max	= 1e6;

// 	/// Performance
// 	std::srand(std::time(NULL));
// 	const int NRUN 		= 1e4;
// 	double	time_for	= 0.0, ///< fortran timing
// 			time_cpp	= 0.0; ///< C++ timing
// 	auto start 			= std::chrono::high_resolution_clock::now();
// 	auto stop 			= std::chrono::high_resolution_clock::now();
// 	auto duration 		= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

// 	QCDORDER::F2ORDER.set(2);
// 	foralpsrenorm_.kordf2_ = 2;

// 	for(int i = 0; i < NRUN; i++)	{
// 		{
// 			std::cout << "\33[2K\r" << std::flush;
// 			std::cout << (int)(100*(double)(i+1)/(double)NRUN) << "% done" << std::flush;
// 		}
// 		double r1	= (double)std::rand()/(double)RAND_MAX;
// 		double r2	= (double)std::rand()/(double)RAND_MAX;

// 		double x	= xmin + (xmax-xmin)*r1;
// 		double q2	= q2min + (q2max-q2min)*r2;

// 		start 			= std::chrono::high_resolution_clock::now();
// 		F2(x, q2);
// 		stop 			= std::chrono::high_resolution_clock::now();
// 		duration		= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
// 		time_cpp		+= (double)duration.count();

// 		start 		= std::chrono::high_resolution_clock::now();
// 		f2qcd_(3,1,22,x,q2);
// 		stop 		= std::chrono::high_resolution_clock::now();
// 		duration	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
// 		time_for	+= (double)duration.count();
// 	}

// 	std::cout	<< std::endl
// 				<< std::setw(WIDTH) << "cpp:"	
// 				<< std::setw(WIDTH) << time_cpp/(double)NRUN << " ns"
// 				<< std::endl
// 				<< std::setw(WIDTH) << "fortran:"	
// 				<< std::setw(WIDTH) << time_for/(double)NRUN << " ns"
// 				<< std::endl;

// 	return 0;
// }