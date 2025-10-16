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
	forpreccontrol_.nf2qcd1	= 10;
	forpreccontrol_.nf2qcd2	= 10;
	
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
	const int PREC = 7;
	const int WIDTH = 23;

	std::cout	<< std::scientific << std::setprecision(PREC);
	std::cout	<< std::setw(WIDTH) << "Q2"
				<< std::setw(WIDTH) << "x"
				<< std::setw(WIDTH) << "F2@nlo(fortran/cpp)"
				<< std::setw(WIDTH) << "F2@nnlo(fortran/cpp)"
				<< std::endl;

	for(int i = 0; i < 9; i++)	{
		double Q2 = Q2_arr[i];

		int NX = 10;
		double DLNX = (std::log(XMAX_arr[i])-std::log(XMIN_arr[i]))/(NX-1);
		for(int j = 0; j < NX; j++)	{
			double x = XMIN_arr[i] * std::exp(DLNX * j);

			QCDORDER::F2ORDER.set(1);
			foralpsrenorm_.kordf2_ 	= 1;
			double F2_fortran_nlo	= f2qcd_(3,1,22,x,Q2);
			double F2_cpp_nlo		= F2(x,Q2);
			QCDORDER::F2ORDER.set(2);
			foralpsrenorm_.kordf2_ 	= 2;
			double F2_fortran_nnlo	= f2qcd_(3,1,22,x,Q2);
			APPROX::LEVEL.set(APPROX::APPR1);
			double F2_cpp_nnlo_apr1	= F2(x,Q2);
			APPROX::LEVEL.set(APPROX::APPR2);
			double F2_cpp_nnlo_apr2	= F2(x,Q2);
			APPROX::LEVEL.set(APPROX::EXACT);
			double F2_cpp_nnlo_ex	= F2(x,Q2);

			std::cout	<< std::setw(WIDTH) << Q2
						<< std::setw(WIDTH) << x
						<< std::setw(WIDTH) << F2_fortran_nlo
						<< std::setw(WIDTH) << F2_fortran_nnlo
						<< std::endl
						<< std::setw(WIDTH)	<< " "
						<< std::setw(WIDTH)	<< " "
						<< std::setw(WIDTH) << F2_cpp_nlo
						<< std::setw(WIDTH) << F2_cpp_nnlo_apr1
						<< std::endl
						<< std::setw(WIDTH)	<< " "
						<< std::setw(WIDTH)	<< " "
						<< std::setw(WIDTH) << " "
						<< std::setw(WIDTH) << F2_cpp_nnlo_apr2
						<< std::endl
						<< std::setw(WIDTH)	<< " "
						<< std::setw(WIDTH)	<< " "
						<< std::setw(WIDTH) << " "
						<< std::setw(WIDTH) << F2_cpp_nnlo_ex
						<< std::endl << std::endl;
		}
		std::cout << std::endl;
	}

	///
	/// Performance
	std::srand(std::time(NULL));
	const int NRUN 		= 1e3;
	double	time_for	= 0.0, ///< fortran timing
	time_cpp	= 0.0; ///< C++ timing
	auto start 			= std::chrono::high_resolution_clock::now();
	auto stop 			= std::chrono::high_resolution_clock::now();
	auto duration 		= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	
	double 	xmin	= 5e-6,
			xmax	= 0.999,
			q2min	= 1,
			q2max	= 1e6;

	///
	std::cout << "Performance test@NLO" << std::endl;
	QCDORDER::F2ORDER.set(1);
	foralpsrenorm_.kordf2_ 	= 1;
	for(int i = 0; i < NRUN; i++)	{
		{
			std::cout << "\33[2K\r" << std::flush;
			std::cout << (int)(100*(double)(i+1)/(double)NRUN) << "% done" << std::flush;
		}

		double	r1	= (double)std::rand()/(double)RAND_MAX,
				r2	= (double)std::rand()/(double)RAND_MAX;
		double	x	= xmin + r1*(xmax-xmin),
				q2	= q2min + r2*(q2max-q2min);

		start 		= std::chrono::high_resolution_clock::now();
		f2qcd_(3,1,22,x,q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_for 	+= (double)duration.count()/(double)NRUN;

		start 		= std::chrono::high_resolution_clock::now();
		F2(x,q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_cpp 	+= (double)duration.count()/(double)NRUN;
	}

	std::cout	<< std::endl
				<< std::setw(26) << "Fortran runtime avg (ns):"
				<< std::setw(8) << time_for
				<< std::endl
				<< std::setw(26) << "C++ runtime avg (ns):"
				<< std::setw(8) << time_cpp
				<< std::endl
				<< std::endl;

	
	///
	std::cout << "Performance test@NNLO" << std::endl;
	QCDORDER::F2ORDER.set(2);
	foralpsrenorm_.kordf2_ 	= 2;
	for(int i = 0; i < NRUN; i++)	{
		{
			std::cout << "\33[2K\r" << std::flush;
			std::cout << (int)(100*(double)(i+1)/(double)NRUN) << "% done" << std::flush;
		}

		double	r1	= (double)std::rand()/(double)RAND_MAX,
				r2	= (double)std::rand()/(double)RAND_MAX;
		double	x	= xmin + r1*(xmax-xmin),
				q2	= q2min + r2*(q2max-q2min);

		start 		= std::chrono::high_resolution_clock::now();
		f2qcd_(3,1,22,x,q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_for 	+= (double)duration.count()/(double)NRUN;

		start 		= std::chrono::high_resolution_clock::now();
		F2(x,q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_cpp 	+= (double)duration.count()/(double)NRUN;
	}

	std::cout	<< std::endl
				<< std::setw(26) << "Fortran runtime avg (ns):"
				<< std::setw(8) << time_for
				<< std::endl
				<< std::setw(26) << "C++ runtime avg (ns):"
				<< std::setw(8) << time_cpp
				<< std::endl
				<< std::endl;

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