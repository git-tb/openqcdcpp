#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>

#include "LHAPDF/LHAPDF.h"
#include "structfunc.h"
#include "fortransymbols.h"

/// @brief quick redifiniton such that we can pass literal constants (e.g. 1.234) ass parameters
double f2qcd_(int nb, int nt, int ni, double xb, double q2) 	{
	return f2qcd_(&nb, &nt, &ni, &xb, &q2); ///< as defined in fortransymbols.h
}

int main()	{
	Pdf::initialize("ABMP16_3_nnlo", 0);
	Pdf::setSampling(SAMPLINGMETHOD::fromLHAPDF);
	Pdf::printLHAPDFinfo();

	PRECISION::EPSABS.set(1e-5);
	PRECISION::EPSREL.set(1e-5);
	PRECISION::ITER.set(1000);
	forpreccontrol_.nf2qcd1	= 3;
	forpreccontrol_.nf2qcd2	= 3;
	
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

	double	Q2min	= 2,
			Q2max	= 1e6,
			xmin	= 1e-5,
			xmax	= 0.999;

	int	NQ2full	=	50,	///< for file output
		NQ2red	=	10,	///< fir console output
		Nxfull	=	50,	///< for file output
		Nxred	=	10;	///< fir console output

	double	logQ2min	= std::log(Q2min),
			logQ2max	= std::log(Q2max),
			logxmin		= std::log(xmin),
			logxmax		= std::log(xmax);


	/// output formatting
	const int PREC = 7;
	const int WIDTH = 23;

	std::ofstream fileout("output.dat");
	fileout << std::scientific << std::setprecision(PREC);
	fileout << "Q2;x;F2@nlo(fortran);F2@nlo(cpp);F2@nnlo(fortran);F2@nnlo(cpp,apr1);F2@nnlo(cpp,apr2);F2@nnlo(cpp,ex)"
			<< std::endl;

	std::cout	<< std::scientific << std::setprecision(PREC);
	std::cout	<< std::setw(WIDTH) << "Q2"
				<< std::setw(WIDTH) << "x"
				<< std::setw(WIDTH) << "F2@nlo"
				<< std::setw(WIDTH) << "F2@nnlo"
				<< std::endl
				<< std::setw(WIDTH) << " "
				<< std::setw(WIDTH) << " "
				<< std::setw(WIDTH) << "fortran"
				<< std::setw(WIDTH) << "fortran"
				<< std::endl
				<< std::setw(WIDTH) << " "
				<< std::setw(WIDTH) << " "
				<< std::setw(WIDTH) << "cpp"
				<< std::setw(WIDTH) << "cpp (apr1)"
				<< std::endl
				<< std::setw(WIDTH) << " "
				<< std::setw(WIDTH) << " "
				<< std::setw(WIDTH) << " "
				<< std::setw(WIDTH) << "cpp (apr2)"
				<< std::endl
				<< std::setw(WIDTH) << " "
				<< std::setw(WIDTH) << " "
				<< std::setw(WIDTH) << " "
				<< std::setw(WIDTH) << "cpp (exct)"
				<< std::endl;

	for(int iQ2 = 0; iQ2 < NQ2full; iQ2++)	{
		double logQ2	= logQ2min + (double)iQ2/(double)NQ2full * (logQ2max - logQ2min);
		double Q2		= std::exp(logQ2);
		for(int ix = 0; ix < Nxfull; ix++)	{
			double logx	= logxmin + (double)ix/(double)Nxfull * (logxmax - logxmin);
			double x	= std::exp(logx);

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

			fileout		<< Q2				<< ";"
						<< x				<< ";"
						<< F2_fortran_nlo	<< ";"
						<< F2_cpp_nlo		<< ";"
						<< F2_fortran_nnlo	<< ";"
						<< F2_cpp_nnlo_apr1	<< ";"
						<< F2_cpp_nnlo_apr2	<< ";"
						<< F2_cpp_nnlo_ex	<< std::endl;

			if(	iQ2%(int)std::ceil((double)NQ2full/(double)NQ2red) == 0
				&& ix%(int)std::ceil((double)Nxfull/(double)Nxred) == 0)	{
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
		}
		if(	iQ2%(int)std::ceil((double)NQ2full/(double)NQ2red) == 0 ) std::cout << std::endl;
	}
	fileout.close();

	///
	/// Performance
	std::srand(std::time(NULL));
	const int NRUN 		= 1e3;
	double	time_for	= 0.0, ///< fortran timing
			time_cppa1	= 0.0, ///< C++ timing, approximation 1
			time_cppa2	= 0.0, ///< C++ timing, approximation 2
			time_cppex	= 0.0; ///< C++ timing, exact
	auto start 			= std::chrono::high_resolution_clock::now();
	auto stop 			= std::chrono::high_resolution_clock::now();
	auto duration 		= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

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
				Q2	= Q2min + r2*(Q2max-Q2min);

		start 		= std::chrono::high_resolution_clock::now();
		f2qcd_(3,1,22,x,Q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_for 	+= (double)duration.count()/(double)NRUN;

		APPROX::LEVEL.set(APPROX::EXACT);
		start 		= std::chrono::high_resolution_clock::now();
		F2(x,Q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_cppex 	+= (double)duration.count()/(double)NRUN;
	}

	std::cout	<< std::endl
				<< std::setw(26) << "Fortran runtime avg (ns):"
				<< std::setw(15) << time_for
				<< std::endl
				<< std::setw(26) << "C++ runtime avg (ns):"
				<< std::setw(15) << time_cppex
				<< std::endl
				<< std::endl;

	
	///
	std::cout << "Performance test@NNLO" << std::endl;
	APPROX::LEVEL.set(APPROX::APPR1);
	QCDORDER::F2ORDER.set(2);
	foralpsrenorm_.kordf2_ 	= 2;
	time_for	= 0.0, ///< fortran timing
	time_cppa1	= 0.0, ///< C++ timing, approximation 1
	time_cppa2	= 0.0, ///< C++ timing, approximation 2
	time_cppex	= 0.0; ///< C++ timing, exact
	for(int i = 0; i < NRUN; i++)	{
		{
			std::cout << "\33[2K\r" << std::flush;
			std::cout << (int)(100*(double)(i+1)/(double)NRUN) << "% done" << std::flush;
		}

		double	r1	= (double)std::rand()/(double)RAND_MAX,
				r2	= (double)std::rand()/(double)RAND_MAX;
		double	x	= xmin + r1*(xmax-xmin),
				Q2	= Q2min + r2*(Q2max-Q2min);

		start 		= std::chrono::high_resolution_clock::now();
		f2qcd_(3,1,22,x,Q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_for 	+= (double)duration.count()/(double)NRUN;

		APPROX::LEVEL.set(APPROX::APPR1);
		start 		= std::chrono::high_resolution_clock::now();
		F2(x,Q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_cppa1 	+= (double)duration.count()/(double)NRUN;

		APPROX::LEVEL.set(APPROX::APPR2);
		start 		= std::chrono::high_resolution_clock::now();
		F2(x,Q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_cppa2 	+= (double)duration.count()/(double)NRUN;

		APPROX::LEVEL.set(APPROX::EXACT);
		start 		= std::chrono::high_resolution_clock::now();
		F2(x,Q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_cppex 	+= (double)duration.count()/(double)NRUN;
	}

	std::cout	<< std::endl
				<< std::setw(26) << "Fortran runtime avg (ns):"
				<< std::setw(15) << time_for
				<< std::endl
				<< std::setw(26) << "C++ runtime avg (ns):"
				<< std::setw(15) << time_cppa1
				<< std::setw(15) << "(appr1)"
				<< std::endl
				<< std::setw(26) << " "
				<< std::setw(15) << time_cppa2
				<< std::setw(15) << "(appr2)"
				<< std::endl
				<< std::setw(26) << " "
				<< std::setw(15) << time_cppex
				<< std::setw(15) << "(exct)"
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