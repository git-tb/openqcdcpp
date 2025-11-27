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

/// @brief quick redifiniton such that we can pass literal constants (e.g. 1.234) as parameters
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
	forpreccontrol_.nf2qcd1	= 5;
	forpreccontrol_.nf2qcd2	= 5;
	int NF					= 3;
	int ORDERALPS			= 1;

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
			xmax	= 0.99;

	int	NQ2	=	10,
		Nx	=	10;

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

	for(int i = 0; i < Nx * NQ2; i++)	{
		int ix	= i%Nx;
		int iQ2	= i/Nx;
		{
			std::cout << "\33[2K\r" << std::flush;
			std::cout << "[\33[32m";
			for(int j = 0; j < 100; j++)	{
				if((int)(100*(double)(i+1)/(double)(Nx*NQ2)) >= j)	{
					std::cout << "\u2589";
				} else	{
					std::cout << "\u2591";
				}
			}
			std::cout << "\33[0m]" << std::flush;
		}

		double logQ2	= logQ2min + (double)iQ2/(double)(NQ2-1) * (logQ2max - logQ2min);
		double Q2		= std::exp(logQ2);
		double logx	= logxmin + (double)ix/(double)(Nx-1) * (logxmax - logxmin);
		double x	= std::exp(logx);

		ORDERALPS				= 2;
		foralpsrenorm_.kordf2 	= 1;
		double F2_fortran_nlo	= f2qcd_(3,1,22,x,Q2);
		double F2_cpp_nlo		= F2(x,Q2,ORDERALPS, NF);
		ORDERALPS				= 2;
		foralpsrenorm_.kordf2 	= 2;
		double F2_fortran_nnlo	= f2qcd_(3,1,22,x,Q2);
		APPROX::LEVEL.set(APPROX::APPR1);
		double F2_cpp_nnlo_apr1	= F2(x,Q2,ORDERALPS, NF);
		APPROX::LEVEL.set(APPROX::APPR2);
		double F2_cpp_nnlo_apr2	= F2(x,Q2,ORDERALPS, NF);
		APPROX::LEVEL.set(APPROX::EXACT);
		double F2_cpp_nnlo_ex	= F2(x,Q2,ORDERALPS, NF);

		fileout		<< Q2				<< ";"
					<< x				<< ";"
					<< F2_fortran_nlo	<< ";"
					<< F2_cpp_nlo		<< ";"
					<< F2_fortran_nnlo	<< ";"
					<< F2_cpp_nnlo_apr1	<< ";"
					<< F2_cpp_nnlo_apr2	<< ";"
					<< F2_cpp_nnlo_ex	<< std::endl;
	}
	std::cout << std::endl;
	fileout.close();

	///
	/// Performance
	std::srand(std::time(NULL));
	const int NRUN 		= 1e4;
	double	time_for	= 0.0, ///< fortran timing
			time_cppa1	= 0.0, ///< C++ timing, approximation 1
			time_cppa2	= 0.0, ///< C++ timing, approximation 2
			time_cppex	= 0.0; ///< C++ timing, exact
	auto start 			= std::chrono::high_resolution_clock::now();
	auto stop 			= std::chrono::high_resolution_clock::now();
	auto duration 		= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	///
	std::cout << "Performance test@NLO" << std::endl;
	ORDERALPS				= 1;
	foralpsrenorm_.kordf2 	= 1;
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
		F2(x,Q2, ORDERALPS, NF);
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
	ORDERALPS				= 2;
	foralpsrenorm_.kordf2 	= 2;
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
		F2(x,Q2, ORDERALPS, NF);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_cppa1 	+= (double)duration.count()/(double)NRUN;

		APPROX::LEVEL.set(APPROX::APPR2);
		start 		= std::chrono::high_resolution_clock::now();
		F2(x,Q2, ORDERALPS, NF);
		stop 		= std::chrono::high_resolution_clock::now();
		duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		time_cppa2 	+= (double)duration.count()/(double)NRUN;

		APPROX::LEVEL.set(APPROX::EXACT);
		start 		= std::chrono::high_resolution_clock::now();
		F2(x,Q2, ORDERALPS, NF);
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