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
double f2charm_ffn(double xb, double q2, int nq)	{
	return f2charm_ffn_(&xb, &q2, &nq); ///< as defined in fortransymbols.h
}

int main()	{
	Pdf::initialize("ABMP16_3_nnlo", 0);
	Pdf::setSampling(SAMPLINGMETHOD::fromOPENQCDRAD);
	Pdf::printLHAPDFinfo();

	PRECISION::EPSABS.set(1e-5);
	PRECISION::EPSREL.set(1e-5);
	PRECISION::ITER.set(1000);
	forpreccontrol_.nf2qcd1	= 5;
	forpreccontrol_.nf2qcd2	= 5;
	
	QCDORDER::F2ORDER.set(1);
	foralpsrenorm_.kordhq 	= 1;
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
	fileout << "Q2"	<< ";"
			<< "x"	<< ";"
			<< "F2heavy@NLO(Fortran)"	<< ";"
			<< "F2heavy@NLO(C++)"		<< ";"
			<< "F2heavy@NNLO(Fortran)"	<< ";"
			<< "F2heavy@NNLO(C++)"		<< std::endl;

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
		double logx		= logxmin + (double)ix/(double)(Nx-1) * (logxmax - logxmin);
		double x		= std::exp(logx);

		QCDORDER::F2ORDER.set(1);
		foralpsrenorm_.kordhq 	= 1;
		double F2_fortran_nlo	= f2charm_ffn(x,Q2,8/*=charm*/);	///< in openQCDrad parton flavors are
																	///		{1,2,3,4,5,6,7,8,...} -> 
																	///		{g,d,db,u,ub,s,sb,c,...}
		double F2_cpp_nlo		= F2heavy(x,Q2,3);
		QCDORDER::F2ORDER.set(2);
		foralpsrenorm_.kordhq	= 2;
		double F2_fortran_nnlo	= f2charm_ffn(x,Q2,8/*=charm*/);
		double F2_cpp_nnlo		= F2heavy(x,Q2,3);

		fileout		<< Q2				<< ";"
					<< x				<< ";"
					<< F2_fortran_nlo	<< ";"
					<< F2_cpp_nlo		<< ";"
					<< F2_fortran_nnlo	<< ";"
					<< F2_cpp_nnlo		<< std::endl;
	}
	std::cout << std::endl;
	fileout.close();
	Pdf::destroy();
}
