#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <ctime>

#include "LHAPDF/LHAPDF.h"
#include "structfunc.h"

#define NXPLIM 201
#define NXMLIM 201
#define NSPLIM 241
#define NSMLIM 121
#define NFLIM 11

extern "C"	{
	struct commonblock1 {
		double	q20,
				q2rep,
				q2s,
				q20alphas,
				alphas0,
				alpsz,
				alpss,
				alpsc,
				alpsb,
				alpst,
				tscale,
				rscale,
				fscale,
				hqscale1,
				hqscale2;
		int 	nfeff,
				kordalps,
				kfeff,
				kordhq,
				kordf2_,
				kordfl,
				kordf3;
		double	almszl;	
	};
	extern commonblock1 foralpsrenorm_;

	struct commonblock2	{
		double	cf,
				cg,
				tr,
				qsum[6],
				qsum0[6],
				vfnth[6],
				vqu,
				aqu,
				vqd,
				aqd,
				vaq2u,
				vaq2d,
				vaq2sum[6],
				vaqsum,
				vlu,
				alu,
				vld,
				ald,
				val2u,
				val2d;
		int		nc,
				nf,
				nfe,
				nfc;
	};
	extern commonblock2 qcdpar_;

	struct commonblock3 {
		double	delder,
				alphastol;
		int		nmthq,
				nflhq,
				nf2hq,
				nf3hq,
				nflqcd,
				nf2qcd1,
				nf2qcd2,
				nf3qcd;
		double	omeint;
		int		lpcdint;
	};
	extern commonblock3 forpreccontrol_;

	void initgridconst_();
	void mypdffillgrid_witharg_(const char* arg, int arg_len);
	double f2qcd_(int* nb, int* nt, int* ni, double* xb, double* q2);
}

/// @brief quick redifiniton such that we can pass literal constants (e.g. 1.234) ass parameters
double f2qcd_(int nb, int nt, int ni, double xb, double q2) 	{
	return f2qcd_(&nb, &nt, &ni, &xb, &q2);
}

int main()	{
	Pdf::initialize("ABMP16_3_nnlo", 0);
	Pdf::printLHAPDFinfo();

	PRECISION::EPSABS.set(10e-6);
	PRECISION::EPSREL.set(10e-6);
	
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
	forpreccontrol_.nf2qcd1	= 100;
	forpreccontrol_.nf2qcd2	= 100;

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
			double F2_cpp_nnlo		= F2(x,Q2);

			std::cout	<< std::setw(WIDTH) << Q2
						<< std::setw(WIDTH) << x
						<< std::setw(WIDTH) << F2_fortran_nlo
						<< std::setw(WIDTH) << F2_fortran_nnlo
						<< std::endl
						<< std::setw(WIDTH)	<< " "
						<< std::setw(WIDTH)	<< " "
						<< std::setw(WIDTH) << F2_cpp_nlo
						<< std::setw(WIDTH) << F2_cpp_nnlo
						<< std::endl << std::endl;
		}
		std::cout << std::endl;
	}
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
// 	double	total_for	= 0.0, ///< fortran timing
// 			total_cpp	= 0.0; ///< C++ timing
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
// 		total_cpp		+= (double)duration.count();

// 		start 		= std::chrono::high_resolution_clock::now();
// 		f2qcd_(3,1,22,x,q2);
// 		stop 		= std::chrono::high_resolution_clock::now();
// 		duration	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
// 		total_for	+= (double)duration.count();
// 	}

// 	std::cout	<< std::endl
// 				<< std::setw(WIDTH) << "cpp:"	
// 				<< std::setw(WIDTH) << total_cpp/(double)NRUN << " ns"
// 				<< std::endl
// 				<< std::setw(WIDTH) << "fortran:"	
// 				<< std::setw(WIDTH) << total_for/(double)NRUN << " ns"
// 				<< std::endl;

// 	return 0;
// }