#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#include <chrono>
#include <cstdlib>
#include <ctime>

#include "LHAPDF/LHAPDF.h"

#define NXPLIM 201
#define NXMLIM 201
#define NSPLIM 241
#define NSMLIM 121
#define NFLIM 11

extern "C"	{
	struct gridset_type {
		double	delx1,
				delx2,
				delxp,
				dels1[8],
				dels2[8],
				xlog1,
				xlog2,
				x1,
				q2ini[8],
				q2min,
				q2max,
				xbmin,
				xbmax;
		int		nxmgrid,
				nxpgrid,
				nspgrid,
				nsmgrid,
				khalf;
	};
	extern gridset_type gridset_;

	struct gridtype {
		double 	xgrid[(NXMLIM+NXPLIM+1)], 
				sgrid[(NSMLIM+NSPLIM+1)*8],
				y[8*(NFLIM+1)*(NXMLIM+NXPLIM+1)*(NSMLIM+NSPLIM+1)],
				yhalf[2*NFLIM*(NXMLIM+NXPLIM+1)],
				pgrid[2*8*4*NFLIM*NFLIM*(NXMLIM+NXPLIM+1)],
				xpgrid[(NXMLIM+NXPLIM+1)];
	};
	extern gridtype grid_;

	void initgridconst_();
	void mypdffillgrid_witharg_(const char* arg, int arg_len);
	double xqg_(int* iq, double* xx, double* q2, int* kp);
	void evolvepdf_(double* x, double* Q, double* f);
}

double xqg_(int iq, double xx, double q2, int kp)	{
	return xqg_(&iq, &xx, &q2, &kp);
}

/**
 * @brief
 */
int main(int argc, char** argv)	{
	/// output formatting
	const int PREC = 15;
	const int WIDTH = 22;

	///
	std::string pdfset = "ABMP16_3_nnlo";
	const LHAPDF::PDF* currentpdf = LHAPDF::mkPDF(pdfset, 0);

	/// initialization of openQCDrad grid variables
	initgridconst_();
	mypdffillgrid_witharg_(pdfset.c_str(), pdfset.length());

	/// print the openQCDrad sampling points where the interpolation is exact
	std::cout << std::setprecision(PREC) << std::endl;
	for(int n = NXMLIM-gridset_.nxmgrid; n <= NXMLIM+gridset_.nxpgrid; n++)	{
		std::cout << std::setw(WIDTH) << grid_.xgrid[n] << std::endl;
	}

	const int 		NFINE	= 10*(gridset_.nxmgrid + gridset_.nxpgrid + 1);
	const double 	xmin	= grid_.xgrid[NXMLIM-gridset_.nxmgrid],
					xmax	= grid_.xgrid[NXMLIM+gridset_.nxpgrid];
	const double 	dx		= (xmax-xmin)/NFINE;
	// const double	q2		= 100;
	const double 	q2min	= 1.0,
					q2max	= 1e5;

	std::ofstream fileout("somepdf.dat");
	fileout << std::setprecision(PREC);
	fileout << std::setw(WIDTH) << "q2" 
			<< std::setw(WIDTH) << ";"
			<< std::setw(WIDTH) << "x" 
			<< std::setw(WIDTH) << ";"
			<< std::setw(WIDTH) << "xqg"
			<< std::setw(WIDTH) << ";"
			<< std::setw(WIDTH) << "evolvePDF"
			<< std::setw(WIDTH) << ";"
			<< std::setw(WIDTH) << "pdf->xfxQ2"
			<< std::endl;
	for(int iq = 0; iq <= 5; iq++)	{
		double q2 = 1.0 * std::pow(10,iq);
		for(int n = 0; n < NFINE; n++)	{
			double x = xmin + n * dx;
			double q = std::sqrt(q2);
			double farr[13];			///< parton flavors (-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6)
										///					(tb,bb,cb,sb,ub,db,g,d,u,s,c,b,t)
			evolvepdf_(&x, &q, farr);

			fileout 	<< std::setw(WIDTH) << q2
						<< std::setw(WIDTH) << ";"
						<< std::setw(WIDTH) << x 
						<< std::setw(WIDTH) << ";"
						<< std::setw(WIDTH) << xqg_(1, x, q2, 0)			///< openQCDrad interpolated values
																			///  openQCDrad uses 	(0,		1,2, 3,4, 5,...)
																			///						(alphas,g,d,db,u,ub,...)
						<< std::setw(WIDTH) << ";"
						<< std::setw(WIDTH) << farr[6]						///< LHAPDF fortran interface
						<< std::setw(WIDTH) << ";"
						<< std::setw(WIDTH) << currentpdf->xfxQ2(0,x, q2)	///< LHAPDF C++ interface
						<< std::endl;
		}
	}
	fileout.close();



	/// Performance
	std::srand(std::time(NULL));
	const int NRUN 		= 1e6;
	double	total_xqg	= 0.0,
			total_xfQ2	= 0.0,
			total_evPDF	= 0.0;
	auto start 			= std::chrono::high_resolution_clock::now();
	auto stop 			= std::chrono::high_resolution_clock::now();
	auto duration 		= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	for(int i = 0; i < NRUN; i++)	{
		double r1	= (double)std::rand()/(double)RAND_MAX;
		double r2	= (double)std::rand()/(double)RAND_MAX;

		double x	= xmin + (xmax-xmin)*r1;
		double q2	= q2min + (q2max-q2min)*r2;
		double q = std::sqrt(q2);

		start 		= std::chrono::high_resolution_clock::now();
		xqg_(1, x, q2, 0);
		stop 		= std::chrono::high_resolution_clock::now();
		duration	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		total_xqg	+= (double)duration.count();

		start 		= std::chrono::high_resolution_clock::now();
		currentpdf->xfxQ2(0,x, q2);
		stop 		= std::chrono::high_resolution_clock::now();
		duration	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		total_xfQ2	+= (double)duration.count();
		
		double farr[13];			///< parton flavors (-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6)
									///					(tb,bb,cb,sb,ub,db,g,d,u,s,c,b,t)
		start 		= std::chrono::high_resolution_clock::now();
		evolvepdf_(&x, &q, farr);
		stop 		= std::chrono::high_resolution_clock::now();
		duration	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
		total_evPDF	+= (double)duration.count();
	}

	std::cout	<< std::setw(WIDTH) << "xqg:"	
				<< std::setw(WIDTH) << total_xqg/(double)NRUN << " ns"
				<< std::endl
				<< std::setw(WIDTH) << "xfQ2:"	
				<< std::setw(WIDTH) << total_xfQ2/(double)NRUN << " ns"
				<< std::endl
				<< std::setw(WIDTH) << "evolvePDF:"	
				<< std::setw(WIDTH) << total_evPDF/(double)NRUN << " ns"
				<< std::endl;

	delete currentpdf;

	return 0;
}