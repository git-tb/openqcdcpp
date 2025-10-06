#include <iostream>
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

	const double 	xmin	= grid_.xgrid[NXMLIM-gridset_.nxmgrid],
					xmax	= grid_.xgrid[NXMLIM+gridset_.nxpgrid];
	const double 	q2min	= 1.0,
					q2max	= 1e5;

	
	


	delete currentpdf;

	return 0;
}