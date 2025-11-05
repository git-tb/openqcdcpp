#include <iostream>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <fstream>

#include "LHAPDF/LHAPDF.h"
#include "hoppet.h"
#include "structfunc.h"
#include "fortransymbols.h"

/// @brief quick redifiniton such that we can pass literal constants (e.g. 1.234) as parameters
double f2qcd_(int nb, int nt, int ni, double xb, double q2) 	{
	return f2qcd_(&nb, &nt, &ni, &xb, &q2); ///< as defined in fortransymbols.h
}

/// some LHAPDF code needed by hoppet
LHAPDF::PDF *pdf = nullptr;
 
void lhapdf_interface(const double & x, const double & Q, double * res)  {
  std::vector<double> data(13, 0.0); // Pre-allocate and zero;
  pdf->xfxQ2(x, Q*Q, data);
  std::copy(data.begin(), data.end(), res); // Fast copy to the output array
};

void load_lhapdf_assign_hoppet(const std::string & pdfname, int imem=0) {
  // Start by loading a PDF set from LHAPDF
  pdf = LHAPDF::mkPDF(pdfname, imem);

  // Now let's access some basic information about the PDF to set up hoppet
  int nloop = pdf->orderQCD() + 1; // LHAPDF is zero indexed
  double xmin = pdf->xMin();
  double xmax = pdf->xMax();
  double Qmin = pdf->qMin();
  double Qmax = pdf->qMax();
  double mc = pdf->quarkMass(4);
  double mb = pdf->quarkMass(5);
  double mt = pdf->quarkMass(6);
  if(!pdf->hasFlavor(6)) mt = 2.0*Qmax; // If no top mass defined, use 173 GeV
 
  std::cout << "LHAPDF set: " << pdfname << " loaded successfully" << std::endl;
 
  // Now let us define some hoppet specific parameters. These are
  // typical values, and should guarantee similar accuracy as can be
  // expected from LHAPDF
  double ymax = static_cast<float>(std::ceil(log(1.0/xmin))); // To get a nice value of ymax that can contain the full LHAPDF grid
  double dy = 0.05;
  double dlnlnQ = dy/4.0;
  if(ymax > 15.0){
    dlnlnQ = dy/8.0; // for large ymax we need a finer grid in Q
  }
  int order = -6; // Default
  int yorder = 2; // Quadratic interpolation in y
  int lnlnQorder = 2; // Quadratic interpolation in lnlnQ
 
  // Print all the relevant parameters that we have passed to hoppet
  std::cout << "Hoppet starting with:" << std::endl;
  std::cout << " ymax:       " << ymax << std::endl;
  std::cout << " dy:         " << dy << std::endl;
  std::cout << " Qmin:       " << Qmin << std::endl;
  std::cout << " Qmax:       " << Qmax << std::endl;
  std::cout << " dlnlnQ:     " << dlnlnQ << std::endl;
  std::cout << " nloop:      " << nloop << std::endl;
  std::cout << " order:      " << order << std::endl;
  std::cout << " yorder:     " << yorder << std::endl;
  std::cout << " lnlnQorder: " << lnlnQorder << std::endl;
 
//   hoppetSetPoleMassVFN(mc,mb,mt); // set the pole masses
  hoppetSetFFN(3); /// or rather use fixed flavor number scheme for abmp16_3_nnlo
  hoppetSetYLnlnQInterpOrders(yorder, lnlnQorder); // Set the interpolation orders
  hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, hoppet::factscheme_MSbar); // Start hoppet
 
  // Now we fill the hoppet grid using the LHAPDF grid directly,
  // rather than evolving ourselves
  double Q0 = Qmin;
  hoppetSetCoupling(pdf->alphasQ(Q0), Q0, nloop);
  hoppetAssign(lhapdf_interface);
 
  // If instead we want to evolve the PDF with hoppet starting from
  // some low scale Q0 (>= Qmin) make a call to hoppetEvolve instead
  // of hoppetAssign
  //hoppetEvolve(pdf->alphasQ(Q0), Q0, nloop, 1.0, lhapdf_interface, Q0);
 
}

//----------------------------------------------------------------------
int main () {
	std::string pdfname = "ABMP16_3_nnlo";
	int imem = 0;

	/// set up HOPPET
	load_lhapdf_assign_hoppet(pdfname, imem);
	int nflav = 3;
	int order_max = 4; ///< 1=LO, 2=NLO, 3=NNLO, 4=N3LO
	int sc_choice = hoppet::scale_choice_Q;
	double zmass = 91.1876;
	double wmass = 80.377;
	bool param_coefs = true;
	hoppetStartStrFctExtended(order_max, nflav, sc_choice, /*mu_fixed(if needed)=*/zmass,param_coefs,wmass,zmass);
	double xmuR   = 1.0;
	double xmuF   = 1.0;	
	hoppetInitStrFct(order_max, /*separate_orders=*/true, xmuR, xmuF);
	
	/// set up openQCD++
	Pdf::initialize(pdfname, imem);
	PRECISION::EPSABS.set(1e-5);
	PRECISION::EPSREL.set(1e-5);
	PRECISION::ITER.set(1000);
	QCDORDER::F2ORDER.set(order_max-1);
	APPROX::LEVEL.set(APPROX::APPR2);

	/// set up openQCDrad
	forpreccontrol_.nf2qcd1 = 5;
	forpreccontrol_.nf2qcd2 = 5;
	foralpsrenorm_.kordf2_ 	= order_max-1;
	qcdpar_.nfc 			= 3;
	qcdpar_.cf				= 4./3.;
	qcdpar_.qsum[0]			= 1./9.;
	qcdpar_.qsum[1]			= 5./9.;
	qcdpar_.qsum[2]			= 6./9.;
	qcdpar_.qsum[3]			= 10./9.;
	qcdpar_.qsum[4]			= 12./9.;
	qcdpar_.qsum[5]			= 16./9.;

	/// comparison of outputs
	/// first, x and Q range to loop over
	double xmin		= 1e-5;
	double xmax		= 0.999;
	double Q2min	= 2;
	double Q2max	= 1e6;
	int Nx			= 20;
	int NQ2			= 500;
	
	double lnxmin	= std::log(xmin);
	double lnxmax	= std::log(xmax);
	double lnQ2min	= std::log(Q2min);
	double lnQ2max	= std::log(Q2max);

	/// ...and some preparation of output to disk.
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    std::stringstream timestamp_sstr;
    timestamp_sstr << std::put_time(&tm, "%Y%m%d_%H%M%S");
    std::string timestamp = timestamp_sstr.str();
	std::ofstream fileout_f2("f2_"+timestamp+".dat");
	std::ofstream fileout_pdf("pdf_"+timestamp+".dat");
	std::ofstream fileout_alphas("alphas_"+timestamp+".dat");
	fileout_f2	<< "# H = Hoppet" << std::endl
				<< "# C = my C++ program" << std::endl
				<< "# O = openQCDrad" << std::endl;
	fileout_pdf	<< "# H = Hoppet" << std::endl
				<< "# C = my C++ program" << std::endl
				<< "# O = openQCDrad" << std::endl;
	fileout_alphas	<< "# H = Hoppet" << std::endl
					<< "# C = my C++ program" << std::endl
					<< "# O = openQCDrad" << std::endl;
	fileout_f2 	<< "Q2;"
				<< "x;"
				<< "F2H@lo;"
				<< "F2C@lo;"
				<< "F2O@lo;"
				<< "F2H@nlo;"
				<< "F2C@nlo;"
				<< "F2O@nlo;"
				<< "F2H@nnlo;"
				<< "F2C@nnlo;"
				<< "F2O@nnlo;"
				<< std::endl;
	fileout_pdf	<< "Q2;"
				<< "x;"
				<< "tbH;tbC;tbO;"
				<< "bbH;bbC;bbO;"
				<< "cbH;cbC;cbO;"
				<< "sbH;sbC;sbO;"
				<< "ubH;ubC;ubO;"
				<< "dbH;dbC;dbO;"
				<< "gH;gC;gO;"
				<< "dH;dC;dO;"
				<< "uH;uC;uO;"
				<< "sH;sC;sO;"
				<< "cH;cC;cO;"
				<< "bH;bC;bO;"
				<< "tH;tC;tO"
				<< std::endl;
	fileout_alphas	<< "Q2;"
					<< "alphasH;"
					<< "alphasC;"
					<< "alphasO"
					<< std::endl;

	///
	///
	///
	for(int i = 0; i < Nx * NQ2; i++)	{
		int ix	= i%Nx;
		int iQ2	= i/Nx;
		// if(ix == 0 && i > 0) std::cout << std::endl;
		{
			std::cout << "\33[2K\r" << std::flush;
			// std::cout << (int)(100*(double)(i+1)/(double)(Nx*NQ2)) << "% done" << std::flush;
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

		double lnx	= lnxmin + (double)ix/(double)(Nx-1) * (lnxmax - lnxmin);
		double lnQ2	= lnQ2min + (double)iQ2/(double)(NQ2-1) * (lnQ2max - lnQ2min);
		double x	= std::exp(lnx);
		double Q2	= std::exp(lnQ2);
		double Q	= std::sqrt(Q2);
		

		///
		/// structure functions
		double	f2H_lo, f2H_nlo, f2H_nnlo,	///< Hoppet results
				f2C_lo, f2C_nlo, f2C_nnlo,	///< my C++ results
				f2O_lo, f2O_nlo, f2O_nnlo;	///< openQCDrad results

		/// H
		double StrFct_lo[14], StrFct_nlo[14], StrFct_nnlo[14];
		hoppetStrFctLO(x,Q,Q,Q,StrFct_lo);
		hoppetStrFctNLO(x,Q,Q,Q,StrFct_nlo);
		hoppetStrFctNNLO(x,Q,Q,Q,StrFct_nnlo);
		f2H_lo = StrFct_lo[hoppet::iF2EM];
		f2H_nlo = StrFct_nlo[hoppet::iF2EM] + f2H_lo;
		f2H_nnlo = StrFct_nnlo[hoppet::iF2EM] + f2H_nlo;
		/// C
		Pdf::setSampling(SAMPLINGMETHOD::fromOPENQCDRAD);
		QCDORDER::F2ORDER.set(0);
		f2C_lo = F2(x,Q2);
		QCDORDER::F2ORDER.set(1);
		f2C_nlo = F2(x,Q2);
		QCDORDER::F2ORDER.set(2);
		f2C_nnlo = F2(x,Q2);
		/// O
		foralpsrenorm_.kordf2_ = 0;
		f2O_lo = f2qcd_(3,1,22,x,Q2);
		foralpsrenorm_.kordf2_ = 1;
		f2O_nlo = f2qcd_(3,1,22,x,Q2);
		foralpsrenorm_.kordf2_ = 2;
		f2O_nnlo = f2qcd_(3,1,22,x,Q2);

		fileout_f2 	<< Q2 << ";"
					<< x << ";"
					<< f2H_lo << ";"
					<< f2C_lo << ";"
					<< f2O_lo << ";"
					<< f2H_nlo << ";"
					<< f2C_nlo << ";"
					<< f2O_nlo << ";"
					<< f2H_nnlo << ";"
					<< f2C_nnlo << ";"
					<< f2O_nnlo << ";"
					<< std::endl;

		///
		/// pdfs
		double hoppetpdf[13];
		hoppetEval(x,Q,hoppetpdf);

		fileout_pdf	<< Q2 << ";"
					<< x << ";";
		for(int id = -6; id <= 6; id++)	{
			double 	xfH,
					xfC,
					xfO;

			/// H
			xfH = hoppetpdf[id+6];
			/// C
			Pdf::setSampling(SAMPLINGMETHOD::fromLHAPDF);
			xfC = Pdf::xf(id, x, Q2);
			/// O
			Pdf::setSampling(SAMPLINGMETHOD::fromOPENQCDRAD);
			xfO = Pdf::xf(id, x, Q2);

			/// output
			fileout_pdf	<< xfH << ";"
						<< xfC << ";"
						<< xfO;
			if(id != 6) fileout_pdf << ";";
		}
		fileout_pdf << std::endl;

		///
		/// alphas
		if(ix != 0) continue;

		double	asH,
				asC,
				asO;

		/// H
		asH = hoppetAlphaS(Q);
		/// C
		Pdf::setSampling(SAMPLINGMETHOD::fromLHAPDF);
		asC = Pdf::alphas(Q2);
		/// O
		Pdf::setSampling(SAMPLINGMETHOD::fromOPENQCDRAD);
		asO = Pdf::alphas(Q2);

		/// output
		fileout_alphas	<< Q2 << ";"
						<< asH << ";"
						<< asC << ";"
						<< asO
						<< std::endl;		

	}
	std::cout << std::endl;

	fileout_f2.close();
	fileout_pdf.close();
	fileout_alphas.close();
}