/// An example in C++ showing how to replace LHAPDF interpolation with
/// (faster) hoppet interpolation
///
/// Usage (from CMake build directory):
///
///   example_cpp/fast_pdf_evaluation [PDFsetName]
///
#include "hoppet.h"
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
#include <iostream>
#include <cmath>
#include <cstdio>
#include <chrono>

using namespace std;
using namespace hoppet; // To access factscheme_MSbar
using namespace LHAPDF;


// Global PDF pointer
PDF *pdf = nullptr;

/// Interface to LHAPDF as needed by hoppetAssign
void lhapdf_interface(const double & x, const double & Q, double * res)  {
  vector<double> data(13, 0.0); // Pre-allocate and zero;
  pdf->xfxQ2(x, Q*Q, data);
  copy(data.begin(), data.end(), res); // Fast copy to the output array
};

/// Routine that loads an LHAPDF set, extracts some information from it
/// and transfers the PDF to hoppet. It needs the lhapdf_interface
/// defined above.
void load_lhapdf_assign_hoppet(const string & pdfname, int imem=0) {
  // Start by loading a PDF set from LHAPDF
  auto t1 = chrono::high_resolution_clock::now();
  pdf = mkPDF(pdfname, imem);
  auto t2 = chrono::high_resolution_clock::now();
  cout << "Time to load LHAPDF set: " << chrono::duration<double,milli>(t2 - t1).count() << " ms" << endl;

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

  cout << "LHAPDF set: " << pdfname << " loaded successfully" << endl;

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
  cout << "Hoppet starting with:" << endl;
  cout << " ymax:       " << ymax << endl;
  cout << " dy:         " << dy << endl;
  cout << " Qmin:       " << Qmin << endl;
  cout << " Qmax:       " << Qmax << endl;
  cout << " dlnlnQ:     " << dlnlnQ << endl;
  cout << " nloop:      " << nloop << endl;
  cout << " order:      " << order << endl;
  cout << " yorder:     " << yorder << endl;
  cout << " lnlnQorder: " << lnlnQorder << endl;

  hoppetSetPoleMassVFN(mc,mb,mt); // set the pole masses
  hoppetSetYLnlnQInterpOrders(yorder, lnlnQorder); // Set the interpolation orders
  t1 = chrono::high_resolution_clock::now();
  hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, factscheme_MSbar); // Start hoppet
  t2 = chrono::high_resolution_clock::now();
  cout << "Time to start hoppet: " << chrono::duration<double,milli>(t2 - t1).count() << " ms" << endl;

  // Now we fill the hoppet grid using the LHAPDF grid directly,
  // rather than evolving ourselves
  double Q0 = Qmin;
  hoppetSetCoupling(pdf->alphasQ(Q0), Q0, nloop);
  t1 = chrono::high_resolution_clock::now();
  hoppetAssign(lhapdf_interface);
  t2 = chrono::high_resolution_clock::now();
  cout << "Time to fill hoppet grid from LHAPDF: " << chrono::duration<double,milli>(t2 - t1).count() << " ms" << endl;

  // If instead we want to evolve the PDF with hoppet starting from
  // some low scale Q0 (>= Qmin) make a call to hoppetEvolve instead
  // of hoppetAssign
  //hoppetEvolve(pdf->alphasQ(Q0), Q0, nloop, 1.0, lhapdf_interface, Q0);

}

//----------------------------------------------------------------------
int main (int argc, char *argv[]) {
  string pdfname = "CT18NNLO";
  if (argc > 1) {
    pdfname = argv[1];
  }
  cout << "Using PDF set: " << pdfname << endl;
  int imem = 0;
  load_lhapdf_assign_hoppet(pdfname, imem);

  double x = 0.01;
  double Q = 13.0;

  const LHAPDF::GridPDF& gridpdf = * dynamic_cast<const LHAPDF::GridPDF*>(pdf);

  std::cout << std::setprecision(7) << std::scientific;
  for(auto& Q2: gridpdf.q2Knots()) {
    for(auto& x: gridpdf.xKnots()) {
		// Standard LHAPDF call
  		vector<double> lhapdf;
  		pdf->xfxQ2(x, Q2, lhapdf); // Get the PDF from LHAPDF

		// Equivalent hoppet call
  		double hoppetpdf[13];
  		hoppetEval(x, std::sqrt(Q2), hoppetpdf);

		std::cout 	<< std::setw(15) << std::sqrt(Q2)
					<< std::setw(15) << x
					<< std::setw(15) << (lhapdf[6+2] - hoppetpdf[6+2])/lhapdf[6+2];
		std::cout	<< std::endl;
	}
  }

  hoppetDeleteAll();

}



