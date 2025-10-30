/// An example in C++ similar to the one found in examples_f90
///
#include "LHAPDF/LHAPDF.h"
#include "hoppet.h"
#include <iostream>
#include <cmath>
#include <cstdio>
#include <chrono>

// Global PDF pointer
LHAPDF::PDF *pdf = nullptr;
 
void lhapdf_interface(const double & x, const double & Q, double * res)  {
  std::vector<double> data(13, 0.0); // Pre-allocate and zero;
  pdf->xfxQ2(x, Q*Q, data);
  std::copy(data.begin(), data.end(), res); // Fast copy to the output array
};

void load_lhapdf_assign_hoppet(const std::string & pdfname, int imem=0) {
  // Start by loading a PDF set from LHAPDF
  auto t1 = std::chrono::high_resolution_clock::now();
  pdf = LHAPDF::mkPDF(pdfname, imem);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Time to load LHAPDF set: " << std::chrono::duration<double,std::milli>(t2 - t1).count() << " ms" << std::endl;
 
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
 
  hoppetSetPoleMassVFN(mc,mb,mt); // set the pole masses
  hoppetSetYLnlnQInterpOrders(yorder, lnlnQorder); // Set the interpolation orders
  t1 = std::chrono::high_resolution_clock::now();
  hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, hoppet::factscheme_MSbar); // Start hoppet
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Time to start hoppet: " << std::chrono::duration<double,std::milli>(t2 - t1).count() << " ms" << std::endl;
 
  // Now we fill the hoppet grid using the LHAPDF grid directly,
  // rather than evolving ourselves
  double Q0 = Qmin;
  hoppetSetCoupling(pdf->alphasQ(Q0), Q0, nloop);
  t1 = std::chrono::high_resolution_clock::now();
  hoppetAssign(lhapdf_interface);
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Time to fill hoppet grid from LHAPDF: " << std::chrono::duration<double,std::milli>(t2 - t1).count() << " ms" << std::endl;
 
  // If instead we want to evolve the PDF with hoppet starting from
  // some low scale Q0 (>= Qmin) make a call to hoppetEvolve instead
  // of hoppetAssign
  //hoppetEvolve(pdf->alphasQ(Q0), Q0, nloop, 1.0, lhapdf_interface, Q0);
 
}

//----------------------------------------------------------------------
int main () {
  std::string pdfname = "ABMP16_3_nnlo";
  int imem = 0;
  load_lhapdf_assign_hoppet(pdfname, imem);
  
  int nflav = -5;
  int order_max = 4;
  int sc_choice = hoppet::scale_choice_Q;
  double zmass = 91.1876;
  double wmass = 80.377;
  bool param_coefs = true;

  hoppetStartStrFctExtended(order_max, nflav, sc_choice,zmass,param_coefs,wmass,zmass);
    
  double asQ      = 0.35;
  double Q0       = sqrt(2.0);
  double muR_Q    = 1.0;
  double xmuR   = 1.0;
  double xmuF   = 1.0;

  hoppetInitStrFct(order_max,param_coefs, xmuR, xmuF);
    
  // output the results
  double pdf[13];
  double xvals[9]={1e-5,1e-4,1e-3,1e-2,0.1,0.3,0.5,0.7,0.9};
  double Q = 100;
  double StrFct[14];
  printf("                                Evaluating PDFs and structure functions at Q = %8.3f GeV\n",Q);
  printf("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon       F1γ         F2γ         F1Z         F2Z         F3Z\n");
  for (int ix = 0; ix < 9; ix++) {
    hoppetEval(xvals[ix], Q, pdf);
    hoppetStrFct(xvals[ix],Q,Q,Q,StrFct);
    printf("%7.1E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E\n",xvals[ix],
           pdf[6+2]-pdf[6-2], 
           pdf[6+1]-pdf[6-1], 
           2*(pdf[6-1]+pdf[6-2]),
           (pdf[6-4]+pdf[6+4]),
           pdf[6+0],
	   StrFct[hoppet::iF1EM],
	   StrFct[hoppet::iF2EM],
	   StrFct[hoppet::iF1Z],
	   StrFct[hoppet::iF2Z ],
	   StrFct[hoppet::iF3Z ]);
  }
  
}