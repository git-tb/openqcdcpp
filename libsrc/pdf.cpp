#include "pdf.h"
#include <iostream>
#include "LHAPDF/LHAPDF.h"
#include "fortransymbols.h"

double xqg(int iq, double xx, double q2, int kp)	{
	return xqg_(&iq, &xx, &q2, &kp);
}

void Pdf::initialize(const std::string pdfname, const int pdfmem)	{
	delete lhapdfobject;
	lhapdfobject = LHAPDF::mkPDF(pdfname, pdfmem);
	Pdf::pdfname = pdfname;
	Pdf::pdfmem = pdfmem;

	/// ...in case openQCDrad code should be used...
	initgridconst_();
	/// I think there is some funny business going on in openQCDrads "initgridconst"
	/// function and that the sampling is strongly affected by some artificial
	/// threshold in the q2 grid. Let's check this: In openQCDrad this threshold is
	///	at Q2=9. Change it so some other value and see if the benchmarks become worse
	/// at this value.
	for(auto& v: gridset_.q2ini) v=12;
	mypdffillgrid_witharg_(pdfname.c_str(), pdfname.length());
}

void Pdf::destroy()	{
	delete lhapdfobject;
	lhapdfobject = NULL;
}


LHAPDF::PDF* Pdf::get() {
	if(not lhapdfobject) {
		std::cout << "Initialize PDF object before using it!" << std::endl;
		abort();
	}
	return lhapdfobject;
}

double Pdf::alphas(double Q2)	{
	if(not lhapdfobject) {
		std::cout << "Initialize PDF object before using it!" << std::endl;
		abort();
	}

	if(samplingmethod == fromLHAPDF)	{
		return lhapdfobject->alphasQ2(Q2);
	} else if(samplingmethod == fromOPENQCDRAD) {
		/// @todo This always samples the three flavor Pdf from openQCDrad
		return xqg(0, 0.1, Q2, /*kschemepdf=*/0);
	}

	std::cout << "WARNING: alphas sampling method is undefined!" << std::endl;
	return -1e10;
}

double Pdf::xf(int pID, double x, double Q2)	{
	if(not lhapdfobject) {
		std::cout << "Initialize PDF object before using it!" << std::endl;
		abort();
	}

	if(samplingmethod == fromLHAPDF)	{
		return lhapdfobject->xfxQ2(pID, x, Q2);
	} else if(samplingmethod == fromOPENQCDRAD) {
		if(pID == G || pID == 0) return xqg(1,x,Q2,pdfmem);
		int idx = 2 * std::abs(pID) + (pID < 0 ? 1 : 0);
		/// @todo This always samples the three flavor Pdf from openQCDrad
		return xqg(idx, x, Q2, /*kschemepdf=*/0);
	}

	std::cout << "WARNING: pdf sampling method is undefined!" << std::endl;
	return -1e10;
}

void Pdf::printLHAPDFinfo()	{
	std::cout << "@@@ LHAPDF paths: ";
	for (const std::string& p : LHAPDF::paths())
		std::cout << p << ", ";
	std::cout << std::endl;
	std::cout << "@@@ available PDFs: ";
	for (const std::string& s : LHAPDF::availablePDFSets())
		std::cout << s << ", ";
	std::cout << std::endl;
	
	if(lhapdfobject)	{
		std::cout << "@@@ currently loaded PDFset:\n" << pdfname << " (mem: " << pdfmem << ")" << std::endl;
		for(const auto &key: lhapdfobject->info().keys()) std::cout << "@@@\t" << key << "\t\t" << lhapdfobject->info().get_entry(key) << std::endl;
	} else	{
		std::cout << "@@@ currently no PDFset is loaded!" << std::endl;
	}
}

void Pdf::setSampling(SAMPLINGMETHOD method)	{
	samplingmethod = method;
}


double xfiQi2sum(double x, double Q2)	{
	double result(0.0);
	result += 	4.0/9.0 * ( Pdf::xf(U, x, Q2) + Pdf::xf(UB, x, Q2) ) + \
			+	1.0/9.0 * ( Pdf::xf(D, x, Q2) + Pdf::xf(DB, x, Q2) ) + \
			+	1.0/9.0 * ( Pdf::xf(S, x, Q2) + Pdf::xf(SB, x, Q2) );
	return result;
}

double xfiSingletSum(double x, double Q2)	{
	double result(0.0);
	result +=	( Pdf::xf(U, x, Q2) + Pdf::xf(UB, x, Q2) ) + \
			+	( Pdf::xf(D, x, Q2) + Pdf::xf(DB, x, Q2) ) + \
			+	( Pdf::xf(S, x, Q2) + Pdf::xf(SB, x, Q2) );
	return result;
}