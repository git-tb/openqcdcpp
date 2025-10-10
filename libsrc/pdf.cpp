#include "pdf.h"
#include <iostream>
#include "LHAPDF/LHAPDF.h"
#include "fortransymbols.h"

double xqg(int iq, double xx, double q2, int kp)	{
	return xqg_(&iq, &xx, &q2, &kp);
}

/**
 * @brief The idea of this class is to create a singleton, i.e. a class that can
 * only have 1 instance at runtime. Once the pdf is initialized via calls to LHAPDF
 * the use simply calls pdf(...). The advantage of having a singleton class for that
 * is that there is no global variable shared across files that needs to be created
 * somewhere, while still having a unique object that encapsulates all the pdf info.
 */
void Pdf::initialize(const std::string pdfname, const int pdfmem)	{
	delete lhapdfobject;
	lhapdfobject = LHAPDF::mkPDF(pdfname, pdfmem);
	Pdf::pdfname = pdfname;
	Pdf::pdfmem = pdfmem;

	/// ...in case openQCDrad code should be used...
	initgridconst_();
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
	} else {
		return xqg(0, 0.1, Q2, pdfmem);
	}
}

double Pdf::xf(int pID, double x, double Q2)	{
	if(not lhapdfobject) {
		std::cout << "Initialize PDF object before using it!" << std::endl;
		abort();
	}

	if(samplingmethod == fromLHAPDF)	{
		return lhapdfobject->xfxQ2(pID, x, Q2);
	} else {
		if(pID == G) return xqg(1,x,Q2,pdfmem);
		int idx = 2 * std::abs(pID) + (pID < 0 ? 1 : 0);
		return xqg(idx, x, Q2, pdfmem);
	}
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