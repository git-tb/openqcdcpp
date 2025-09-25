#ifndef PDF_H
#define PDF_H

#include <iostream>
#include "LHAPDF/LHAPDF.h"

/// @brief PDGID codes for enhanced redability
const int 	G(1),		/// < gluon
			D(1),		/// < down
			DB(-1),		/// < anti down
			U(2),		/// < up
			UB(-2),		/// < anti up
			S(3),		/// < strange
			SB(-3),		/// < anti strange
			C(4),		/// < charm
			CB(-4),		/// < anti charm
			B(5),		/// < bottom
			BB(-5);		/// < anti bottom


/**
 * @brief The idea of this class is to create a singleton, i.e. a class that can
 * only have 1 instance at runtime. Once the pdf is initialized via calls to LHAPDF
 * the use simply calls pdf(...). The advantage of having a singleton class for that
 * is that there is no global variable shared across files that needs to be created
 * somewhere, while still having a unique object that encapsulates all the pdf info.
 */
class Pdf {
    public:
		static void initialize(const std::string pdfname, const int pdfmem)	{
			delete lhapdfobject;
			lhapdfobject = LHAPDF::mkPDF(pdfname, pdfmem);
			Pdf::pdfname = pdfname;
			Pdf::pdfmem = pdfmem;
		}

		inline static LHAPDF::PDF* get() {
			if(not lhapdfobject) {
				std::cout << "Initialize PDF object before using it!" << std::endl;
				abort();
			}

			return lhapdfobject;
		}

		static void printLHAPDFinfo()	{
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
				// for(const auto key: lhapdfobject->info().keys()) std::cout << "@@@\t" << key << "\t\t" << lhapdfobject->info().get_entry(key) << std::endl;

				// std::cout << "@@@ available flavors in current pdf: ";
				// for(int pid : lhapdfobject->flavors())	
				// 	std::cout << pid << "(" << LHAPDF::to_str(pid) << "), ";
				// std::cout << std::endl;
			} else	{
				std::cout << "@@@ currently no PDFset is loaded!" << std::endl;
			}
		}
        
    private:
		/// Delete the copy constructor and assignment operator
        Pdf(const Pdf&) = delete;
        Pdf& operator=(const Pdf&) = delete;

        Pdf();
        ~Pdf();
		
		/// With "inline" we dont need to define the members elsewhere (in some .cpp file).
		inline static LHAPDF::PDF* lhapdfobject{nullptr};
		inline static std::string pdfname;
		inline static int pdfmem;
};

/**
 * @brief Computes $\sum_i Q_i^2 x f(x,Q^2)$
 * 
 * @param x 
 * @param Q2 
 * @return double 
 * 
 * @todo Check whether a requested flavor pid is actually available in the pdf set
 */
double xfiQi2sum(double x, double Q2)	{
	double result(0.0);
	result += 	4.0/9.0 * ( Pdf::get()->xfxQ2(U, x, Q2) + Pdf::get()->xfxQ2(UB, x, Q2) ) + \
			+	1.0/9.0 * ( Pdf::get()->xfxQ2(D, x, Q2) + Pdf::get()->xfxQ2(DB, x, Q2) ) + \
			+	1.0/9.0 * ( Pdf::get()->xfxQ2(S, x, Q2) + Pdf::get()->xfxQ2(SB, x, Q2) );
	return result;
}


#endif