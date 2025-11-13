#ifndef PDF_H
#define PDF_H

#include <iostream>
#include "LHAPDF/LHAPDF.h"

/// add here the PDF sampling functionality from openQCDrad such that we can compare predicitions
#include "fortransymbols.h"

double xqg(int iq, double xx, double q2, int kp);

/// PDGID codes for enhanced redability
/// gluon
const int 	G(21);	
/// down
const int	D(1);	
/// anti down
const int	DB(-1);	
/// up
const int	U(2);	
/// anti up
const int	UB(-2);	
/// strange
const int	S(3);	
/// anti strange
const int	SB(-3);	
/// charm
const int	C(4);	
/// anti charm
const int	CB(-4);	
/// bottom
const int	B(5);	
/// anti bottom
const int	BB(-5);	
/// bottom
const int	T(6);	
/// anti bottom
const int	TB(-6);	

/// @brief determine how PDFs are sampled
enum SAMPLINGMETHOD {
	fromLHAPDF,			///< use LHAPDFs sampling method xfxQ2
	fromOPENQCDRAD		///< use openQCDrads sampling method xqg
};

/**
 * @brief The idea of this class is to create a singleton, i.e. a class that can
 * only have 1 instance at runtime. Once the pdf is initialized via calls to LHAPDF
 * the use simply calls pdf(...). The advantage of having a singleton class for that
 * is that there is no global variable shared across files that needs to be created
 * somewhere, while still having a unique object that encapsulates all the pdf info.
 */
class Pdf {
    public:
		static void initialize(const std::string pdfname, const int pdfmem);

		static void destroy();

		static LHAPDF::PDF* get();

		static double alphas(double Q2);

		static double xf(int pID, double x, double Q2);

		/**
		 * @brief $\sum_i Q_i^2 x f(x,Q^2)$
		 * 
		 * @param x 
		 * @param Q2 
		 * @return double 
		 * 
		 * @todo check whether a requested flavor pid is actually available in the pdf set
		 * @todo variable flavor number
		 */
		static double xfiQi2sum(double x, double Q2);

		/**
		 * @brief x times singlet combination of quark PDFs
		 * 
		 * @param x 
		 * @param Q2 
		 * @return double 
		 * 
		 * @todo variable flavor number
		 */
		static double xfiSingletSum(double x, double Q2);

		static void printLHAPDFinfo();

		static void setSampling(SAMPLINGMETHOD method);
        
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
		inline static SAMPLINGMETHOD samplingmethod{fromLHAPDF};
};


#endif