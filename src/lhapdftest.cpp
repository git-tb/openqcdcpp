#include <iostream>
#include <boost/program_options.hpp>
#include <vector>
#include <string>
#include <fstream>

#include "LHAPDF/LHAPDF.h"

int main(int argc, char** argv) {
	///
	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);
	std::stringstream timestamp_sstr;
	timestamp_sstr << std::put_time(&tm, "%Y%m%d_%H%M%S");
	std::string timestamp = timestamp_sstr.str();

	/// declare cmd line options
	namespace po = boost::program_options;
	po::options_description desc("Program to test LHAPDF integration. Usage:");
	desc.add_options()
		("help", 			"produce help message")
		("pdfset",			po::value<std::string>()->required(),			"requested pdf set")	
		("pdfmem",			po::value<int>()->required(),					"requested pdf set member");	
	std::cout << desc << "\n";

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	/// process cmdline options
	if (vm.count("help"))    {
		std::cout << desc << "\n";
		return 1;
	}

	std::string pdfset 		= vm["pdfset"].as<std::string>();
	int pdfmem				= vm["pdfmem"].as<int>();

	///
	std::cout << "Available PDFsets:\n";
	for (const std::string& s : LHAPDF::availablePDFSets())
		std::cout << " " << s << ", ";
	std::cout << std::endl;



	/// actual pdf interfacing
	const LHAPDF::PDF* currentpdf = LHAPDF::mkPDF(pdfset, pdfmem);
	std::vector<int> pIDs = currentpdf->flavors();			/// < available parton IDs for this pdf set

	const double LOG10XMIN	= -5;
	const double LOG10XMAX	= -0.05;
	const int NX			= 10;
	const double DX			= (LOG10XMAX - LOG10XMIN)/(NX-1);

	const double LOG10Q2MIN	= 0;
	const double LOG10Q2MAX	= 5;
	const int NQ2			= 10;
	const double DQ2		= (LOG10Q2MAX - LOG10Q2MIN)/(NQ2-1);

	std::cout << pIDs.size() << " flavors available" << std::endl;
	for (int pID: pIDs)	{
		const std::string strpID = LHAPDF::lexical_cast<std::string>(pID);
		const std::string filename =  pdfset + "_" + std::to_string(pdfmem) + "_" + strpID + ".dat";
		std::ofstream fileout(filename.c_str());

		for(int ix = 0; ix < NX; ++ix) {
			const double log10x = (LOG10XMIN + ix*DX < -1e-3) ? LOG10XMIN + ix*DX : 0;
			const double x = std::pow(10, log10x);

			for(int iq2 = 0; iq2 < NQ2; ++iq2)	{
				const double log10q2 = LOG10Q2MIN + iq2*DQ2;
				const double q2 = std::pow(10, log10q2);

				const double xf = currentpdf->xfxQ2(pID, x, q2);

				fileout << x << " " << q2 << " " << xf << std::endl;
			}
		}

		fileout.close();
	}

	const int WIDTH = 20;	///< width per column for structured console output
	std::cout << std::scientific << std::setprecision(8);
	std::cout << std::setw(WIDTH) << "Q2" << std::setw(WIDTH) << "x" << std::setw(WIDTH) << "alphas";
	for (int pID: pIDs)	{	std::cout << std::setw(WIDTH) << "xf("+std::to_string(pID)+")";	}
	std::cout << std::endl;
	for(int iq2 = 0; iq2 < NQ2; ++iq2)	{
		const double log10q2 = LOG10Q2MIN + iq2*DQ2;
		const double q2 = std::pow(10, log10q2);
		for(int ix = 0; ix < NX; ++ix) {
			const double log10x = (LOG10XMIN + ix*DX < -1e-3) ? LOG10XMIN + ix*DX : 0;
			const double x = std::pow(10, log10x);

			std::cout << std::setw(WIDTH) << q2 << std::setw(WIDTH) << x;
			std::cout << std::setw(WIDTH) << currentpdf->alphasQ2(q2);
			for (int pID: pIDs)	{
				const double xf = currentpdf->xfxQ2(pID, x, q2);
				std::cout << std::setw(WIDTH) << xf;
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	// for(double log10q2 = LOG10Q2MIN; log10q2 <= LOG10Q2MAX; log10q2 += 0.2)	{
	// 	const double q2 = std::pow(10, log10q2);
	// 	std::cout << "alpha_s(" << std::setprecision(1) << std::fixed << std::sqrt(q2) << " GeV) = " << std::setprecision(5) << currentpdf->alphasQ2(q2) << std::endl;
	// }

	/// delete all pointers
	delete currentpdf;

	return 0;
}