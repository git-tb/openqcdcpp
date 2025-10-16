#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <ctime>

#include "LHAPDF/LHAPDF.h"
#include "structfunc.h"
#include "fortransymbols.h"
#include "chaplin.h"

int main()	{
	double 	xmin	=	5e-6,
			xmax	=	0.999;
	int		NX		=	10;

	/// output formatting
	const int PREC = 5;
	const int WIDTH = 30;

	std::cout	<< std::scientific << std::setprecision(PREC);
	std::cout	<< std::setw(WIDTH) << "x"
				<< std::setw(WIDTH) << "c2ns2int(apr/apr2/ext)"
				<< std::setw(WIDTH) << "c2ns2loc(apr/apr2/ext)"
				<< std::setw(WIDTH) << "c2ps2int(apr/apr2/ext)"
				<< std::setw(WIDTH) << "c2g2int(apr/apr2/ext)"
				<< std::endl;

	QCD::NF.set(3);

	for(int i = 0; i < NX; i++)	{
		double x = xmin + (double)i/(double)(NX-1)*(xmax-xmin);
		std::cout	<< std::setw(WIDTH) << x
					<< std::setw(WIDTH) << c2q_ns_2_0_plus_approx(x) + QCD::NF * c2q_ns_2_1_plus_approx(x) 
										 + c2q_ns_2_0_reg_approx(x) + QCD::NF * c2q_ns_2_1_reg_approx(x)
					<< std::setw(WIDTH) << c2q_ns_2_0_local_approx() + QCD::NF * c2q_ns_2_1_local_approx() + c2q_ns_2_0_localplus_approx(x) + QCD::NF * c2q_ns_2_1_localplus_approx(x)
					<< std::setw(WIDTH) << c2q_ps_2_0_reg_approx(x)
					<< std::setw(WIDTH) << c2g_2_0_reg_approx(x)
					<< std::endl
					<< std::setw(WIDTH) << " "
					<< std::setw(WIDTH) << c2q_ns_2_0_plus_approx2(x) + QCD::NF * c2q_ns_2_1_plus_approx2(x)
										 + c2q_ns_2_0_reg_approx2(x) + QCD::NF * c2q_ns_2_1_reg_approx2(x)
					<< std::setw(WIDTH) << c2q_ns_2_0_local_approx2() + QCD::NF * c2q_ns_2_1_local_approx2() + c2q_ns_2_0_localplus_approx2(x) + QCD::NF * c2q_ns_2_1_localplus_approx2(x)
					<< std::setw(WIDTH) << c2q_ps_2_0_reg_approx2(x)
					<< std::setw(WIDTH) << c2g_2_0_reg_approx2(x)
					<< std::endl
					<< std::setw(WIDTH) << " "
					<< std::setw(WIDTH) << c2q_ns_2_01_plus_exact(x) + c2q_ns_2_01_reg_exact(x)
					<< std::setw(WIDTH) << c2q_ns_2_01_local_exact() + c2q_ns_2_01_localplus_exact(x)
					<< std::setw(WIDTH) << c2q_ps_2_0_reg_exact(x) 
					<< std::setw(WIDTH) << c2g_2_0_reg_exact(x) 
					<< std::endl
					<< std::endl;
	}


	// ///
	// /// Performance
	// std::srand(std::time(NULL));
	// const int NRUN 		= 1e4;
	// double	time_appr	= 0.0, ///< approximation
	// 		time_ex		= 0.0; ///< exact calculation
	// auto start 			= std::chrono::high_resolution_clock::now();
	// auto stop 			= std::chrono::high_resolution_clock::now();
	// auto duration 		= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	// ///
	// std::cout << "Performance test" << std::endl;
	
	// for(int i = 0; i < NRUN; i++)	{
	// 	{
	// 		std::cout << "\33[2K\r" << std::flush;
	// 		std::cout << (int)(100*(double)(i+1)/(double)NRUN) << "% done" << std::flush;
	// 	}

	// 	double	r1	= (double)std::rand()/(double)RAND_MAX;
	// 	double	x	= xmin + r1*(xmax-xmin);

	// 	start 		= std::chrono::high_resolution_clock::now();
	// 	c2q_ns_2_0_plus_approx(x) + QCD::NF * c2q_ns_2_1_plus_approx(x) + c2q_ns_2_0_reg_approx(x) + QCD::NF * c2q_ns_2_1_reg_approx(x);
	// 	stop 		= std::chrono::high_resolution_clock::now();
	// 	duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	// 	time_appr 	+= (double)duration.count()/(double)NRUN;

	// 	start 		= std::chrono::high_resolution_clock::now();
	// 	c2q_ns_2_01_plus_exact(x) + c2q_ns_2_01_reg_exact(x);
	// 	stop 		= std::chrono::high_resolution_clock::now();
	// 	duration 	= std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	// 	time_ex 	+= (double)duration.count()/(double)NRUN;
	// }

	// std::cout	<< std::endl
	// 			<< std::setw(40) << "approximation runtime avg (ns):"
	// 			<< std::setw(10) << time_appr
	// 			<< std::endl
	// 			<< std::setw(40) << "exact calculation runtime avg (ns):"
	// 			<< std::setw(10) << time_ex
	// 			<< std::endl
	// 			<< std::endl;

	return 0;
}