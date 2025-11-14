#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>

#include "LHAPDF/LHAPDF.h"
#include "structfunc.h"
#include "fortransymbols.h"

int main()	{
	Pdf::initialize("ABMP16_3_nnlo", 0);
	Pdf::setSampling(SAMPLINGMETHOD::fromOPENQCDRAD);
	Pdf::printLHAPDFinfo();

	///// check interpolation and tabulation for 1 specific grid knot
	// int ieta = 0,
	// 	ichi = 1;
	// double	logeta = ch_g_2_logetalist[ieta],
	// 		logchi	= ch_g_2_logchilist[ichi];
	// double	eta = std::pow(10,logeta),
	// 		chi = std::pow(10,logchi);
	// double sclca;
	// sclca_(&eta, &chi, &sclca);
	// std::cout << logeta << std::endl;
	// std::cout << logchi << std::endl;
	// std::cout << chL_g_2_0_A_table[ichi][ieta] << std::endl;
	// std::cout << sclca << std::endl;
	// std::cout << chL_g_2_0_A_interp(eta,chi) << std::endl;
	// std::cout << chL_g_2_0_A_interper(logeta,logchi) << std::endl;
	// std::cout << myInterp2D(
	// 	logchi,
	// 	logeta,
	// 	ch_g_2_logchilist,Nchi,
	// 	ch_g_2_logetalist,Neta,
	// 	chL_g_2_0_A_table[0]
	// ) << std::endl;

	// logeta = 2.523;
	// logchi = 2.523;
	// eta = std::pow(10,logeta);
	// chi = std::pow(10,logchi);
	// sclca_(&eta, &chi, &sclca);
	// std::cout << sclca << std::endl;
	// std::cout << chL_g_2_0_A_interp(eta,chi) << std::endl;
	// std::cout << chL_g_2_0_A_interper(logeta,logchi) << std::endl;
	// std::cout << myInterp2D(
	// 	logchi,
	// 	logeta,
	// 	ch_g_2_logchilist,Nchi,
	// 	ch_g_2_logetalist,Neta,
	// 	chL_g_2_0_A_table[0]
	// ) << std::endl;
	// return 0;

	qcdpar_.cf			= 4./3.;

	double	logetamin	= -6.0,
			logetamax	= 6.0,
			logchimin	= -3.0,
			logchimax	= 5.0;

	/// The hard coded grid in openQCD rad uses Neta=73 and Nchi=49,
	/// so in order to approximately include these points choose
	/// Neta = 1 + n*72, Nchi = 1 + m*48,
	/// where n,m specify the number of sample points between grid knots.
	/// The grid knots might not align perfectly, since the grid points 
	/// taken from openQCDrad are given by hard coded decimals
	/// e.g. -5.83333333d0 instead of -35d0/6d0.
	int	Neta	=	1+72,
		Nchi	=	1+48;


	/// output formatting
	const int PREC = 7;
	const int WIDTH = 23;

	std::ofstream fileout("output.dat");
	fileout << std::scientific << std::setprecision(PREC);

	int J = 0;
	for(int i = 0; i < Neta * Nchi; i++)	{
		int ieta	= i%Neta;
		int ichi	= i/Neta;
		while(J < (int)(100*(double)i/(double)(Neta*Nchi-1))) {
			std::cout << "\33[2K\r" << std::flush;
			std::cout << "[\33[32m";
			for(int j = 0; j < 100; j++)	{
				if(j < (int)(100*(double)i/(double)(Neta*Nchi-1)))	{
					std::cout << "\u2589";
				} else	{
					std::cout << "\u2591";
				}
			}
			J++;
			std::cout << "\33[0m]" << J << "%" << std::flush;
		}

		double logeta	= logetamin + (double)ieta/(double)(Neta-1) * (logetamax - logetamin);
		double logchi	= logchimin + (double)ichi/(double)(Nchi-1) * (logchimax - logchimin);

		double eta		= std::pow(10, logeta);
		double chi		= std::pow(10, logchi);

		////// gluon, longitudinal
		// double sclca;
		// double sclcf;
		// double sclbar;

		// sclca_(&eta, &chi, &sclca);
		// sclcf_(&eta, &chi, &sclcf);
		// sclbar_(&eta, &chi, &sclbar);

		// fileout		<< eta							<< ";" 	/// 0
		// 			<< chi							<< ";"	/// 1
		// 			///
		// 			<< clnlog_(&eta,&chi)			<< ";"	/// 2
		// 			<< asymp_l_(&chi)				<< ";"	/// 3
		// 			<< thresha_l_(&eta, &chi)		<< ";"	/// 4
		// 			<< threshf_l_(&eta, &chi)		<< ";"	/// 5
		// 			<< sclca						<< ";"	/// 6
		// 			<< sclcf						<< ";"	/// 7
		// 			///
		// 			<< clnlobarg_(&eta,&chi)		<< ";"	/// 8
		// 			<< asympbar_l_(&chi)			<< ";"	/// 9
		// 			<< threshbar_l_(&eta, &chi)		<< ";"	/// 10
		// 			<< sclbar						<< ";"	/// 11
		// 			///
		// 			<< chL_g_2_0(eta,chi)			<< ";"	/// 12
		// 			<< chL_g_2_0_asympG(eta,chi)	<< ";"	/// 13
		// 			<< chL_g_2_0_A_asympE(eta,chi)	<< ";"	/// 14
		// 			<< chL_g_2_0_F_asympE(eta,chi)	<< ";"	/// 15
		// 			<< chL_g_2_0_A_interp(eta,chi)	<< ";"	/// 16
		// 			<< chL_g_2_0_F_interp(eta,chi)	<< ";"	/// 17
		// 			///
		// 			<< chL_g_2_1(eta,chi)			<< ";"	/// 18
		// 			<< chL_g_2_1_asympG(eta,chi)	<< ";"	/// 19
		// 			<< chL_g_2_1_A_asympE(eta,chi)	<< ";"	/// 20
		// 			<< chL_g_2_1_interp(eta,chi)	<< std::endl;	/// 21

		////// gluon, transversal
		// double sctca;
		// double sctcf;
		// double sctbar;

		// sctca_(&eta, &chi, &sctca);
		// sctcf_(&eta, &chi, &sctcf);
		// sctbar_(&eta, &chi, &sctbar);

		// fileout		<< eta							<< ";" 	/// 0
		// 			<< chi							<< ";"	/// 1
		// 			///
		// 			<< ctnlog_(&eta,&chi)			<< ";"	/// 2
		// 			<< asymp_t_(&chi)				<< ";"	/// 3
		// 			<< thresha_t_(&eta, &chi)		<< ";"	/// 4
		// 			<< threshf_t_(&eta, &chi)		<< ";"	/// 5
		// 			<< sctca						<< ";"	/// 6
		// 			<< sctcf						<< ";"	/// 7
		// 			///
		// 			<< ctnlobarg_(&eta,&chi)		<< ";"	/// 8
		// 			<< asympbar_t_(&chi)			<< ";"	/// 9
		// 			<< threshbar_t_(&eta, &chi)		<< ";"	/// 10
		// 			<< sctbar						<< ";"	/// 11
		// 			///
		// 			<< chT_g_2_0(eta,chi)			<< ";"	/// 12
		// 			<< chT_g_2_0_asympG(eta,chi)	<< ";"	/// 13
		// 			<< chT_g_2_0_A_asympE(eta,chi)	<< ";"	/// 14
		// 			<< chT_g_2_0_F_asympE(eta,chi)	<< ";"	/// 15
		// 			<< chT_g_2_0_A_interp(eta,chi)	<< ";"	/// 16
		// 			<< chT_g_2_0_F_interp(eta,chi)	<< ";"	/// 17
		// 			///
		// 			<< chT_g_2_1(eta,chi)			<< ";"	/// 18
		// 			<< chT_g_2_1_asympG(eta,chi)	<< ";"	/// 19
		// 			<< chT_g_2_1_A_asympE(eta,chi)	<< ";"	/// 20
		// 			<< chT_g_2_1_interp(eta,chi)	<< std::endl;	/// 21

		//// QUARK
		double schql;
		double schqt;
		double sclql;
		double sclqt;
		double sqlbar;
		double sqtbar;

		schql_(&eta, &chi, &schql);
		schqt_(&eta, &chi, &schqt);
		sclql_(&eta, &chi, &sclql);
		sclqt_(&eta, &chi, &sclqt);
		sqlbar_(&eta, &chi, &sqlbar);
		sqtbar_(&eta, &chi, &sqtbar);

		fileout		<< eta							<< ";" 	/// 0
					<< chi							<< ";"	/// 1
					///
					<< clnloq_(&eta,&chi)			<< ";"	/// 2
					<< schql						<< ";"	/// 3
					<< ctnloq_(&eta,&chi)			<< ";"	/// 4
					<< schqt						<< ";"	/// 5
					<< dlnloq_(&eta,&chi)			<< ";"	/// 6
					<< sclql						<< ";"	/// 7
					<< dtnloq_(&eta,&chi)			<< ";"	/// 8
					<< sclqt						<< ";"	/// 9
					///
					<< clnlobarq_(&eta,&chi)		<< ";"	/// 10
					<< sqlbar						<< ";"	/// 11
					<< ctnlobarq_(&eta,&chi)		<< ";"	/// 12
					<< sqtbar						<< ";"	/// 13
					///
					<< chL_q_2_0_Hcoupl(eta,chi)		<< ";"	/// 14
					<< chL_q_2_0_Hcoupl_interp(eta,chi)	<< ";"	/// 15
					<< chT_q_2_0_Hcoupl(eta,chi)		<< ";"	/// 16
					<< chT_q_2_0_Hcoupl_interp(eta,chi)	<< ";"	/// 17
					<< chL_q_2_0_Lcoupl(eta,chi)		<< ";"	/// 18
					<< chL_q_2_0_Lcoupl_interp(eta,chi)	<< ";"	/// 19
					<< chT_q_2_0_Lcoupl(eta,chi)		<< ";"	/// 20
					<< chT_q_2_0_Lcoupl_interp(eta,chi)	<< ";"	/// 21
					///
					<< chL_q_2_1_Hcoupl(eta,chi)		<< ";"	/// 22
					<< chL_q_2_1_Hcoupl_interp(eta,chi)	<< ";"	/// 23
					<< chT_q_2_1_Hcoupl(eta,chi)		<< ";"	/// 24
					<< chT_q_2_1_Hcoupl_interp(eta,chi)	<< std::endl; /// 25
	}
	std::cout << std::endl;
	fileout.close();
	Pdf::destroy();


	///// check interpolation and tabulation for 1 specific grid knot
	// int ieta = 1,
	// 	ichi = 1;
	// double	logeta = ch_g_2_logetalist[ieta],
	// 		logchi	= ch_g_2_logchilist[ichi];
	// double	eta = std::pow(10,logeta),
	// 		chi = std::pow(10,logchi);
	// double sclca;
	// sclca_(&eta, &chi, &sclca);
	// std::cout << logeta << std::endl;
	// std::cout << logchi << std::endl;
	// std::cout << chL_g_2_0_A_table[ichi][ieta] << std::endl;
	// std::cout << sclca << std::endl;
	// std::cout << chL_g_2_0_A_interper(logeta,logchi) << std::endl;
	// std::cout << chL_g_2_0_A_interp(eta,chi) << std::endl;
}
