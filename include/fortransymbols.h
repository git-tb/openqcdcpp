#ifndef FORTRANSYMBOLS_H
#define FORTRANSYMBOLS_H

/**
 * @brief Here we import some symbols from the old fortran code
 * in order to use and modify them in our code. COMMON blocks from
 * Fortran can be imported by mapping them to structs with the same
 * order of member variables as defined in the Fortran code.
 * 
 * Depending on the compiler the symbol names from the Fortran source
 * code are changed slightly in the corresponding object files that
 * we link to the C++ program, hence the trailing underscore "name_"
 * at the end of each variable name.
 * 
 * Check for the precise symbols names with
 * >>	nm libMyLib.so | grep symbolname
 *  
 */

#define NXPLIM 201
#define NXMLIM 201
#define NSPLIM 241
#define NSMLIM 121
#define NFLIM 11

extern "C"	{
	struct FORALPSRENORM_COMMON {
		double	q20,
				q2rep,
				q2s,
				q20alphas,
				alphas0,
				alpsz,
				alpss,
				alpsc,
				alpsb,
				alpst,
				tscale,
				rscale,
				fscale,
				hqscale1,
				hqscale2;
		int 	nfeff,
				kordalps,
				kfeff,
				kordhq,
				kordf2,
				kordfl,
				kordf3;
		double	almszl;	
	};
	extern FORALPSRENORM_COMMON foralpsrenorm_;

	struct QCDPAR_COMMON	{
		double	cf,
				cg,
				tr,
				qsum[6],
				qsum0[6],
				vfnth[6],
				vqu,
				aqu,
				vqd,
				aqd,
				vaq2u,
				vaq2d,
				vaq2sum[6],
				vaqsum,
				vlu,
				alu,
				vld,
				ald,
				val2u,
				val2d;
		int		nc,
				nf,
				nfe,
				nfc;
	};
	extern QCDPAR_COMMON qcdpar_;

	struct FORPRECCONTROL_COMMON {
		double	delder,
				alphastol;
		int		nmthq,
				nflhq,
				nf2hq,
				nf3hq,
				nflqcd,
				nf2qcd1,
				nf2qcd2,
				nf3qcd;
		double	omeint;
		int		lpcdint;
	};
	extern FORPRECCONTROL_COMMON forpreccontrol_;

	struct GRIDSET_COMMON {
		double	delx1,
				delx2,
				delxp,
				dels1[8],
				dels2[8],
				xlog1,
				xlog2,
				x1,
				q2ini[8],
				q2min,
				q2max,
				xbmin,
				xbmax;
		int 	nxmgrid,
				nxpgrid,
				nspgrid,
				nsmgrid,
				khalf;
	};
	extern GRIDSET_COMMON gridset_;

	struct FORF2CHARM_COMMON	{
		double	xb0,
				q2ss,
				rm2,
				qqs,
				rmu2,
				an,
				rq;
		int		nb0,
				nt0,
				ni0,
				nq0;
	};
	extern FORF2CHARM_COMMON forf2charm_;

	struct FORSCHEMEDEF_COMMON	{
		double	ddnlohq;
		int	msbarm,
				hqnons;
		double	bmsnfopt,
				bmsnnlo,
				vloop;
	};
	extern FORSCHEMEDEF_COMMON forschemedef_;	

	void initgridconst_();
	void mypdffillgrid_witharg_(const char* arg, int arg_len);
	double xqg_(int* iq, double* xx, double* q2, int* kp);
	double f2qcd_(int* nb, int* nt, int* ni, double* xb, double* q2);
	double f2charm_ffn_(double* xb, double* q2, int* nq);
	double f2charmi_(double* t);
	double clnlog_(double* eta, double* xi);
	double clnloq_(double* eta, double* xi);
	double dlnloq_(double* eta, double* xi);
	double ctnlog_(double* eta, double* xi);
	double ctnloq_(double* eta, double* xi);
	double dtnloq_(double* eta, double* xi);
	double clnlobarg_(double* eta, double* xi);
	double clnlobarq_(double* eta, double* xi);
	double ctnlobarg_(double* eta, double* xi);
	double ctnlobarq_(double* eta, double* xi);
	void sclca_(double* eta, double* xi, double* xsclca);
	void sctca_(double* eta, double* xi, double* xsclca);
	void sclcf_(double* eta, double* xi, double* xsclcf);
	void sctcf_(double* eta, double* xi, double* xsclcf);
	void schql_(double* eta, double* xi, double* xschql);
	void schqt_(double* eta, double* xi, double* xschqt);
	void sclql_(double* eta, double* xi, double* xsclql);
	void sclqt_(double* eta, double* xi, double* xsclqt);
	void sclbar_(double* eta, double* xi, double* xsclcbar);
	void sctbar_(double* eta, double* xi, double* xsclcbar);
	void sqlbar_(double* eta, double* xi, double* xsqlcbar);
	void sqtbar_(double* eta, double* xi, double* xsqlcbar);
	double asymp_l_(double* xi);
	double asymp_t_(double* xi);
	double asympbar_l_(double* xi);
	double asympbar_t_(double* xi);
	double thresha_l_(double* eta, double* xi);
	double thresha_t_(double* eta, double* xi);
	double threshf_l_(double* eta, double* xi);
	double threshf_t_(double* eta, double* xi);
	double threshbar_l_(double* eta, double* xi);
	double threshbar_t_(double* eta, double* xi);
}

#endif