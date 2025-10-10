#ifndef FORTRANSYMBOLS_H
#define FORTRANSYMBOLS_H

#define NXPLIM 201
#define NXMLIM 201
#define NSPLIM 241
#define NSMLIM 121
#define NFLIM 11

extern "C"	{
	struct commonblock1 {
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
				kordf2_,
				kordfl,
				kordf3;
		double	almszl;	
	};
	extern commonblock1 foralpsrenorm_;

	struct commonblock2	{
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
	extern commonblock2 qcdpar_;

	struct commonblock3 {
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
	extern commonblock3 forpreccontrol_;

	void initgridconst_();
	void mypdffillgrid_witharg_(const char* arg, int arg_len);
	double xqg_(int* iq, double* xx, double* q2, int* kp);
	double f2qcd_(int* nb, int* nt, int* ni, double* xb, double* q2);
}

#endif