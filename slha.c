#include "src/include.h"


/*--------------------------------------------------------*/
/* Calculation of the observables using a given SLHA file */
/*--------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[500];
	
	int test;
	double obs[Nobs_BKsll+1],omega;

  	if(argc<2) 
  	{ 
    		printf(" This program needs 1 parameter:\n"
           	"   name    name of the SLHA file\n");
      		exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%s",name);
  	}


	int filesOK=1;
	if(!filesOK) return 1;
 
	printf("\n");
	
	printf("SuperIso v4.0 - F. Mahmoudi\n\n");
	printf("SLHA input file\n\n");
 
	test=test_slha(name);
	
	if(test>0)
	{
		if(test==2) printf("WARNING: only tested in the MFV scenario!\n\n");

		printf("Observable\t\t\tValue\n\n");

		printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma_calculator(name));
		printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0_calculator(name));

		printf("BR(B0->K* gamma)\t\t%.3e\n",BR_BKstargamma_calculator(0,name));
		printf("BR(B+->K* gamma)\t\t%.3e\n\n",BR_BKstargamma_calculator(1,name));

		printf("BR(Bs->mu mu)\t\t\t%.3e\n",Bsmumu_calculator(name));
		printf("BR(Bs->mu mu)_untag\t\t%.3e\n",Bsmumu_untag_calculator(name));
		printf("BR(Bd->mu mu)\t\t\t%.3e\n\n",Bdmumu_calculator(name));
	
		printf("BR(B->K* mu mu)_low\t\t%.3e\n",BRobs_BKstarmumu_lowq2_calculator(name,obs));
		printf("AFB(B->K* mu mu)_low\t\t%.3e\n",obs[1]);
		printf("FL(B->K* mu mu)_low\t\t%.3e\n",obs[2]);
		printf("P1(B->K* mu mu)_low\t\t%.3e\n",obs[5]);
		printf("P2(B->K* mu mu)_low\t\t%.3e\n",obs[14]);
		printf("P4'(B->K* mu mu)_low\t\t%.3e\n",obs[17]);
		printf("P5'(B->K* mu mu)_low\t\t%.3e\n",obs[18]);
		printf("P6'(B->K* mu mu)_low\t\t%.3e\n",obs[19]);
		printf("P8'(B->K* mu mu)_low\t\t%.3e\n",obs[21]);
		printf("AI(B->K* mu mu)_low\t\t%.3e\n\n",AI_BKstarmumu_lowq2_calculator(name));
	
		printf("BR(B->K* mu mu)_high\t\t%.3e\n",BRobs_BKstarmumu_highq2_calculator(name,obs));
		printf("AFB(B->K* mu mu)_high\t\t%.3e\n",obs[1]);
		printf("FL(B->K* mu mu)_high\t\t%.3e\n",obs[2]);
		printf("P1(B->K* mu mu)_high\t\t%.3e\n",obs[5]);
		printf("P2(B->K* mu mu)_high\t\t%.3e\n",obs[14]);
		printf("P4'(B->K* mu mu)_high\t\t%.3e\n",obs[17]);
		printf("P5'(B->K* mu mu)_high\t\t%.3e\n",obs[18]);
		printf("P6'(B->K* mu mu)_high\t\t%.3e\n",obs[19]);
		printf("P8'(B->K* mu mu)_high\t\t%.3e\n",obs[21]);
		printf("AI(B->K* mu mu)_high\t\t%.3e\n\n",AI_BKstarmumu_highq2_calculator(name));

		printf("BR(B->Xs mu mu)_low\t\t%.3e\n",BRBXsmumu_lowq2_calculator(name));
		printf("BR(B->Xs mu mu)_high\t\t%.3e\n",BRBXsmumu_highq2_calculator(name));
		printf("q0^2(AFB(B->Xs mu mu)\t\t%.3e\n",A_BXsmumu_zero_calculator(name));
		printf("BR(B->Xs tau tau)_high\t\t%.3e\n\n",BRBXstautau_highq2_calculator(name));
	
		printf("BR(B->tau nu)\t\t\t%.3e\n",Btaunu_calculator(name));
		printf("R(B->tau nu)\t\t\t%.3e\n",RBtaunu_calculator(name));
		printf("BR(B->D tau nu)\t\t\t%.3e\n",BDtaunu_calculator(name));
		printf("BR(B->D tau nu)/BR(B->D e nu)\t%.3e\n",BDtaunu_BDenu_calculator(name));
		printf("BR(Ds->tau nu)\t\t\t%.3e\n",Dstaunu_calculator(name));
		printf("BR(Ds->mu nu)\t\t\t%.3e\n",Dsmunu_calculator(name));
		printf("BR(D->mu nu)\t\t\t%.3e\n",Dmunu_calculator(name));
		printf("BR(K->mu nu)/BR(pi->mu nu)\t%.3e\n",Kmunu_pimunu_calculator(name));
		printf("Rmu23(K->mu nu)\t\t\t%.3e\n\n",Rmu23_calculator(name));

		printf("a_muon\t\t\t\t%.3e\n\n",muon_gm2_calculator(name));

 		if(test==3)
		{ 	
			printf("theory_excluded\t\t\t%d\n",NMSSM_theory_excluded(name));
		}
		else printf("excluded_LEP/Tevatron_mass\t%d\n",excluded_mass_calculator(name));	

		printf("charged_LSP\t\t\t%d\n\n",charged_LSP_calculator(name));


		flha_generator(name,"output.flha");
		printf("output.flha generated\n\n");	
	}
	else if(test==-1) printf("Invalid point\n\n");
	else if(test==-2) printf("Model not yet implemented\n\n");
	else if(test==-3) printf("Invalid SLHA file\n\n");
	else if(test==-4) printf("SLHA file absent\n\n");
	
	return 1;
}
