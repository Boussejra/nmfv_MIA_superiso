#include "src/include.h"
#include "src/isajet.h"
#include "src/softsusy.h"
#include "src/suspect.h"
#include "src/spheno.h"


#define USE_SOFTSUSY /* to be commented if SOFTSUSY is unavailable */
//#define USE_ISAJET /* to be commented if ISAJET (>=7.80) is unavailable */
//#define USE_SUSPECT /* to be commented if SUSPECT is unavailable */
//#define USE_SPHENO /* to be commented if SPHENO is unavailable */


/*--------------------------------------------------------------------*/
/* Calculation of the observables, corresponding to a GMSB point generated by SOFTSUSY and/or ISAJET and/or SUSPECT and/or SPHENO */
/*--------------------------------------------------------------------*/

int main(int argc,char** argv)
{
	char name[50];
	double Lambda,Mmess,tanb,cGrav,sgnmu,mtop,mbot,alphas_mz,delta0;
	int N5;
	double obs[Nobs_BKsll+1],omega;

  	if(argc<5) 
  	{ 
    		printf(" This program needs 4 parameters:\n"
           	"   Lambda  scale of the SUSY breaking (10000-100000 GeV)\n"
           	"   Mmess   messenger mass scale > Lambda\n"
	   	"   N5      equivalent number of 5+5bar messenger fields\n"
           	"   tanb    tan(beta) \n");
    		printf(" Auxiliary parameters are:\n"
           	"   sgnmu   +/-1,  sign of Higgsino mass term (default 1)\n" 
	   	"   cGrav   >=1, ratio of the gravitino mass to its value for a breaking scale of F_m (not available for SUSPECT) \n"
           	"   mtop     top quark pole mass\n"
           	"   mbot     Mb(Mb) scale independent b-quark mass (not available for ISAJET)\n"
           	"   alphas_mz  strong coupling at MZ (not available for ISAJET)\n");                        
      		exit(1); 
  	} 
	else  
  	{  
		sscanf(argv[1],"%lf",&Lambda);
     		sscanf(argv[2],"%lf",&Mmess);
     		sscanf(argv[3],"%d",&N5);
     		sscanf(argv[4],"%lf",&tanb);
     		if(argc>5) sscanf(argv[5],"%lf",&sgnmu); else sgnmu=1.;
     		if(argc>6) sscanf(argv[6],"%lf",&cGrav); else cGrav=1.;   
     		if(argc>7) sscanf(argv[7],"%lf",&mtop); else mtop=173.34;   
     		if(argc>8) sscanf(argv[8],"%lf",&mbot); else mbot=4.18;
     		if(argc>9) sscanf(argv[9],"%lf",&alphas_mz); else alphas_mz=0.1184;
  	}

	
	int filesOK=1;
#ifdef USE_ISAJET
	if(!test_file(ISAJET)) 
	{
		printf("\"%s\" absent. Please check the ISAJET path or comment \"#define USE_ISAJET\" in gmsb.c\n",ISAJET);
		filesOK=0;
	}
#endif
#ifdef USE_SOFTSUSY
	if(!test_file(SOFTSUSY)) 
	{
		printf("\"%s\" absent. Please check the SOFTSUSY path or comment \"#define USE_SOFTSUSY\" in gmsb.c\n",SOFTSUSY);
		filesOK=0;
	}
#endif
#ifdef USE_SUSPECT
	if(!test_file(SUSPECT)) 
	{
		printf("\"%s\" absent. Please check the SUSPECT path or comment \"#define USE_SUSPECT\" in gmsb.c\n",SUSPECT);
		filesOK=0;
	}
#endif
#ifdef USE_SPHENO
	if(!test_file(SPHENO)) 
	{
		printf("\"%s\" absent. Please check the SPHENO path or comment \"#define USE_SPHENO\" in gmsb.c\n",SPHENO);
		filesOK=0;
	}
#endif
	if(!filesOK) return 1;


	if(Lambda>Mmess) printf("Lambda=%.0f must be smaller than Mmess=%.0f\n\n",Lambda,Mmess);
	else
	{

		if(!test_file("tmp")) system("mkdir tmp");
		chdir("tmp");

		printf("\n");
	
	printf("SuperIso v4.0 - F. Mahmoudi\n\n");

#ifdef USE_SOFTSUSY
		sprintf(name,"gmsb_softsusy%d.tmplha",getpid());
		softsusy_gmsb(Lambda, Mmess, tanb, N5, cGrav, sgnmu, mtop, mbot, alphas_mz, name);
		printf("---------------------------------------\n\n");
		printf("GMSB - SLHA file generated by SOFTSUSY\n\n");

		delta0=delta0_calculator(name);
		if(delta0 !=0.)
		{

			printf("Observable\t\t\tValue\n\n");

			printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma_calculator(name));
			printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0);

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

 			printf("excluded_LEP/Tevatron_mass\t%d\n",excluded_mass_calculator(name));
		
		
			flha_generator(name,"../output1.flha");
			printf("output1.flha generated\n\n");	
		}
		else printf("Invalid point\n\n");
		sprintf(name,"rm gmsb_softsusy%d.tmplha",getpid());
		system(name);		
#endif

#ifdef USE_ISAJET
		sprintf(name,"gmsb_isajet%d.tmplha",getpid());
		isajet_gmsb(Lambda, Mmess, tanb, N5, cGrav, sgnmu, mtop, name);
		printf("---------------------------------------\n\n");
		printf("GMSB - SLHA file generated by ISAJET\n\n");

		delta0=delta0_calculator(name);
		if(delta0 !=0.)
		{

			printf("Observable\t\t\tValue\n\n");

			printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma_calculator(name));
			printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0);

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

 			printf("excluded_LEP/Tevatron_mass\t%d\n",excluded_mass_calculator(name));
		
		
			flha_generator(name,"../output2.flha");
			printf("output2.flha generated\n\n");	
		}
		else printf("Invalid point\n\n");
		sprintf(name,"rm gmsb_isajet%d.tmplha",getpid());
		system(name);		
#endif
	
#ifdef USE_SUSPECT
		sprintf(name,"gmsb_suspect%d.tmplha",getpid());
		suspect_gmsb(Lambda, Mmess, tanb, N5, sgnmu, mtop, mbot, alphas_mz, name);
		printf("---------------------------------------\n\n");
		printf("GMSB - SLHA file generated by SUSPECT\n\n");

		delta0=delta0_calculator(name);
		if(delta0 !=0.)
		{

			printf("Observable\t\t\tValue\n\n");

			printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma_calculator(name));
			printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0);

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

 			printf("excluded_LEP/Tevatron_mass\t%d\n",excluded_mass_calculator(name));
		
		
			flha_generator(name,"../output3.flha");
			printf("output3.flha generated\n\n");	
		}
		else printf("Invalid point\n\n");
		sprintf(name,"rm gmsb_suspect%d.tmplha",getpid());
		system(name);		
#endif

#ifdef USE_SPHENO
		sprintf(name,"gmsb_spheno%d.tmplha",getpid());
		spheno_gmsb(Lambda, Mmess, tanb, N5, sgnmu, mtop, mbot, alphas_mz, name);
		printf("---------------------------------------\n\n");
		printf("GMSB - SLHA file generated by SPHENO\n\n");

		delta0=delta0_calculator(name);
		if(delta0 !=0.)
		{

			printf("Observable\t\t\tValue\n\n");

			printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma_calculator(name));
			printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0);

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

 			printf("excluded_LEP/Tevatron_mass\t%d\n",excluded_mass_calculator(name));

		
			flha_generator(name,"../output4.flha");
			printf("output4.flha generated\n\n");	
		}
		else printf("Invalid point\n\n");
		sprintf(name,"rm gmsb_spheno%d.tmplha",getpid());
		system(name);		
#endif
	}
	return 1;
}
