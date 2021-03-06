#include "include.h"

double Dslnu(int gen, struct parameters* param)
/* computes the branching ratio of Ds -> ell nu */
{
	double Vcs=cabs(param->Vcs);

	double ml,lambda_l;
	if(gen==1) {ml=param->mass_e; lambda_l=param->lambda_l[1][1];}
	else if(gen==2) {ml=param->mass_mu; lambda_l=param->lambda_l[2][2];}
	else {ml=param->mass_tau; lambda_l=param->lambda_l[3][3];}

	if(param->SM==1) return param->m_Ds/8./pi*pow(param->Gfermi*Vcs*ml*param->f_Ds*(1.-ml*ml/param->m_Ds/param->m_Ds),2.)*param->life_Ds/hbar;
	
	double mc=running_mass(param->mass_c,param->mass_c,param->m_Ds,param->mass_top_pole,param->mass_b,param);
	double ms=running_mass(param->mass_s,2.,param->m_Ds,param->mass_top_pole,param->mass_b,param);

	if(param->THDM_model>0) return param->m_Ds/8./pi*pow(param->Gfermi*Vcs*ml*param->f_Ds*(1.-ml*ml/param->m_Ds/param->m_Ds)*(1.+param->m_Ds*param->m_Ds/param->mass_H/param->mass_H*(mc*param->lambda_u[2][2]-ms*param->lambda_d[2][2])*lambda_l/(ms+mc)),2.)*param->life_Ds/hbar;

	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);
	double epsilon0=-2./3.*alphas_MSOFT/pi*param->mu_Q/param->mass_gluino*H2(param->MqL2_Q*param->MqL2_Q/param->mass_gluino/param->mass_gluino,param->McR_Q*param->McR_Q/param->mass_gluino/param->mass_gluino);
	
	return param->m_Ds/8./pi*pow(param->Gfermi*Vcs*ml*param->f_Ds*(1.-ml*ml/param->m_Ds/param->m_Ds)*(1.+param->m_Ds*param->m_Ds/param->mass_H/param->mass_H*(mc-ms*param->tan_beta*param->tan_beta/(1.+epsilon0*param->tan_beta))/(ms+mc)),2.)*param->life_Ds/hbar;
}

/*--------------------------------------------------------------------*/

double Dslnu_calculator(int gen, char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(Ds -> ell nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Dslnu(gen,&param);
}

/*--------------------------------------------------------------------*/


double Dstaunu(struct parameters* param)
/* computes the branching ratio of Ds -> tau nu */
{
	return Dslnu(3,param);
}

/*--------------------------------------------------------------------*/

double Dstaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(Ds -> tau nu) */
{
	return Dslnu_calculator(3,name);
}

/*--------------------------------------------------------------------*/

double Dsmunu(struct parameters* param)
/* computes the branching ratio of Ds -> mu nu */
{
	return Dslnu(2,param);
}

/*--------------------------------------------------------------------*/

double Dsmunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(Ds -> mu nu) */
{
	return Dslnu_calculator(2,name);
}
