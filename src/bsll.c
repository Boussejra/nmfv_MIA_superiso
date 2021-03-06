#include "include.h"

/*----------------------------------------------------------------------*/

double complex g_bsll(double z, double s)
{
	double z2=z*z;

	if(s==0.) return -4./9.*log(z2)+8./27.-4./9.;
	
	if(z==0.) return 8./27.-4./9.*(log(s)-I*pi);
	
	if(4.*z2<s) return -4./9.*log(z2)+8./27.+16./9.*z2/s -2./9.*sqrt(1.-4.*z2/s)*(2.+4.*z2/s)*(log((1.+sqrt(1.-4.*z2/s))/(1.-sqrt(1.-4.*z2/s)))-I*pi);
		
	else if(4.*z2>s) return -4./9.*log(z2)+8./27.+16./9.*z2/s -4./9.*sqrt(4.*z2/s-1.)*(2.+4.*z2/s)*atan(1./sqrt(4.*z2/s-1.));
	
	else return -4./9.*log(z2)+8./27.+16./9.*z2/s;
}

/*----------------------------------------------------------------------*/

double Rcchad(double s, struct parameters* param)
{
	int ie;
	double BR[7],GammaTot[7],GammaHad[7],mV[7];
	
	mV[1]=3.096917;
	BR[1]=5.93e-2;
	GammaHad[1]=0.08147e-3;
	GammaTot[1]=0.0929e-3;
	
	mV[2]=3.68609;
	BR[2]=7.7e-3;
	GammaHad[2]=0.29746e-3;
	GammaTot[2]=0.304e-3;
	
	mV[3]=3.77292;
	BR[3]=1.1e-5;
	GammaHad[3]=23.6e-3;
	GammaTot[3]=27.3e-3;
	
	mV[4]=4.039;
	BR[4]=1.4e-5;
	GammaHad[4]=52.e-3;
	GammaTot[4]=80.e-3;
	
	mV[5]=4.153;
	BR[5]=1.0e-5;
	GammaHad[5]=78.e-3;
	GammaTot[5]=103.e-3;
	
	mV[6]=4.421;
	BR[6]=1.1e-5;
	GammaHad[6]=43.e-3;
	GammaTot[6]=62.e-3;
	
	double Rccres=0.;
	for(ie=1;ie<=6;ie++) Rccres+=BR[ie]*GammaTot[ie]*GammaHad[ie]/(pow(s*param->mass_b_1S*param->mass_b_1S-mV[ie]*mV[ie],2.)+mV[ie]*mV[ie]*GammaTot[ie]*GammaTot[ie]);
	Rccres*=9.*s*param->mass_b_1S*param->mass_b_1S/param->inv_alpha_em/param->inv_alpha_em;
	
	double Rcccont=0.;
	if((s>=0.6)&&(s<0.69)) Rcccont=-6.80+11.33*s;
	if((s>=0.69)&&(s<=1.)) Rcccont=1.02;
	
	return Rccres+Rcccont;
}

/*----------------------------------------------------------------------*/

double complex g_bsll_parametrized(double z, double s, struct parameters* param)
{    
	double integ1=0.;
	double integ2=0.;
	int ie;
	int nb1=25;
	int nb2=25;
	double epsilon=1.e-2;
	double smax=1.;
		
	double sp=4.*param->m_D*param->m_D/param->mass_b_1S/param->mass_b_1S;
	double dsp=(s-epsilon-sp)/nb1;
	integ1+=Rcchad(sp,param)/(sp-s)/sp/2.;
	for(ie=1;ie<nb1;ie++)
	{
		sp+=dsp;
		integ1+=Rcchad(sp,param)/(sp-s)/sp;
	}
	sp=s-epsilon;
	integ1+=Rcchad(sp,param)/(sp-s)/sp/2.;
	integ1*=dsp;
	
	sp=s+epsilon;
	dsp=(smax-sp)/nb2;
	integ2+=Rcchad(sp,param)/(sp-s)/sp/2.;
	for(ie=1;ie<nb2;ie++)
	{
		sp+=dsp;
		integ2+=Rcchad(sp,param)/(sp-s)/sp;
	}
	sp=smax;
	integ2+=Rcchad(sp,param)/(sp-s)/sp/2.;
	integ2*=dsp;
	
	return g_bsll(z,0.)+s/3.*(integ1+integ2)+I*pi/3.*Rcchad(s,param);
}

/*----------------------------------------------------------------------*/

double glambda_bsll(double z)
{
	return 3.-8.*z+24.*z*z-24.*z*z*z+5.*pow(z,4.)+12.*z*z*log(z);
}

/*----------------------------------------------------------------------*/

double grho_bsll(double z)
{
	return 77.-88.*z+24.*z*z-8.*z*z*z+5.*pow(z,4.)+48.*log(z)+36.*z*z*log(z);
}

/*----------------------------------------------------------------------*/

double f_bsll(double x)
{
	return 1.-8.*x*x+8.*pow(x,6.)-pow(x,8.)-24.*pow(x,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double h_bsll(double z)
{
	return -(1.-z*z)*(25./4.-239./3.*z+25./4.*z*z)
	+z*log(z)*(20.+90.*z-4./3.*z*z+17./3.*z*z*z)
	+z*z*log(z)*log(z)*(36.+z*z)
	+(1.-z*z)*(17./3.-64./3.*z+17./3.*z*z)*log(1.-z)
	-4.*(1.+30.*z*z+pow(z,4.))*log(z)*log(1.-z)
	-(1.+16.*z*z+pow(z,4.))*(6.*Li2(z)-pi*pi)
	-32.*pow(z,1.5)*(1.+z)*(pi*pi-4.*Li2(sqrt(z))+4.*Li2(-sqrt(z))-2.*log(z)*log((1.-sqrt(z))/(1.+sqrt(z))));
}

/*----------------------------------------------------------------------*/

double kappa_bsll(double x, double alfas)
{
	double z=x*x;
	return 1.-2./3./pi*alfas*h_bsll(z)/f_bsll(x);
}

/*----------------------------------------------------------------------*/

double sigma_bsll(double s)
{
	return -4./3.*Li2(s)-2./3.*log(s)*log(1.-s)-2./9.*pi*pi-log(1.-s)-2./9.*(1.-s)*log(1.-s);
}

/*----------------------------------------------------------------------*/

double sigma7_bsll(double s, double L)
{
	return sigma_bsll(s)+1./6.-8./3.*L;
}

/*----------------------------------------------------------------------*/

double sigma9_bsll(double s)
{
	return sigma_bsll(s)+1.5;
}

/*----------------------------------------------------------------------*/

double f7_bsll(double s)
{
	return 1./6./(s-1.)/(s-1.)*(24.*(1.+13.*s-4.*s*s)*Li2(sqrt(s))+12.*(1.-17.*s+6.*s*s)*Li2(s)+6.*s*(6.-7.*s)*log(s)
	+24.*(1.-s)*(1.-s)*log(s)*log(1.-s)+12.*(-13.+16.*s-3.*s*s)*(log(1.-sqrt(s))-log(1.-s))
	+39.-2.*pi*pi+252.*s-26.*pi*pi*s+21.*s*s+8.*pi*pi*s*s-180.*sqrt(s)-132.*s*sqrt(s));
}

/*----------------------------------------------------------------------*/

double f9_bsll(double s)
{
	return -1./6./(s-1.)/(s-1.)*(48.*s*(-5.+2.*s)*Li2(sqrt(s))+24.*(-1.+7.*s-3.*s*s)*Li2(s)+6.*s*(-6.+7.*s)*log(s)
	-24.*(1.-s)*(1.-s)*log(s)*log(1.-s)+24.*(5.-7.*s+2.*s*s)*(log(1.-sqrt(s))-log(1.-s))
	-21.-156.*s+20.*pi*pi*s+9.*s*s-8.*pi*pi*s*s+120.*sqrt(s)+48.*s*sqrt(s));
}

/*----------------------------------------------------------------------*/

double complex Gm1_bsll(double t)
{
	if(t>4.) return -2.*I*pi*log((sqrt(t)+sqrt(t-4.))/2.)-pi*pi/2.+2.*pow(log((sqrt(t)+sqrt(t-4.))/2.),2.);
	else return 2.*pi*atan(sqrt((4.-t)/t))-pi*pi/2.-2.*pow(atan(sqrt((4.-t)/t)),2.);
}

/*----------------------------------------------------------------------*/

double complex G0_bsll(double t)
{
	if(t>4.) return -I*pi*sqrt((t-4.)/t)-2.+2.*sqrt((t-4.)/t)*log((sqrt(t)+sqrt(t-4.))/2.);

	else return pi*sqrt((4.-t)/t)-2.-2.*sqrt((4.-t)/t)*atan(sqrt((4.-t)/t));
}

/*----------------------------------------------------------------------*/

double complex Di23_bsll(double s, double w, double z)
{
	return -2.+4./(w-s)*(z*Gm1_bsll(s/z)-z*Gm1_bsll(w/z)-s/2.*G0_bsll(s/z)+s/2.*G0_bsll(w/z));
}

/*----------------------------------------------------------------------*/

double complex Di27_bsll(double s, double w, double z)
{
	return 2.*(G0_bsll(s/z)-G0_bsll(w/z));
}

/*----------------------------------------------------------------------*/

double tau77_bsll(double s)
{
	return -2./9./(2.+s)*(2.*(1.-s)*(1.-s)*log(1.-s)+6.*s*(2.-2.*s-s*s)/(1.-s)/(1.-s)*log(s)+(11.-7.*s-10.*s*s)/(1.-s));
}

/*----------------------------------------------------------------------*/

double tau99_bsll(double s)
{
	return -4./9./(1.+2.*s)*(2.*(1.-s)*(1.-s)*log(1.-s)+3.*s*(1.+s)*(1.-2.*s)/(1.-s)/(1.-s)*log(s)+3.*(1.-3.*s*s)/(1.-s));
}

/*----------------------------------------------------------------------*/

double tau79_bsll(double s)
{
	return -4.*(1.-s)*(1.-s)/9./s*log(1.-s)-4.*s*(3.-2.*s)*log(s)/9./(1.-s)/(1.-s)-2./9.*(5.-3.*s)/(1.-s);
}

/*----------------------------------------------------------------------*/

double tau710_bsll(double s)
{
	return -5./2.+1./3./(1.-3.*s)-1./3.*s*(6.-7.*s)*log(s)/(1.-s)/(1.-s)-1./9.*(3.-7.*s+4.*s*s)*log(1.-s)/s+f7_bsll(s)/3.;
}

/*----------------------------------------------------------------------*/

double tau910_bsll(double s)
{
	return -5./2.+1./3./(1.-s)-1./3.*s*(6.-7.*s)*log(s)/(1.-s)/(1.-s)-2./9.*(3.-5.*s+2.*s*s)*log(1.-s)/s+f9_bsll(s)/3.;
}

/*----------------------------------------------------------------------*/

double tau22_bsll(double w, double s, double z)
{
	return 8./27.*(w-s)*(1.-w)*(1.-w)/s/w/w/w*((3.*w*w+2.*s*s*(2.+w)-s*w*(5.-2.*w))*pow(cabs(Di23_bsll(s,w,z)),2.)
	+(2.*s*s*(2.+w)+s*w*(1.+2.*w))*pow(cabs(Di27_bsll(s,w,z)),2.)
	+4.*s*(w*(1.-w)-s*(2.+w))*creal(Di23_bsll(s,w,z)*conj(Di27_bsll(s,w,z))));
}

/*----------------------------------------------------------------------*/

double integ_tau22(double s, double z)
{
	double integ=0.;
	int ie;
	int nb=20;
	double w=s;
	double dw=(1.-s)/nb;
	
	integ=tau22_bsll(s*1.01,s,z)/2.;
	for(ie=1;ie<nb;ie++)
	{
		w+=dw;
		integ+=tau22_bsll(w,s,z);
	}
	integ+=tau22_bsll(0.99,s,z)/2.;
	
	return integ*dw;
}

/*----------------------------------------------------------------------*/

double tau78_bsll(double s)
{
	return 8./9./s*(25.-2.*pi*pi-27.*s+3.*s*s-s*s*s+12.*(s+s*s)*log(s)
	+6.*pow((pi/2.-atan((2.-4.*s+s*s)/(2.-s)*sqrt(s)*sqrt(4.-s))),2.)
	-24.*creal(CLi2((s-I*sqrt(s)*sqrt(4.-s))/2.))
	-12.*((1.-s)*sqrt(s)*sqrt(4.-s)-atan((sqrt(s)*sqrt(4.-s))/(2.-s)))
	*(atan(sqrt((4.-s)/s))-atan((sqrt(s)*sqrt(4.-s))/(2.-s))));
}

/*----------------------------------------------------------------------*/

double tau88_bsll(double s)
{
	return 4./27./s*(-8.*pi*pi+(1.-s)*(77.-s-4.*s*s)-24.*Li2(1.-s)
	+3.*(10.-4.*s-9.*s*s+8.*log((sqrt(s))/(1.-s)))*log(s)
	+48.*creal(CLi2((3.-s)/2.+I*(1.-s)*sqrt(4.-s)/2./sqrt(s)))
	-6.*((20.*s+10.*s*s-3.*s*s*s)/sqrt(s)/sqrt(4.-s)-8.*pi+8.*atan(sqrt((4.-s)/s)))
	*(atan(sqrt((4.-s)/s))-atan((sqrt(s)*sqrt(4.-s))/(2.-s))));
}

/*----------------------------------------------------------------------*/

double tau89_bsll(double s)
{
	return 2./3.*(s*(4.-s)-3.-4.*log(s)*(1.-s-s*s)
	-8.*creal(CLi2(s/2.+I*sqrt(s)*sqrt(4.-s)/2.)-CLi2((-2.+s*(4.-s))/2.+I*((2.-s)*sqrt(s)*sqrt(4.-s))/2.))
	+4.*(s*s*sqrt((4.-s)/s)+2.*atan(sqrt(s)*sqrt(4.-s)/(2.-s)))
	*(atan(sqrt((4.-s)/s))-atan((sqrt(s)*sqrt(4.-s))/(2.-s))));
}

/*----------------------------------------------------------------------*/

double complex tau27_bsll(double w, double s, double z)
{
	return 8./3./s/w*(((1.-w)*(4.*s*s-s*w+w*w)+s*w*(4.+s-w)*log(w))*Di23_bsll(s,w,z)
	-4.*s*s*(1.-w)+s*w*(4.+s-w)*log(w)*Di27_bsll(s,w,z));
}

/*----------------------------------------------------------------------*/

double complex integ_tau27(double s, double z)
{
	double complex integ=0.;
	int ie;
	int nb=20.;
	double w=s;
	double dw=(1.-s)/nb;
	
	integ=tau27_bsll(s*1.01,s,z)/2.;
	for(ie=1;ie<nb;ie++)
	{
		w+=dw;
		integ+=tau27_bsll(w,s,z);
	}
	integ+=tau27_bsll(0.99,s,z)/2.;
	
	return integ*dw;
}

/*----------------------------------------------------------------------*/

double complex tau28_bsll(double w, double s, double z)
{
	return 8./9./s/w/(w-s)*((pow(w-s,2.)*(2.*s-w)*(1.-w))*Di23_bsll(s,w,z)
	-(2.*s*pow(w-s,2.)*(1.-w))*Di27_bsll(s,w,z)
	+s*w*((1.+2.*s-2.*w)*Di23_bsll(s,w,z)-2.*(1.+s-w)*Di27_bsll(s,w,z))*log(s/((1.+s-w)*(w*w+s*(1.-w)))));
}

/*----------------------------------------------------------------------*/

double complex integ_tau28(double s, double z)
{
	double complex integ=0.;
	int ie;
	int nb=20;
	double w=s;
	double dw=(1.-s)/nb;
	
	integ=tau28_bsll(s*1.01,s,z)/2.;
	for(ie=1;ie<nb;ie++)
	{
		w+=dw;
		integ+=tau28_bsll(w,s,z);
	}
	integ+=tau28_bsll(0.99,s,z)/2.;
	
	return integ*dw;
}

/*----------------------------------------------------------------------*/

double complex tau29_bsll(double w, double s, double z)
{
	return 4./3./w*((2.*s*(1.-w)*(s+w)+4.*s*w*log(w))*Di23_bsll(s,w,z)
	-(2.*s*(1.-w)*(s+w)+w*(3.*s+w)*log(w))*Di27_bsll(s,w,z));
}

/*----------------------------------------------------------------------*/

double complex integ_tau29(double s, double z)
{
	double complex integ=0.;
	int ie;
	int nb=20;
	double w=s;
	double dw=(1.-s)/nb;
	
	integ=tau29_bsll(s*1.01,s,z)/2.;
	for(ie=1;ie<nb;ie++)
	{
		w+=dw;
		integ+=tau29_bsll(w,s,z);
	}
	integ+=tau29_bsll(0.99,s,z)/2.;
	
	return integ*dw;
}

/*----------------------------------------------------------------------*/

double complex tau810_bsll(double s)
{
	return 1./6./(1.-s)/(1.-s)*(3.*((1.-sqrt(s))*(1.-sqrt(s))*(23.-6.*sqrt(s)-s)+4.*(1.-s)*(7.+s)*log(1.+sqrt(s))
	+2.*s*(1.+s-log(s))*log(s))+2.*(-3.*pi*pi*(1.+2.*s)+6.*(3.-s)*s*log(2.-sqrt(s))
	-36.*(1.+2.*s)*CLi2(-sqrt(s))-6.*sqrt(s/(4.-s))*(2.*(-3.+s)*s*atan((2.+sqrt(s))/sqrt(4.-s))
	+2.*pi*log(2.-sqrt(s))-atan(sqrt((4.-s)/s))*((-3.+s)*s+4.*log(2.-sqrt(s)))
	-atan(sqrt(s*(4.-s))/(2.-s))*((-3.+s)*s-log(s))+4.*creal(I*CLi2((-2.+I*sqrt(4.-s)+sqrt(s))*sqrt(s)/(I*sqrt(4.-s)-sqrt(s))))
	-2.*creal(I*CLi2(I/2.*sqrt(4.-s)*(1.-s)*sqrt(s)+(3.-s)*s/2.)))));
}

/*----------------------------------------------------------------------*/

double complex dtau210_bsll(double w, double s, double z)
{
	return -s/(s-w)/(1.-s)/(1.-s)*((4.*(1.-s)*(1.+w)-2.*fabs(s-w*w)*(w*(3.+w)-s*(1.-w))/w/w
	+(2.+5.*w+2.*w*w+s*(3.+4.*w))*log((s+w*w+fabs(s-w*w))/2./w)-(s-w)/s/sqrt((1.+w)*(1.+w)-4.*s)
	*(w*(2.-w)-s*(6.-5.*w))*(log(1.+w-s*(3.-w)+(1.-s)*sqrt((1.+w)*(1.+w)-4.*s))
	-log(s*(1.-3.*w)+w*w*(1.+w)+fabs(s-w*w)*sqrt((1.+w)*(1.+w)-4.*s))))*Di23_bsll(s,w,z)
	-(2.*(1.-s)*(1.+2.*w)-2.*fabs(s-w*w)*(w*(2.+w)-s*(1.-w))/w/w
	+2.*(s*(1.+2.*w)+w*(2.+w))*log((s+w*w+fabs(s-w*w))/2./w)
	+4.*(1.-w)*(s-w)/sqrt((1.+w)*(1.+w)-4.*s)*(log(1.+w-s*(3.-w)+(1.-s)*sqrt((1.+w)*(1.+w)-4.*s))
	-log(s*(1.-3.*w)+w*w*(1.+w)+fabs(s-w*w)*sqrt((1.+w)*(1.+w)-4.*s))))*Di27_bsll(s,w,z));
}

/*----------------------------------------------------------------------*/

double complex tau210_bsll(double s, double z)
{
	double complex integ=0.;
	int ie;
	int nb=20;
	double w=s;
	double dw=(1.-s)/nb;

	integ=dtau210_bsll(s*1.01,s,z)/2.;
	for(ie=1;ie<nb;ie++)
	{
		w+=dw;
		integ+=dtau210_bsll(w,s,z);
	}
	integ+=dtau210_bsll(0.99,s,z)/2.;

	return integ*dw;
}
/*----------------------------------------------------------------------*/

double complex F_bsll(double r)
{
	if(r<1.) return 3./2./r*(1./sqrt(r*(1.-r))*atan(sqrt(r/(1.-r)))-1.);
	else return 3./2./r*(1./2./sqrt(r*(r-1.))*(log((1.-sqrt(1.-1./r))/(1.+sqrt(1.-1./r)))+I*pi)-1.);
}

/*----------------------------------------------------------------------*/

double complex F87_bsll(double s, double L)
{
	return 4.*pi*pi/27.*(2.+s)/pow(1.-s,4.)-4./9.*(11.-16.*s+8.*s*s)/(1.-s)/(1.-s)
	-8./9.*sqrt(s)*sqrt(4.-s)/pow(1.-s,3.)*(9.-5.*s+2.*s*s)*asin(sqrt(s)/2.)
	-16./3.*(2.+s)/pow(1.-s,4.)*pow(asin(sqrt(s)/2.),2.)
	-8.*s/9./(1.-s)*log(s)-32./9.*L-I*8./9.*pi;
}

/*----------------------------------------------------------------------*/

double F89_bsll(double s)
{
	return -8.*pi*pi/27.*(4.-s)/pow(1.-s,4.)+8./9.*(5.-2.*s)/(1.-s)/(1.-s)
	+16./9.*sqrt(4.-s)/sqrt(s)/pow(1.-s,3.)*(4.+3.*s-s*s)*asin(sqrt(s)/2.)
	+32./3.*(4.-s)/pow(1.-s,4.)*pow(asin(sqrt(s)/2.),2.)+16./9./(1.-s)*log(s);
}

/*----------------------------------------------------------------------*/

double w77em(double s, double L)
{
	return L*(s/2./(1.-s)/(2.+s)+log(1.-s)-s*(-3.+2.*s*s)/2./(1.-s)/(1.-s)/(2.+s)*log(s));
}

/*----------------------------------------------------------------------*/

double w79em(double s, double L)
{
	return L*(-0.5/(1.-s)+log(1.-s)+(-1.+2.*s-2.*s*s)/2./(1.-s)/(1.-s)*log(s));
}

/*----------------------------------------------------------------------*/

double w710em(double s, double L)
{
	return L*((7.-16.*sqrt(s)+9.*s)/4./(1.-s)+log(1.-sqrt(s))+(1.+3.*s)/(1.-s)*log((1.+sqrt(s))/2.)-s*log(s)/(1.-s));
}

/*----------------------------------------------------------------------*/

double w99em(double s, double L)
{
	return L*(-(1.+4.*s-8.*s*s)/6./(1.-s)/(1.+2.*s)+log(1.-s)-(1.-6.*s*s+4.*s*s*s)*log(s)/2./(1.-s)/(1.-s)/(1.+2.*s))
	-Li2(s)/9.+4./27.*pi*pi-(37.-3.*s-6.*s*s)/72./(1.-s)/(1.+2.*s)-((41.+76.*s)*log(1.-s))/36./(1.+2.*s)
	+((6.-10.*s-17.*s*s+14.*s*s*s)/18./(1.-s)/(1.-s)/(1.+2.*s)+17.*log(1.-s)/18.)*log(s)-(1.-6.*s*s+4.*s*s*s)*log(s)*log(s)/2./(1.-s)/(1.-s)/(1.+2.*s);
}

/*----------------------------------------------------------------------*/

double w910em(double s, double L)
{
	return L*(-(5.-16.*sqrt(s)+11.*s)/4./(1.-s)+log(1.-sqrt(s))+(1.-5.*s)/(1.-s)*log((1.+sqrt(s))/2.)-(1.-3.*s)*log(s)/(1.-s));

}

/*----------------------------------------------------------------------*/

double w1010em(double s, double L)
{
	return L*(-(1.+4.*s-8.*s*s)/6./(1.-s)/(1.+2.*s)+log(1.-s)-(1.-6.*s*s+4.*s*s*s)*log(s)/2./(1.-s)/(1.-s)/(1.+2.*s));

}

/*----------------------------------------------------------------------*/

double w22em(double s, double L, double mub)
{
	double Sigma1,Sigma2;
	
	if(s<=0.4)
	{
		Sigma1=23.787-120.948*s+365.373*s*s-584.206*s*s*s;
		Sigma2=11.488-36.987*s+255.330*s*s-812.388*s*s*s+1011.791*s*s*s*s;
	}
	else
	{
		double d=(1.-s);
		Sigma1=-148.061*d*d+492.539*d*d*d-1163.847*pow(d,4.)+1189.528*pow(d,5.);
		Sigma2=-221.904*d*d+900.822*d*d*d-2031.620*pow(d,4.)+1984.303*pow(d,5.);
	}
	
	return L*(Sigma2/8./(1.-s)/(1.-s)/(1.+2.*s)+Sigma1/9./(1.-s)/(1.-s)/(1.+2.*s)*log(mub/5.))
	+64./81.*w1010em(s,L)*log(mub/5.)*log(mub/5.);
}

/*----------------------------------------------------------------------*/

double complex w27em(double s, double L, double mub)
{
	double Sigma3,Sigma3I;
	
	if(s<=0.4)
	{
		Sigma3=109.311-846.039*s+2890.115*s*s-4179.072*s*s*s;
		Sigma3I=4.606+17.650*s-53.244*s*s+348.069*s*s*s;
	}
	else
	{
		double d=(1.-s);
		Sigma3=-298.730*d*d+828.0675*d*d*d-2217.6355*pow(d,4.)+2241.792*pow(d,5.);
		Sigma3I=-528.759*d*d+2095.723*d*d*d-4681.843*pow(d,4.)+5036.677*pow(d,5.);
	}
	
	return L*((Sigma3+I*Sigma3I)/96./(1.-s)/(1.-s))+8./9.*w79em(s,L)*log(mub/5.);
}

/*----------------------------------------------------------------------*/

double complex w29em(double s, double L, double mub)
{
	double Sigma1,Sigma1I;
	
	if(s<=0.4)
	{
		Sigma1=23.787-120.948*s+365.373*s*s-584.206*s*s*s;
		Sigma1I=1.653+6.009*s-17.080*s*s+115.880*s*s*s;
	}
	else
	{
		double d=(1.-s);
		Sigma1=-148.061*d*d+492.539*d*d*d-1163.847*pow(d,4.)+1189.528*pow(d,5.);
		Sigma1I=-261.287*d*d+1170.856*d*d*d-2546.948*pow(d,4.)+2540.023*pow(d,5.);
	}
	
	return L*((Sigma1+I*Sigma1I)/8./(1.-s)/(1.-s)/(1.+2.*s))+16./9.*w1010em(s,L)*log(mub/5.);
}

/*----------------------------------------------------------------------*/

double complex w210em(double s, double L, double a, double mub)
{
	double Sigma7,Sigma7I;
	
	if(s<=0.4)
	{
		Sigma7=-0.259023-28.424*s+205.533*s*s-603.219*s*s*s+722.031*s*s*s*s;
		Sigma7I=(-12.20658-215.8208*(s-a)+412.1207*(s-a)*(s-a))*(s-a)*(s-a)*(s>a);
	}
	else
	{
		double d=(1.-s);
		Sigma7=77.0256*d*d-264.705*d*d*d+595.814*pow(d,4.)-610.1637*pow(d,5.);
		Sigma7I=135.858*d*d-618.990*d*d*d+1325.040*pow(d,4.)-1277.170*pow(d,5.);
	}


	return L*(-(Sigma7+I*Sigma7I)/24./s/(1.-s)/(1.-s))+8./9.*w910em(s,L)*log(mub/5.);
}

/*----------------------------------------------------------------------*/

double dBR_BXsll_dshat(int gen, double shat, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	double ml;
	if(gen==1) ml=param->mass_e;
	else if(gen==3) ml=param->mass_tau;
	else ml=param->mass_mu;
			
	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	
	int ie;
	
	double complex Cmub[11];
	
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*(C1b[ie]+alphas_mub/4./pi*C2b[ie]);
	
	double mchat=param->mass_c/param->mass_b_1S;
	
	double z=mchat*mchat;
 
 	double complex C7eff=Cmub[7];
			
	double complex g_bsll_mchat=g_bsll_parametrized(mchat,shat,param);
	double complex g_bsll_1=g_bsll(1.,shat);
	double complex g_bsll_0=g_bsll(0.,shat);
	
	double complex C9eff=Cmub[9]
	+(-32./27.*Cmub[1]-8./9.*Cmub[2]-16./9.*Cmub[3]+32./27.*Cmub[4]-112./9.*Cmub[5]+512./27.*Cmub[6])*log(param->mass_b_1S/mu_b)
	+4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
	+g_bsll_mchat*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
	+g_bsll_1*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
	+g_bsll_0*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);
	
	double complex C90eff=C0b[9] +(-32./27.*C0b[1]-8./9.*C0b[2]-16./9.*C0b[3]+32./27.*C0b[4]-112./9.*C0b[5]+512./27.*C0b[6])*log(param->mass_b_1S/mu_b)
	+4./3.*C0b[3]+64./9.*C0b[5]+64./27.*C0b[6]
	+g_bsll_mchat*(4./3.*C0b[1]+C0b[2]+6.*C0b[3]+60.*C0b[5])
	+g_bsll_1*(-7./2.*C0b[3]-2./3.*C0b[4]-38.*C0b[5]-32./3.*C0b[6])
	+g_bsll_0*(-1./2.*C0b[3]-2./3.*C0b[4]-8.*C0b[5]-32./3.*C0b[6]);
	
	double complex C10eff=Cmub[10];

	double complex sigma7=sigma7_bsll(shat,log(mu_b/param->mass_b_1S));
	double complex F17=F17_bsll(shat,z,log(mu_b/param->mass_b_1S));
	double complex F27=F27_bsll(shat,z,log(mu_b/param->mass_b_1S));
	double complex F87=F87_bsll(shat,log(mu_b/param->mass_b_1S));
	double complex sigma9=sigma9_bsll(shat);
	double complex F19=F19_bsll(shat,z,log(mu_b/param->mass_b_1S));
	double complex F29=F29_bsll(shat,z,log(mu_b/param->mass_b_1S));
	double complex F89=F89_bsll(shat);
	
 	double complex C7new=(1.+alphas_mub/pi*sigma7)*C7eff
	-alphas_mub/4./pi*(C0b[1]*F17+C0b[2]*F27+C0b[8]*F87);
	
	double complex C9new=(1.+alphas_mub/pi*sigma9)*C9eff
	-alphas_mub/4./pi*(C0b[1]*F19+C0b[2]*F29+C0b[8]*F89);

	double complex C10new=(1.+alphas_mub/pi*sigma9)*C10eff;
 
	double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
	double complex CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
	
	double t=ml/param->mass_b_1S;
	
	double lambda2=param->lambda2;
	
	double f=f_bsll(param->mass_c/param->mass_b_1S);
	double kappa=kappa_bsll(param->mass_c/param->mass_b_1S,alphas_mub);
	double Fbsll=F_bsll(shat/4./z);

	double rho1=param->rho1;

	double dR_mb2=0.;
	double dR_mb3=0.;
	double dR_mc2=0.;
	
	if(gen!=3)
	{
		dR_mb2=3.*lambda2/2./param->mass_b_1S/param->mass_b_1S*(1./param->inv_alpha_em/param->inv_alpha_em/4./pi/pi*pow(cabs(param->Vts/param->Vcb),2.)/f/kappa*(-(6.+3.*shat-5.*shat*shat*shat)*4.*C7new*conj(C7new)/shat
		+(1.-15.*shat*shat+10.*shat*shat*shat)*(C9new*conj(C9new)+C10new*conj(C10new))
		-4.*(5.+6.*shat-7.*shat*shat)*creal(C7new*conj(C9new))));
		
		if(shat<0.4) dR_mb3=-rho1/pow(param->mass_b_1S,3.)*(1./param->inv_alpha_em/param->inv_alpha_em/4./pi/pi*pow(cabs(param->Vts/param->Vcb),2.)/f/kappa*(
		(5.*pow(shat,4.)+19.*pow(shat,3.)+9.*shat*shat-7.*shat+22.)/6./(1.-shat)*4.*C7new*conj(C7new)/shat
		+(10.*pow(shat,4.)+23.*pow(shat,3.)-9.*shat*shat+13.*shat+11.)/6./(1.-shat)*(C9new*conj(C9new)+C10new*conj(C10new))
		+4.*(-3.*pow(shat,3.)+17.*shat*shat-shat+3.)/2./(1.-shat)*creal(C7new*conj(C9new))));
	
		dR_mc2=8.*lambda2/9./param->mass_c/param->mass_c/param->inv_alpha_em/param->inv_alpha_em/4./pi/pi*cabs(param->Vts*conj(param->Vcs)/param->Vtb/conj(param->Vcb))/f/kappa*(1.-shat)*(1.-shat)
		*creal((1.+6.*shat-shat*shat)/shat*Fbsll*Cmub[2]*conj(C7new)+(2.+shat)*Fbsll*Cmub[2]*conj(C9new));
	}
		
	double complex c78=4./3.*C0b[7]*conj(C0b[8]);
	double complex c89=4./3.*C0b[8]*conj(C90eff);
	double complex c88=4./3.*C0b[8]*conj(C0b[8]);
	
	double tau78=tau78_bsll(shat);
	double tau89=tau89_bsll(shat);
	double tau88=tau88_bsll(shat);
	double integtau22=integ_tau22(shat,z);
	double complex integtau27=integ_tau27(shat,z);
	double complex integtau28=integ_tau28(shat,z);
	double complex integtau29=integ_tau29(shat,z);
	
	double dR_bremA=4.*pow(1./(4.*pi*param->inv_alpha_em),2.)*alphas_mub/4./pi/f/kappa*pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*(2.*creal(c78*tau78+c89*tau89)+c88*tau88);	

	double complex c11=C0b[1]*conj(C0b[1])/27.;
	double complex c12=-4./9.*C0b[1]*conj(C0b[2]);
	double complex c22=4./3.*C0b[2]*conj(C0b[2]);
	double complex c17=-2./9.*C0b[1]*conj(C0b[7]);
	double complex c27=4./3.*C0b[2]*conj(C0b[7]);
	double complex c18=-2./9.*C0b[1]*conj(C0b[8]);
	double complex c28=4./3.*C0b[2]*conj(C0b[8]);
	double complex c19=-2./9.*C0b[1]*conj(C90eff);
	double complex c29=4./3.*C0b[2]*conj(C90eff);
	
	double dR_bremB=4.*pow(1./(4.*pi*param->inv_alpha_em),2.)*alphas_mub/4./pi/f/kappa*pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*
	(creal(c11+c12+c22)*integtau22+2.*creal((c17+c27)*integtau27+creal(c18+c28)*integtau28+creal(c19+c29)*integtau29));

	double L=2.*log(param->mass_b_1S/ml);

	double dR_em=pow(1./param->inv_alpha_em/4./pi,3.)*4.*pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*(1.-shat)*(1.-shat)/f/kappa*
	(8.*(1.+2.*shat)*(Cmub[9]*conj(Cmub[9])*w99em(shat,L)+Cmub[10]*conj(Cmub[10])*w1010em(shat,L)
	+creal((Cmub[2]+4./3.*Cmub[1])*conj(Cmub[9])*w29em(shat,L,mu_b))
	+pow(cabs(Cmub[2]+4./3.*Cmub[1]),2.)*w22em(shat,L,mu_b))
	+96.*(creal(Cmub[7]*conj(Cmub[9]))*w79em(shat,L)
	+creal((Cmub[2]+4./3.*Cmub[1])*conj(Cmub[7])*w27em(shat,L,mu_b)))
	+8.*(4.+8./shat)*Cmub[7]*conj(Cmub[7])*w77em(shat,L));


	double complex C7effp=Cpb[7];
		
	double complex C9effp=Cpb[9]
	+(-32./27.*Cpb[1]-8./9.*Cpb[2]-16./9.*Cpb[3]+32./27.*Cpb[4]-112./9.*Cpb[5]+512./27.*Cpb[6])*log(param->mass_b_1S/mu_b)
	+4./3.*Cpb[3]+64./9.*Cpb[5]+64./27.*Cpb[6]
	+g_bsll_mchat*(4./3.*Cpb[1]+Cpb[2]+6.*Cpb[3]+60.*Cpb[5])
	+g_bsll_1*(-7./2.*Cpb[3]-2./3.*Cpb[4]-38.*Cpb[5]-32./3.*Cpb[6])
	+g_bsll_0*(-1./2.*Cpb[3]-2./3.*Cpb[4]-8.*Cpb[5]-32./3.*Cpb[6]);
	
	double complex C10effp=Cpb[10];

	double complex C7newp=(1.+alphas_mub/pi*sigma7)*C7effp
	-alphas_mub/4./pi*(Cpb[1]*F17+Cpb[2]*F27+Cpb[8]*F87);
	
	double complex C9newp=(1.+alphas_mub/pi*sigma9)*C9effp
	-alphas_mub/4./pi*(Cpb[1]*F19+Cpb[2]*F29+Cpb[8]*F89);

	double complex C10newp=(1.+alphas_mub/pi*sigma9)*C10effp;
 
	double complex CQ1p=CQpb[1];
	double complex CQ2p=CQpb[2];
	
	double dR_mb2p=0.;
	double dR_mb3p=0.;
	double dR_mc2p=0.;
	
	if(gen!=3)
	{
		dR_mb2p=3.*lambda2/2./param->mass_b_1S/param->mass_b_1S*(1./param->inv_alpha_em/param->inv_alpha_em/4./pi/pi*pow(cabs(param->Vts/param->Vcb),2.)/f/kappa*(-(6.+3.*shat-5.*shat*shat*shat)*4.*C7newp*conj(C7newp)/shat
		+(1.-15.*shat*shat+10.*shat*shat*shat)*(C9newp*conj(C9newp)+C10newp*conj(C10newp))
		-4.*(5.+6.*shat-7.*shat*shat)*creal(C7newp*conj(C9newp))));
	
		if(shat<0.4) dR_mb3p=-rho1/pow(param->mass_b_1S,3.)*(1./param->inv_alpha_em/param->inv_alpha_em/4./pi/pi*pow(cabs(param->Vts/param->Vcb),2.)/f/kappa*(
		(5.*pow(shat,4.)+19.*pow(shat,3.)+9.*shat*shat-7.*shat+22.)/6./(1.-shat)*4.*C7newp*conj(C7newp)/shat
		+(10.*pow(shat,4.)+23.*pow(shat,3.)-9.*shat*shat+13.*shat+11.)/6./(1.-shat)*(C9newp*conj(C9newp)+C10newp*conj(C10newp))
		+4.*(-3.*pow(shat,3.)+17.*shat*shat-shat+3.)/2./(1.-shat)*creal(C7newp*conj(C9newp))));
	
		dR_mc2p=8.*lambda2/9./param->mass_c/param->mass_c/param->inv_alpha_em/param->inv_alpha_em/4./pi/pi*cabs(param->Vts*conj(param->Vcs)/param->Vtb/conj(param->Vcb))/f/kappa*(1.-shat)*(1.-shat)
		*creal((1.+6.*shat-shat*shat)/shat*Fbsll*Cpb[2]*conj(C7newp)+(2.+shat)*Fbsll*Cpb[2]*conj(C9newp));
	}
	
	double complex c78p=4./3.*Cpb[7]*conj(Cpb[8]);
	double complex c89p=4./3.*Cpb[8]*conj(C9effp);
	double complex c88p=4./3.*Cpb[8]*conj(Cpb[8]);
		
	double dR_bremAp=4.*pow(1./(4.*pi*param->inv_alpha_em),2.)*alphas_mub/4./pi/f/kappa*pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*(2.*creal(c78p*tau78+c89p*tau89)+c88p*tau88);	

	double complex c11p=Cpb[1]*conj(Cpb[1])/27.;
	double complex c12p=-4./9.*Cpb[1]*conj(Cpb[2]);
	double complex c22p=4./3.*Cpb[2]*conj(Cpb[2]);
	double complex c17p=-2./9.*Cpb[1]*conj(Cpb[7]);
	double complex c27p=4./3.*Cpb[2]*conj(Cpb[7]);
	double complex c18p=-2./9.*Cpb[1]*conj(Cpb[8]);
	double complex c28p=4./3.*Cpb[2]*conj(Cpb[8]);
	double complex c19p=-2./9.*Cpb[1]*conj(C9effp);
	double complex c29p=4./3.*Cpb[2]*conj(C9effp);
	
	double dR_bremBp=4.*pow(1./(4.*pi*param->inv_alpha_em),2.)*alphas_mub/4./pi/f/kappa*pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*(creal(c11p+c12p+c22p)*integtau22+2.*creal((c17p+c27p)*integtau27+creal(c18p+c28p)*integtau28+creal(c19p+c29p)*integtau29));

	double dR_emp=pow(1./param->inv_alpha_em/4./pi,3.)*4.*pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*(1.-shat)*(1.-shat)/f/kappa*
	(8.*(1.+2.*shat)*(Cpb[9]*conj(Cpb[9])*w99em(shat,L)+Cpb[10]*conj(Cpb[10])*w1010em(shat,L)
	+creal((Cpb[2]+4./3.*Cpb[1])*conj(Cpb[9])*w29em(shat,L,mu_b))
	+pow(cabs(Cpb[2]+4./3.*Cpb[1]),2.)*w22em(shat,L,mu_b))
	+96.*(creal(Cpb[7]*conj(Cpb[9]))*w79em(shat,L)
	+creal((Cpb[2]+4./3.*Cpb[1])*conj(Cpb[7])*w27em(shat,L,mu_b)))
	+8.*(4.+8./shat)*Cpb[7]*conj(Cpb[7])*w77em(shat,L));

	return param->BR_BXclnu_exp*(1./param->inv_alpha_em/param->inv_alpha_em/4./pi/pi/f/kappa*(1.-shat)*(1.-shat)*sqrt(1.-4.*t*t/shat)*
	pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*
	(
	pow(cabs(C9new),2.)*(1.+2.*t*t/shat)*(1.+2.*shat)*(1.+alphas_mub/pi*tau99_bsll(shat))
	+4.*pow(cabs(C7new),2.)*(1.+2.*t*t/shat)*(1.+2./shat)*(1.+alphas_mub/pi*tau77_bsll(shat))
	+pow(cabs(C10new),2.)*((1.+2.*shat)+2.*t*t/shat*(1.-4.*shat))*(1.+alphas_mub/pi*tau99_bsll(shat))
	+12.*creal(C7new*conj(C9new))*(1.+2.*t*t/shat)*(1.+alphas_mub/pi*tau79_bsll(shat))
	+1.5*pow(cabs(CQ1),2.)*(shat-4.*t*t)+1.5*pow(cabs(CQ2),2.)*shat+6.*creal(C10new*conj(CQ2))*t
	
	+pow(cabs(C9newp),2.)*(1.+2.*t*t/shat)*(1.+2.*shat)*(1.+alphas_mub/pi*tau99_bsll(shat))
	+4.*pow(cabs(C7newp),2.)*(1.+2.*t*t/shat)*(1.+2./shat)*(1.+alphas_mub/pi*tau77_bsll(shat))
	+pow(cabs(C10newp),2.)*((1.+2.*shat)+2.*t*t/shat*(1.-4.*shat))*(1.+alphas_mub/pi*tau99_bsll(shat))
	+12.*creal(C7newp*conj(C9newp))*(1.+2.*t*t/shat)*(1.+alphas_mub/pi*tau79_bsll(shat))
	+1.5*pow(cabs(CQ1p),2.)*(shat-4.*t*t)+1.5*pow(cabs(CQ2p),2.)*shat+6.*creal(C10newp*conj(CQ2p))*t
	))
	*(1.+3.*lambda2/2./param->mass_b_1S/param->mass_b_1S*glambda_bsll(z)/f	
	-rho1/pow(param->mass_b_1S,3.)/6.*grho_bsll(z)/f)
	
	+(dR_mb2+dR_mb3+dR_mc2+dR_bremA+dR_bremB+dR_em)*param->BR_BXclnu_exp
	+(dR_mb2p+dR_mb3p+dR_mc2p+dR_bremAp+dR_bremBp+dR_emp)*param->BR_BXclnu_exp;
}

/*----------------------------------------------------------------------*/

double dBR_BXsll_dshat_calculator(int gen, double shat, char name[])
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole/2.;
				
	CW_calculator(gen,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(gen,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(gen,Cpb,CQpb,mu_W,mu_b,&param);

	return dBR_BXsll_dshat(gen,shat,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBXsll(int gen, double smin, double smax, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
/* "container" function scanning the SLHA file "name" and calculating BR(B->Xs mu+ mu-) */
{
	int ie;
	int nmax=10;
	if((smin<0.099/param->mass_b_1S/param->mass_b_1S)||(smax-smin>10./param->mass_b_1S/param->mass_b_1S)) nmax=100;
	double BR=0.;
	double s;
			
	double h=(smax-smin)/nmax;	
	s=smin;
	BR=dBR_BXsll_dshat(gen,s,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(ie=1;ie<nmax;ie++)
	{
		s+=h;
		BR+=4.*dBR_BXsll_dshat(gen,s-h/2.,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		BR+=2.*dBR_BXsll_dshat(gen,s,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	}
	s=smax;
	BR+=4.*dBR_BXsll_dshat(gen,s-h/2.,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	BR+=dBR_BXsll_dshat(gen,s,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	
	BR*=h/6.;

	return BR;
}

/*----------------------------------------------------------------------*/

double BRBXsll_lowq2(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
/* "container" function scanning the SLHA file "name" and calculating BR(B->Xs mu+ mu-) */
{
	return BRBXsll(gen,1./param->mass_b_1S/param->mass_b_1S,6./param->mass_b_1S/param->mass_b_1S,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBXsll_highq2(int gen,double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[] , struct parameters* param, double mu_b)
/* "container" function scanning the SLHA file "name" and calculating BR(B->Xs mu+ mu-) */
{
	return BRBXsll(gen,14.2/param->mass_b_1S/param->mass_b_1S,0.999,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBXsmumu_lowq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B->Xs mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole/2.;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBXsll_lowq2(2,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBXstautau_highq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B->Xs tau+ tau-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole/2.;
				
	CW_calculator(3,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(3,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(3,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBXsll_highq2(3,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBXsmumu_highq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B->Xs mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole/2.;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return BRBXsll_highq2(2,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double A_BXsll(int gen,double shat, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	double ml;
	if(gen==1) ml=param->mass_e;
	else if(gen==3) ml=param->mass_tau;
	else ml=param->mass_mu;
	
	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	
	double mchat=param->mass_c/param->mass_b_1S;
	double z=mchat*mchat;

	int ie;
	
	double complex Cmub[11];
	
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*(C1b[ie]+alphas_mub/4./pi*C2b[ie]);
	
 	double complex C7eff=Cmub[7];
		
	double complex g_bsll_mchat=g_bsll_parametrized(mchat,shat,param);
	double complex g_bsll_1=g_bsll(1.,shat);
	double complex g_bsll_0=g_bsll(0.,shat);
	
	double complex C9eff=Cmub[9] +(-32./27.*Cmub[1]-8./9.*Cmub[2]-16./9.*Cmub[3]+32./27.*Cmub[4]-112./9.*Cmub[5]+512./27.*Cmub[6])*log(param->mass_b_1S/mu_b)
	+4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
	+g_bsll_mchat*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
	+g_bsll_1*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
	+g_bsll_0*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);
	
	double complex C10eff=Cmub[10];

	double complex sigma7=sigma7_bsll(shat,log(mu_b/param->mass_b_1S));
	double complex F17=F17_bsll(shat,z,log(mu_b/param->mass_b_1S));
	double complex F27=F27_bsll(shat,z,log(mu_b/param->mass_b_1S));
	double complex F87=F87_bsll(shat,log(mu_b/param->mass_b_1S));
	double complex sigma9=sigma9_bsll(shat);
	double complex F19=F19_bsll(shat,z,log(mu_b/param->mass_b_1S));
	double complex F29=F29_bsll(shat,z,log(mu_b/param->mass_b_1S));
	double complex F89=F89_bsll(shat);

  	double complex C7new=(1.+alphas_mub/pi*sigma7)*C7eff-alphas_mub/4./pi*(C0b[1]*F17+C0b[2]*F27+C0b[8]*F87);
		
	double complex C9new=(1.+alphas_mub/pi*sigma9)*C9eff-alphas_mub/4./pi*(C0b[1]*F19+C0b[2]*F29+C0b[8]*F89);

	double complex C10new=(1.+alphas_mub/pi*sigma9)*C10eff;

	double complex CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
	
	double t=ml/param->mass_b_1S;
		
	double lambda1=-0.15;
	double lambda2=param->lambda2;
	
	double f=f_bsll(param->mass_c/param->mass_b_1S);
	double kappa=kappa_bsll(param->mass_c/param->mass_b_1S,alphas_mub);

	double dA_mb2=param->BR_BXclnu_exp*3.*lambda2/2./param->mass_b_1S/param->mass_b_1S*(1./param->inv_alpha_em/param->inv_alpha_em/4./pi/pi*pow(cabs(param->Vts/param->Vcb),2.)/f/kappa*(shat*creal(C9new*conj(C10new))*(9.+14.*shat-15.*shat*shat)
	+2.*creal(conj(C10new)*C7new)*(7.+10.*shat-9.*shat*shat)));

	double Fbsll=F_bsll(shat/4./z);
	
	double dA_mc2=-param->BR_BXclnu_exp*lambda2/3./param->mass_c/param->mass_c/param->inv_alpha_em/param->inv_alpha_em/4./pi/pi*cabs(param->Vts*conj(param->Vcs)/param->Vtb/conj(param->Vcb))/f/kappa*(1.-shat)*(1.-shat)*creal((1.+3.*shat)*Fbsll*Cmub[2]*conj(C10new));

	double L=2.*log(param->mass_b_1S/ml);

	double complex tau810=tau810_bsll(shat);
	double complex tau210=tau210_bsll(shat,z);
	double complex w710=w710em(shat,L);
	double complex w910=w910em(shat,L);
	double complex w210=w210em(shat,L,pow(4.*param->mass_c*param->mass_c/param->mass_b_1S/param->mass_b_1S,2.),mu_b);

	double dA_brems=param->BR_BXclnu_exp/param->inv_alpha_em/param->inv_alpha_em/2./pi/pi*alphas_mub/3./pi/f/kappa*pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*(1.-shat)*(1.-shat)*(creal(C0b[8]*conj(C10eff))*tau810+creal((C0b[2]-C0b[1]/6.)*conj(C10eff))*tau210);

	double dA_em=pow(1./param->inv_alpha_em/4./pi,3.)*4.*param->BR_BXclnu_exp*pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*(1.-shat)*(1.-shat)/f/kappa*
	(-48.*creal(Cmub[7]*conj(Cmub[10]))*w710
	-24.*shat*(creal(Cmub[9]*conj(Cmub[10]))*w910
	+creal((Cmub[2]+4./3.*Cmub[1])*conj(Cmub[10])*w210)));

	double complex C7effp=Cpb[7];
			
	double complex C9effp=Cpb[9] +(-32./27.*Cpb[1]-8./9.*Cpb[2]-16./9.*Cpb[3]+32./27.*Cpb[4]-112./9.*Cpb[5]+512./27.*Cpb[6])*log(param->mass_b_1S/mu_b)
	+4./3.*Cpb[3]+64./9.*Cpb[5]+64./27.*Cpb[6]
	+g_bsll_mchat*(4./3.*Cpb[1]+Cpb[2]+6.*Cpb[3]+60.*Cpb[5])
	+g_bsll_1*(-7./2.*Cpb[3]-2./3.*Cpb[4]-38.*Cpb[5]-32./3.*Cpb[6])
	+g_bsll_0*(-1./2.*Cpb[3]-2./3.*Cpb[4]-8.*Cpb[5]-32./3.*Cpb[6]);
	
	double C10effp=Cpb[10];

  	double complex C7newp=(1.+alphas_mub/pi*sigma7)*C7effp-alphas_mub/4./pi*(Cpb[1]*F17+Cpb[2]*F27+Cpb[8]*F87);
		
	double complex C9newp=(1.+alphas_mub/pi*sigma9)*C9effp-alphas_mub/4./pi*(Cpb[1]*F19+Cpb[2]*F29+Cpb[8]*F89);

	double complex C10newp=(1.+alphas_mub/pi*sigma9)*C10effp;

	double complex CQ1p=CQpb[1];
		
	double dA_mb2p=param->BR_BXclnu_exp*3.*lambda2/2./param->mass_b_1S/param->mass_b_1S*(1./param->inv_alpha_em/param->inv_alpha_em/4./pi/pi*pow(cabs(param->Vts/param->Vcb),2.)/f/kappa*(shat*creal(C9newp*conj(C10newp))*(9.+14.*shat-15.*shat*shat)
	+2.*creal(conj(C10newp)*C7newp)*(7.+10.*shat-9.*shat*shat)));
	
	double dA_mc2p=-param->BR_BXclnu_exp*lambda2/3./param->mass_c/param->mass_c/param->inv_alpha_em/param->inv_alpha_em/4./pi/pi*cabs(param->Vts*conj(param->Vcs)/param->Vtb/conj(param->Vcb))/f/kappa*(1.-shat)*(1.-shat)*creal((1.+3.*shat)*Fbsll*Cpb[2]*conj(C10newp));

	double dA_bremsp=param->BR_BXclnu_exp/param->inv_alpha_em/param->inv_alpha_em/2./pi/pi*alphas_mub/3./pi/f/kappa*pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*(1.-shat)*(1.-shat)*(creal(Cpb[8]*conj(C10effp))*tau810+creal((Cpb[2]-Cpb[1]/6.)*conj(C10effp))*tau210);

	double dA_emp=pow(1./param->inv_alpha_em/4./pi,3.)*4.*param->BR_BXclnu_exp*pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*(1.-shat)*(1.-shat)/f/kappa*
	(-48.*creal(Cpb[7]*conj(Cpb[10]))*w710
	-24.*shat*(creal(Cpb[9]*conj(Cpb[10]))*w910
	+creal((Cpb[2]+4./3.*Cpb[1])*conj(Cpb[10])*w210)));

	return (-3.*param->BR_BXclnu_exp/param->inv_alpha_em/param->inv_alpha_em/4./pi/pi/f/kappa*(1.-shat)*(1.-shat)*(1.-4.*t*t/shat)*
	pow(cabs(param->Vtb*conj(param->Vts)/param->Vcb),2.)*
	(creal(C9new*conj(C10new))*shat*(1.+alphas_mub/pi*tau910_bsll(shat)) 
	+2.*creal(C7new*conj(C10new))*(1.+alphas_mub/pi*tau710_bsll(shat)) 
 	+creal(C9new*conj(CQ1))*t+2.*creal(C7new*conj(CQ1))*t
 	
 	+creal(C9newp*conj(C10newp))*shat*(1.+alphas_mub/pi*tau910_bsll(shat)) 
	+2.*creal(C7newp*conj(C10newp))*(1.+alphas_mub/pi*tau710_bsll(shat)) 
 	+creal(C9newp*conj(CQ1p))*t+2.*creal(C7newp*conj(CQ1p))*t))
	*(1.+3.*lambda2/2./param->mass_b_1S/param->mass_b_1S*glambda_bsll(z)/f
	+4.*lambda1/3./param->mass_b_1S/param->mass_b_1S*shat/(1.-shat)/(1.-shat))
	+dA_mb2+dA_mc2+dA_brems+dA_em
	+dA_mb2p+dA_mc2p+dA_bremsp+dA_emp;
}

/*----------------------------------------------------------------------*/

double A_BXsll_bin(int gen, double smin, double smax, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
/* "container" function scanning the SLHA file "name" and calculating integrated AFB(B->Xs mu+ mu-) */
{
	int ie;
	int nmax=10;
	if((smin<0.099/param->mass_b_1S/param->mass_b_1S)||(smax-smin>10./param->mass_b_1S/param->mass_b_1S)) nmax=100;
	double AFB=0.;
	double s;
	
	double h=(smax-smin)/nmax;	
	s=smin;
	AFB=A_BXsll(gen,s,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	for(ie=1;ie<nmax;ie++)
	{
		s+=h;
		AFB+=4.*A_BXsll(gen,s-h/2.,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		AFB+=2.*A_BXsll(gen,s,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	}
	s=smax;
	AFB+=4.*A_BXsll(gen,s-h/2.,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	AFB+=A_BXsll(gen,s,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
	
	AFB*=h/6.;

	return AFB;
}

/*----------------------------------------------------------------------*/

double A_BXsll_lowq2(int gen,double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
/* "container" function scanning the SLHA file "name" and calculating integrated AFB(B->Xs mu+ mu-) */
{
	return A_BXsll_bin(gen,1./param->mass_b_1S/param->mass_b_1S,6./param->mass_b_1S/param->mass_b_1S,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double A_BXsll_highq2(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
/* "container" function scanning the SLHA file "name" and calculating integrated AFB(B->Xs mu+ mu-) */
{
	return A_BXsll_bin(gen,14.2/param->mass_b_1S/param->mass_b_1S,0.999,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double A_BXsmumu_lowq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B->Xs mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole/2.;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return A_BXsll_lowq2(2,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double A_BXsmumu_highq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B->Xs mu+ mu-) */
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole/2.;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return A_BXsll_highq2(2,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double A_BXsll_zero(int gen, double complex C0b[], double complex C1b[], double complex C2b[], double complex CQ0b[], double complex CQ1b[], double complex Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	double smin=0.05;
	double smax=0.25;

	double stemp;

	while(fabs(1.-smin/smax)>1.e-4)
	{
		stemp=(smax+smin)/2.;
		if(A_BXsll(gen,stemp,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b)<0.) smin=stemp; else smax=stemp;
	}
	return stemp*param->mass_b_1S*param->mass_b_1S;
}

/*----------------------------------------------------------------------*/

double A_BXsmumu_zero_calculator(char name[])
{
	double complex C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11],CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole/2.;
				
	CW_calculator(2,C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(2,CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(2,Cpb,CQpb,mu_W,mu_b,&param);

	return A_BXsll_zero(2,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}
