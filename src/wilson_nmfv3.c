/*  This C file serves as a dump file for collecting all functions
 */
#ifndef INCLUDE_H
#define INCLUDE_H
#endif
#include "include.h"
#include "isajet.h"
#include "softsusy.h"
#include "suspect.h"
#include "spheno.h"
#include <tgmath.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#pragma GCC diagnostic ignored "-Wunused-result"
// This is temporary, should write initializer code. 
#define N_C 3.
#define EPS_PREC 1E-10
#define MASS_S 0.125
#define inv_alpha_MZ 128
#define SW2 0.23121 // squared weinberg angle at M_Z
#define MU_0 120.0
#define YUKAWA 1 // formula to use to calculate yukawa coupling (default 1)
#define define_CKM 1 // If true, then CKM taken from PDG_2018
#ifndef M_GAMMAl
/** Euler's constant in high precision */
#define M_GAMMAl 0.5772156649015328606065120900824024L
#endif
#ifndef M_LN2l
/** the natural logarithm of 2 in high precision */
#define M_LN2l 0.6931471805599453094172321214581766L
#endif
char** err_msg[500];
void reset_printf () {
      printf("\033[0m");
}

void WARNING(char* msg){
    printf("\033[1;33m"); // yellow
    printf("WARNING : %s\n", msg);
    reset_printf();
}
void ERROR(char* msg){
    printf("\033[1;31m"); // bold red
    printf("ERROR : %s\n", msg);
    reset_printf();
    exit(EXIT_FAILURE);
}

void init_Inputvar( Input_var *var ){
    // init all struct members to i-1
    var->mu = -1  ;    
    var->tan_beta = -1;
    var->m_tr = -1;
    var->mass_gluino = 250;
    var->mass_sneu = -1;
    var->M_sq = 250;
    var->M2 = 50;
}
int sgn( double var ){
    if( var >= 0 ) return 1;
    else return -1;
}
void set_yukawa( Input_var * var, struct parameters *param ){

    if(YUKAWA){

        param->yut[3] = param->mass_t*param->g2 / (sqrt(2.) *(sin(atan(param->tan_beta)))* param->mass_W );
        param->yub[3] = param->mass_b*param->g2 / (sqrt(2.) *(cos(atan(param->tan_beta)))* param->mass_W );
    }
    else{
        param->yut[3] = param->mass_t * sqrt(2.) / (246. * sin(atan(param->tan_beta)));
        param->yub[3] =  param->mass_b * sqrt(2.) / (246. * sin(atan(param->tan_beta)));
        }

}
/* Useful functions */
bool is_int(double x){
    if (fmod(x,1)==0) return true;
    else return false;
}
long double factorial (double x){
    x = (long double) x;
    if (x == 1 || x==0)  return 1;
    else if( x<0 ){
        sprintf(err_msg,"%s called with negative argument.");
        ERROR(err_msg);
    }
    else    return (x * factorial (x - 1));
}
int Krodelta(int i,int j){ /// Kronecker delta function
        if(i == j) return 1 ;
        else return 0;
}
double Beta_E(double i, double j){

    return tgamma(i)*tgamma(j)/tgamma(i+j);
}
long double Pochhammer(double x, int n){
    long double res = 1;
    x = (long double) x;
    if (n == 0) return 1;
    else if( n<0 ){
        sprintf(err_msg, "Pochhammer symbol not defined w.r.t negative integers.");
        ERROR(err_msg);
    }
    else{
        for (int i = 0; i < n; i++)   res *= (x+i) ;
        if (res == 0) printf("Warning in %s, result is null\n", __func__);
        return res; 
    }
}
double digamma(double x){
/// from https://www2.mpia-hd.mpg.de/~mathar/progs/digamma.c and the boost stdlib (https://www.boost.org/doc/libs/1_74_0/libs/math/doc/html/math_toolkit/sf_gamma/digamma.html)
// Tested, seems to work correctly. However, Psi(0) = -inf.
    if(x == 0.) {
        char* issued_warning[100];
        sprintf(issued_warning,"In %s, infinity occured !",__func__);
        WARNING(issued_warning);
    }
    if (x < 0) return digamma(1.-x) - pi/tan(pi*x); // Reflection formula
    else if (x < 1) return digamma(1. + x) - 1./x ;
    else if (x == 1.) return -M_GAMMAl;
    else if (x == 2.) return 1. - M_GAMMAl;
    else if (x == 3.) return 1.5 - M_GAMMAl;
    else if (x > 3.) return 0.5*(digamma(x/2.) + digamma((x+1.)/2.)) +  M_LN2l ;
    else{
        static long double Kncoe[] = { .30459198558715155634315638246624251L,
                .72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
                .27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
                .17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
                .11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
                .83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
                .59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
                .42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
                .304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
                .21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
                .15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
                .11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
                .80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
                .58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
                .41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L } ;

            register long double Tn_1 = 1.0L ;/* T_{n-1}(x),     started at n=1 */
            register long double Tn = x-2.0L ;/* T_{n}(x) , s   tarted at n=1 */
            register long double resul = Kncoe[0] + Kncoe[1]*Tn ;
            x -= 2.0L ;
            for(int n = 2 ; n < sizeof(Kncoe)/sizeof(long double) ;n++){
                const long double Tn1 = 2.0L * x * Tn - Tn_1 ;//Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun
                resul += Kncoe[n]*Tn1 ;
                Tn_1 = Tn ;
                Tn = Tn1 ;                                    
            }
            return resul ;
    }
}

double  hypgeo2F1_transform15_3_13(double a, double b, double c, double z){
   // this return 2F1(a,b,c,z) only in the case z < -1. 
    if( (1+z)>EPS_PREC){
        sprintf(err_msg, "Domain error, something went wrong in %s : z = %f\n", __func__, z);
        ERROR(err_msg);
    }
    if( (c-a)<0 || (c-b)<0 ){// Undefined behaviour (as of now). Check Abramowitz & Stegun, ch15, [15.3]. 
        sprintf(err_msg,"Domain error, in %s something unexpected happenned. : a = %f, b = %f, c = %f, z = %f",__func__,a,b,c,z);
        ERROR(err_msg);
    }
    
    else{
            double m=0; 
            
            if((b-a)>0) m = b - a ; 
            else { m = a - b ; double foo = a; a = b; b = foo;}
            
            long double prefactor1 = pow(-z,-a)/tgamma(a);
            long double sum_n0 = 0;   long double sum_nm1 = 0;
            long double sum_nm2 = 0; long double sum_nm3=0;
            long double sum_n = 0;
            double n = c-b;
            
            if(fabs(m)<0){
                sprintf(err_msg, "m = %f should be positive !",m);
                WARNING(err_msg);
            }
            double sum_n0m1 = 0;
            int conv_flag = 0;
            for (int j = 0; j  < c-b ; j++) {
                sum_n0m1 = sum_n0;
                sum_n0 += Pochhammer(b,j)*pow(-1.,j)*pow(z,-j-m)/(factorial(j)*factorial(j+m)*tgamma(c-b-j)) * (log(-z) + digamma(j+1) + digamma(j+m+1) - digamma(b+j) - digamma(c-b-j));
                if(j>0 && (fabs((sum_n0m1-sum_n0)/sum_n0) < EPS_PREC)) {conv_flag = 1; break;}
            }
            sum_n = sum_n0; // add the term to ensure coherence of convergence criteria.
            int iter = 0;
            do{
                sum_nm3 = sum_nm2 ;
                sum_nm2 = sum_nm1 ;
                sum_nm1 = sum_n;
                sum_n +=  pow(-1,c-b)*Pochhammer(b,n)*pow(z,-n-m) * tgamma(1-c+b+n)/((factorial(n)*factorial(n+m)));
                n = n + 1;
                iter++;

            }while( (iter<4) || ( fabs((sum_n - sum_nm1)/sum_n)>EPS_PREC && fabs((sum_nm1 - sum_nm2)/sum_nm1)>EPS_PREC ));
            if (conv_flag && (fabs((sum_n0 - sum_n)/sum_n)< EPS_PREC)) sum_n = sum_n0; // Go back to 
            if(conv_flag) printf("iter = %d\n conv_flag : %d\n",iter, conv_flag);
            long double prefactor2 = pow(1-z,-a)/(tgamma(b)*tgamma(c-a));
            long double sum_m = 0;
            for (int j = 0; j < m; j++) 
               if(fabs(m)>EPS_PREC) sum_m += Pochhammer(a,j) * factorial(m-j-1)*pow(z,-j)/(factorial(j)*tgamma(c-a-j));
                else{ sum_m = 0; break; }
            return tgamma(c)*(prefactor1*sum_n + prefactor2*sum_m);
        }
}
double hypgeo_reg_B11(double a, double b, double c, double z){
    double m=0; 
    if((b-a)>0) m = b - a ; 
    else { m = a - b ; double foo = a; a = b; b = foo;}
    printf("\nc-a = %f\t m = %f\n", c-a, m) ;
    double prefactor = 0;
    double sum_m = 0;
    if (m>EPS_PREC){
        for (int i = 0; i < m; i++) 
        sum_m+= Pochhammer(a,i)*Pochhammer(b,i)/(Pochhammer(1-m,i)*tgamma(1+i)) * pow(1-z,i);
        prefactor = tgamma(m)/(tgamma(a+m)*tgamma(b+m));
    }

    double res1 = prefactor * sum_m;
    
    double sum_n = 0 ;
    double sum_nm1 = 0; double sum_nm2 = 0; double sum_nm3 = 0;
    double prefactor2 = -pow(-1,m)*pow(1-z,m);
    int n = 0; 
    do {
        sum_nm3 = sum_nm2; 
        sum_nm2 = sum_nm1;
        sum_nm1 = sum_n; 
        long double rho_1 = Pochhammer(a+m,n) *Pochhammer(b+m,n)*Pochhammer(b,m+n)*digamma(b+m+n)/(tgamma(a)+tgamma(b+m+n));
        long double rho_2 = Pochhammer(b+m,n) *Pochhammer(a+m,n)*Pochhammer(a,m+n)*digamma(a+m+n)/(tgamma(b)+tgamma(a+m+n));
       long double sigma = Pochhammer(a+m,n)*Pochhammer(b+m,n)/(tgamma(a)*tgamma(b));
        sum_n += pow(1-z,n)/(tgamma(1+m+n)*tgamma(1+n)) * ( (rho_1+rho_2) + sigma*(log(1-z)-digamma(1+m+n)-digamma(1+n)));
        n++; 
        
    } while ( (n<4) || ( fabs((sum_n - sum_nm1)/sum_n)>EPS_PREC && fabs((sum_nm1 - sum_nm2)/sum_nm1)>EPS_PREC && fabs((sum_nm2-sum_nm3)/sum_nm2)>EPS_PREC) );
    printf("sum_m = %f \t sum_n = %f \nprefactor = %f \t prefactor2 = %f\n", sum_m, sum_n, prefactor, prefactor2);
    
    return (prefactor * sum_m + prefactor2 * sum_n);

}

double hypgeo2F1_Bateman_2_10_eq19(int a, int b, int  c, double z){
    // Analytic continuation taken from Bateman. Conditions : 
    // abs(arg(-z)) < pi , abs(b-a) = m, c = a +m + l, l positive integer.
    int m;
    int flag = 0 ;
    if((b-a)>0) m = b - a ; 
    else { m = a - b ; int foo = a; a = b; b = foo;}    // Now b-a =m>0
    int l = c - a - m;
    // Check : 
    if(l< 0 || m<0){
        sprintf(err_msg, "In %s, c-b=%d or b-a=%d is negative, should be positive integer",__func__, c-b);
        ERROR(err_msg);
    }
    else{
        double res=0;
        // 3 sums : 
        long double prefactor1, prefactor2, prefactor3 ;
        prefactor1 = pow(-1.,m+l)*pow(-z,-a-m);
        prefactor2 = pow(-z,-a-m)/factorial(l+m-1);
        prefactor3 = pow(-z,-a);
        long double sum1 = 0, sum2 = 0, sum3 = 0;
        long double sum1_nm1 =0, sum1_nm2 = 0, sum1_nm3 = 0; // only 2 is ok ? 
        for (int n = l; n< n+1; n++){
            sum1_nm3 = sum1_nm2;
            sum1_nm2 = sum1_nm1;
            sum1_nm1 = sum1;
            sum1 += Pochhammer(a,n+m)*factorial(n-l)/(factorial(n+m)*factorial(n)) * pow(z,-n);
            if ((n>3) && (isnan(sum1))) {
                printf("NAN occured : n = %d \t sum1_nm1 = %.6e\n", n, sum1_nm1);
                flag = 1 ;
                break;
            }
            else if( (n>3) && ( fabs((sum1 - sum1_nm1)/sum1)<EPS_PREC && fabs((sum1_nm1 - sum1_nm2)/sum1_nm1)<EPS_PREC && fabs((sum1_nm2 - sum1_nm3)/sum1_nm2)<EPS_PREC) )break;
            else if(( n > 100 ) && ( fabs((sum1 - sum1_nm1)/sum1)<EPS_PREC && fabs((sum1_nm1 - sum1_nm2)/sum1_nm1)<1e-4) )break;
            
            else if(( n > 300) && ( fabs((sum1 - sum1_nm1)/sum1)<1e-4) )break;
            else if(n>700 && fabs((sum1 - sum1_nm1)/sum1)>1e-3){ flag = 1 ; break; } // If it still doesn't converge, gives NAN
        }
    
   

       if(l!=0){ 
           for (int n = 0; n < l; n++) 
                sum2+= Pochhammer(a,n+m)*Pochhammer(1-m-l,n+m)*pow(z,-n)/(factorial(n)*factorial(n+m)) * (log(-z)+digamma(1+m+n)+digamma(1+n)-digamma(a+m+n)-digamma(l-n)); 
        }
        
        if(m!=0){
            for (int n = 0; n < m; n++) 
                sum3+= factorial(m-n-1)*Pochhammer(a,n)*pow(z,-n)/(factorial(m+l-n-1)*factorial(n));
        }
        
        res = prefactor1*sum1 + prefactor2*sum2 + prefactor3*sum3 ;
        if(flag) {
            sprintf(err_msg, "In %s numerical instability detected.", __func__);
            WARNING(err_msg);
            return NAN;
        }
         return (tgamma(a+m+l)/tgamma(a+m)) * res; // Should check if result should be renorm w/ Gamma(c)
    }

}

double hypgeo2F1(double a, double b, double c, double z){
  // need to handle only the case where z < -1 . z > 1 cannot happen with args. 
  if(z > -1.)  return gsl_sf_hyperg_2F1(a,b,c,z);
  else if ( z <= -1 ){
    int a1 = (int) a;
    int b1 = (int) b;
    int c1 = (int) c;
    return hypgeo2F1_Bateman_2_10_eq19(a1,b1,c1,z);    
  }
  else {
    sprintf(err_msg, "In %s:%d, case not handled :\na = %f, b = %f c = %f \t z= %f", __func__, __LINE__, a,b,c,z);
    ERROR(err_msg);
  }
} 

double P_loop(int i, int j, int k, double a, double b){
    /// Computes the P_ijk(a,b) loop integrals using hypgeo2F1
    double result;
    if(a<0||b<0){
        sprintf(err_msg,"In %s, a = %f and b = %f. They should be positive");
        ERROR(err_msg);
    }
    if ( (i<0)||(j<0)||(k<0)||(a<0)||(b<0)){ 
        sprintf(err_msg,"ERROR : In %s, parameters must be positive !", __func__);
        ERROR(err_msg);
        exit(EXIT_FAILURE);
    }
    else{
            if ( k!=1 ){
                if ( (i==0) && (j==3) && (k==2) && (fabs(a-b)>EPS_PREC) ){
                 return 0; // same as above                                     
                }                                                                  
                else if(fabs(a-b)>EPS_PREC ) return (1./((1-k)*(a-b))) * Beta_E(i,j+1)*( hypgeo2F1(k-1,i,i+j+1,1-a) - hypgeo2F1(k-1,i,i+j+1,1-b) );
                else if (fabs(a-b)<=EPS_PREC) return Beta_E(i+1,j+1)*hypgeo2F1(k,i+1,i+j+2,1-a);
                else{ 
                    ERROR(sprintf(err_msg,"In %s, unknown domain !", __func__)); 
                    return 0;
                }
            }
            else if(k==1){ // Expressions obtained from sympy 
                if ( (i==1) && (j==1) && (k==1) && (fabs(a-b)>EPS_PREC) ){//P_111(a,b)
                  return (1.0/2.0)*(pow(a - 1, 2)*(a - b)*pow(b - 1, 2) + pow(a - 1, 2)*pow(b - 1, 2)*(log(a) - log(b))*(a*b - a - b + 1) + pow(a - 1, 2)*(2*b - 1)*(-log(b/fabs(b - 1)) + log(1.0/fabs(b - 1)))*(a*b - a - b + 1) + (2*a - 1)*pow(b - 1, 2)*(log(a/fabs(a - 1)) - log(1.0/fabs(a - 1)))*(a*b - a - b + 1))/(pow(a - 1, 2)*(a - b)*pow(b - 1, 2)*(a*b - a - b + 1));
                }
                else if( (i==1) && (j==1) && (k==1) && (fabs(a-b)<=EPS_PREC) ){//P_111(a,a)
                    return (1.0/2.0)*(pow(a, 2) - 2*a*log(a) - 1)/(pow(a, 3) - 3.*pow(a, 2) + 3.*a - 1);

                }
                else if ( (i==0) && (j==2 ) && (k==1)&&(fabs(a-b)<=EPS_PREC)) {//P_021(a,a)
                return (1.0/2.0)*(2.*pow(a, 2)*log(a) - 3.*pow(a, 2) + 4.*a - 1)/(pow(a, 3) - 3.*pow(a, 2) + 3.*a - 1); 
                }
                else if ( (i==0) && (j==2) && (k==1) && (fabs(a-b)>EPS_PREC) ){
                    return 0;
                    // This should be called, but not used due to Kroneckers
                }
               
               
                else{
                    sprintf(err_msg,"In %s:%d, unknown domain ! Given args : \n i = %d\nj = %d\nk = %d\na=%f\nb=%f", __func__, __LINE__, i, j , k, a ,b);
                    ERROR(err_msg);
                    exit(EXIT_FAILURE);
                }
            }
    }
}
double f_loop(double a, double b, double c){
    int method =1; // method used to compute f(a,a,c). 0: own calculation, 1: sympy
    if(a<0||b<0||c<0){
        sprintf(err_msg,"In %s, negative masses were give. Contribution is put to zero",__func__);
        WARNING(err_msg);
        return 0;
    }
    else if(a<EPS_PREC||b<EPS_PREC||c<EPS_PREC){
        sprintf(err_msg, "In %s, some input masses are basically zero.", __func__);    
        WARNING(err_msg);                                                         
    }                                          

    if(a<EPS_PREC||fabs(a-1)<EPS_PREC||fabs(a-c)<EPS_PREC||fabs(a-c+1)<EPS_PREC){
        return f_loop(a+2.*EPS_PREC, b,c);
    }
    if(b<EPS_PREC||fabs(b-1)<EPS_PREC||fabs(b-c)<EPS_PREC||fabs(b-c+1)<EPS_PREC){
        return f_loop(a,b+2.*EPS_PREC,c)    ;
    } 
    if ( fabs(c-1)<EPS_PREC||fabs(c-a)<EPS_PREC||fabs(c-b)<EPS_PREC||fabs(c-b-1)<EPS_PREC||fabs(c-a-1)<EPS_PREC ){
        return f_loop(a,b,c+2.*EPS_PREC);
    }
    if(fabs(a-b)<EPS_PREC && method) {

        double res =  0.25*(a*pow(a - 1, 3)*pow(a - c, 2)*pow(c - 1, 2)*(-a + 2.*c - 1)*(a*c - a - c + 1) + pow(c, 2)*pow(a - 1, 3)*log(a/c)*(a*c - a - c + 1)*(-pow(a, 3)*c + pow(a, 3) + pow(a, 2)*pow(c, 2) + pow(a, 2)*c - 2.*pow(a, 2) - 2.*a*pow(c, 2) + a*c + a + pow(c, 2) - c) + pow(a - 1, 3)*pow(a - c, 2)*pow(c - 1, 2)*(pow(a, 3)*c - pow(a, 3) - pow(a, 2)*pow(c, 2) - pow(a, 2)*c + 2.*pow(a, 2) + 2.*a*pow(c, 2) - a*c - a - pow(c, 2) + c) + pow(a - c, 2)*log(a)*(-2.*a*c + a + 1)*(a*c - a - c + 1)*(-pow(a, 3)*c + pow(a, 3) + pow(a, 2)*pow(c, 2) + pow(a, 2)*c - 2.*pow(a, 2) - 2.*a*pow(c, 2) + a*c + a + pow(c, 2) - c))/(pow(a - 1, 3)*pow(a - c, 2)*pow(c - 1, 2)*(a*c - a - c + 1)*(-pow(a, 3)*c + pow(a, 3) + pow(a, 2)*pow(c, 2) + pow(a, 2)*c - 2.*pow(a, 2) - 2.*a*pow(c, 2) + a*c + a + pow(c, 2) - c));
        return res;
    }
    else return (1.0/4.0)*(pow(a - 1, 2)*(a - c)*(-b*(b - 1)*(b - c)*(c - 1) + pow(c, 2)*pow(b - 1, 2)*(log(c/b)) - (b - c)*log(b)*(-2.*b*c + b + c)) + pow(b - 1, 2)*(b - c)*(a*(a - 1)*(a - c)*(c - 1) + pow(c, 2)*pow(a - 1, 2)*(log(a/c)) + (a - c)*(log(a))*(-2*a*c + a + c)))/(pow(a - 1, 2)*pow(a - b, 2)*(a - c)*pow(b - 1, 2)*(b - c)*pow(c - 1, 2)) ;


}
/**************************Loop functions 1205.1500*****************************/
double f_1(double x){
    if( fabs(x) < EPS_PREC){
        return -7./6. ;
    }
    else if (fabs(x-1) < EPS_PREC){
        return -5./12. ;
    }
    else{
        return (-7. + 5.*x + 8.*pow(x,2.))/(6.*pow(1.-x,3.)) - log(x)*(2.*x-3.*pow(x,2.))/(pow(1.-x,4.));
    }
}
double c_term0(double m1, double m2, double m3){
    if(fabs(m1-m2)<EPS_PREC || fabs(m1-m3)<EPS_PREC){
      printf("In %s, something's wrong\nargs : m1 = %f\t m2= %f\t m3=%f\n", __func__, m1, m2, m3);
      exit(EXIT_FAILURE);
    }
    else{
      return (m1 * log(m1/pow(MU_0,2.))) / ((m1 - m2)*(m1 - m3)); 
    }
}
double c_term2(double m1, double m2, double m3){ 
    return pow(m1,2.) * log(m1/pow(MU_0,2.)) / ((m1 - m2)*(m1 - m3)) ;   
}

double d_term2(double m1, double m2, double m3, double m4){
    if( (fabs(m1 - m2) > EPS_PREC) && (fabs(m1 - m3) > EPS_PREC) && (fabs(m1-m4 ) > EPS_PREC)){
        double res= pow(m1,2.) * log(m1/pow(MU_0,2.)) / ((m1 - m2)*(m1 - m3)*(m1 - m4));
        return res;
        }
    
    else{
        sprintf(err_msg,"Something went wrong in %s:%d", __func__,__LINE__);
        ERROR(err_msg);
    }
}

double c_0(double m1, double m2, double m3){
    if (fabs(m2-m3)>EPS_PREC && fabs(m1-m3)>EPS_PREC){
        double res = -1.*( c_term0(m1, m2, m3) + c_term0(m2, m1, m3) + c_term0(m3, m2, m1));
        return res;
    }
    else if ( fabs(m1-m2) < EPS_PREC && fabs(m1-m3)>EPS_PREC ){
        return -1.*(m2 - m3*log(m2) + m3*log(m3) - m3)/(pow(m2, 2) - 2*m2*m3 + pow(m3, 2));
    }
    else if ( fabs(m2 - m3) < EPS_PREC && fabs(m1 - m3)>EPS_PREC)  {
        return -1.* (m1*log(m1) - m1*log(m2) - m1 + m2 )/(pow(m1, 2) - 2*m1*m3 + pow(m3, 2));
    }  
    else if(fabs(m1-m2)<EPS_PREC && fabs(m2-m3)<EPS_PREC) return -1./(2.*m1);
    else{
        printf("In %s, something's wrong\nargs : m1 = %f\t m2= %f\t m3=%f\n", __func__, m1, m2, m3);
        exit(EXIT_FAILURE);
    }
}



double c_2( double m1, double m2, double m3 ){
    if( fabs(m2 - m3) > EPS_PREC && fabs(m1-m3)>EPS_PREC){
        double res = (3./8.) - (1./4) *( c_term2(m1, m2, m3) + c_term2(m2, m1, m3) + c_term2(m3, m2, m1) ); 
    return res;
    }
    else if(fabs(m1)<EPS_PREC||fabs(m2)<EPS_PREC||fabs(m3)<EPS_PREC){
        sprintf(err_msg,"In %s, args cannont <= 0 !",__func__);
        ERROR(err_msg);
    }
    else if (fabs(m2 - m3)< EPS_PREC && fabs(m1-m3)>EPS_PREC){
        double res = c_2(m1, m2+1.2*EPS_PREC,m3);
        return res; 
    }
    else if(fabs(m1-m3)<EPS_PREC && fabs(m2-m3)>EPS_PREC){
        double res = c_2(m1+1.2*EPS_PREC, m2, m3);
        return res;
    }
    else if(fabs(m1-m2)<EPS_PREC && fabs(m2-m3)>EPS_PREC){
        double res = c_2(m1+1.2*EPS_PREC, m2, m3);
        return res ;
    }
    else if (fabs(m1-m2)<EPS_PREC && fabs(m1-m3)<EPS_PREC){
        double res = 1.0*(0.375*m1 - 0.25)/m1;
        return res;
    }
    else{
        sprintf(err_msg,"In %s, something's wrong\nargs : m1 = %f\t m2= %f\t m3=%f\n", __func__, m1, m2, m3);
        ERROR(err_msg);
    }
}


double d_2(double m1, double m2, double m3, double m4){
    
    // All equal
    if (fabs(m1 - m2 - m3 - m4)<EPS_PREC) return -1./4.* 1/(3*m3);
    
    // Three out of four equal
    else if (fabs(m1-m2-m3)<EPS_PREC) return d_2(1.2*EPS_PREC+m1, -1.2*EPS_PREC+m2,m3, m4); // 1,2,3 
    else if (fabs(m1-m2-m4)<EPS_PREC) return d_2(1.2*EPS_PREC+m1,m2,m3,-1.2*EPS_PREC+m4); // 1,2,4
    else if (fabs(m1-m3-m4)<EPS_PREC) return d_2(1.2*EPS_PREC+m1,m2, -1.2*EPS_PREC+m3,m4); // 1,3,4
    else if (fabs(m2-m3-m4)<EPS_PREC) return d_2(m1, m2,1.2*EPS_PREC+m3, -1.2*EPS_PREC+m4); // 2,3,4

    // Two
    else if ( fabs(m1 - m2)<EPS_PREC && fabs(m1-m3)>EPS_PREC && fabs(m1-m4)>EPS_PREC ) // 1,2
        return (-1./4.)*(m2*(2.*log(m2/pow(MU_0,2)) + 1)/((m2-m3)*(m2-m4)) - d_term2(m2, m3, m3, m4) - d_term2(m2,m3, m4, m4) + d_term2(m3, m2,m2,m4) + d_term2(m4,m2, m3, m2) );

    else if( fabs(m1-m2)>EPS_PREC && fabs(m1-m3)<EPS_PREC && fabs(m1-m4)>EPS_PREC ) // 1,3
        return (-1./4.)*(m3*(2.*log(m3/pow(MU_0,2)) + 1)/((m3-m2)*(m3-m4)) - d_term2(m3, m2, m2, m4) - d_term2(m3, m2, m4, m4) + d_term2(m2, m3, m3, m4) + d_term2(m4, m2, m3, m3) );
    
    else if( fabs(m1-m2)>EPS_PREC && fabs(m1-m3)> EPS_PREC && fabs(m1-m4)<EPS_PREC ) // 1,4
        return (-1./4.)*(m4*(2.*log(m4/pow(MU_0,2))+1)/((m4-m2)*(m4-m3)) - d_term2(m4,m2,m2,m3) - d_term2(m4, m2,m3,m3) + d_term2(m2, m4, m3, m4) + d_term2(m3, m2, m4, m4) );

    else if( fabs(m2-m3)<EPS_PREC && fabs(m2-m4) >EPS_PREC && fabs(m2-m1)>EPS_PREC ) // 2,3
        return (-1./4.)*(m3*(2.*log(m3/pow(MU_0,2))+1)/((m3-m1)*(m3-m4)) - d_term2(m3,m1, m1, m4) - d_term2(m3,m4, m4, m1) + d_term2(m1,m3, m3, m4) + d_term2(m4,m3, m3, m1) );
    
    else if( fabs(m2-m4)<EPS_PREC && fabs(m2-m3)>EPS_PREC && fabs(m2-m1)>EPS_PREC) // 2,4
        return (-1./4.)*(m4*(2.*log(m4/pow(MU_0,2))+1)/((m4-m1)*(m4-m3)) - d_term2(m4,m1, m1, m3) - d_term2(m4,m3, m3, m1) + d_term2(m1,m4, m4, m3) + d_term2(m3,m4, m4, m1) );

    else if (fabs(m3-m4)<EPS_PREC)// 3,4 (should be the last possible pair)
        return (-1./4.)*(m4*(2.*log(m4/pow(MU_0,2))+1)/((m4-m1)*(m4-m2)) - d_term2(m4,m1, m1, m2) - d_term2(m4,m2, m2, m1) + d_term2(m1,m4, m4, m2) + d_term2(m2,m4, m4, m1) );
    else if ( fabs(m1-m2)>EPS_PREC &&  fabs(m1-m3)>EPS_PREC &&  fabs(m1-m4)>EPS_PREC && fabs(m2-m3)>EPS_PREC && fabs(m2-m4)>EPS_PREC && fabs(m3-m4)>EPS_PREC ){
        
        double res = (-1./4.)*(d_term2(m1, m2, m3, m4) + d_term2(m2, m1, m3, m4) + d_term2(m3, m2, m1, m4) + d_term2(m4, m2, m3, m1));
        return res;
    }
    else{
        sprintf(err_msg,"In %s:%d, unknown case ! : %f\t%f\t%f\t%f", __func__, __LINE__, m1, m2, m3, m4);
        ERROR(err_msg);
    }
}


double f_7(double x){
    if( fabs(x-1.) > EPS_PREC ){
        return (1.0/6.0)*(43*pow(x, 2) - 101*x + 52)/pow(1 - x, 3) + (2*pow(x, 3) - 9*x + 6)*log(x)/pow(1 - x, 4); 
    }
    else if (fabs(x-1)<=EPS_PREC){
        return -7./4. ;
    }
    else if( fabs(x)<EPS_PREC ){
        sprintf(err_msg,"In %s : infinity occured !", __func__);
        ERROR(err_msg);
    }
    else{
        sprintf(err_msg, "In %s : unknown domain for args : x = %f", __func__, x);
        ERROR(err_msg);
    }
}

double F_Zp(struct parameters * param, double  x[], double x_tr){
    double x_av = (1./6.)*(5.+x_tr);
    double sum = 0;
    if( fabs(x_tr-1)<EPS_PREC) return F_Zp(param, x, x_tr-1.5*EPS_PREC);
    for (int i = 1; i < 3; i++) {
        for (int j = 1; j < 3; j++) {
            sum += param->charg_Vmix[j][1]*conj(param->charg_Vmix[i][2])* ((conj(param->charg_Umix[j][1])*param->charg_Umix[i][1]*sqrt(x[i]*x[j])*( (c_0(x_tr, x[i], x[j]) - c_0(1., x[i], x[j]))/(x_tr - 1.))) - (2.*param->charg_Vmix[i][1]*conj(param->charg_Vmix[j][1])* ( (c_2(x_tr, x[i], x[j])  - c_2(1., x[i], x[j]))/(x_tr - 1.)) ) +( 2.*Krodelta(i,j)* ((c_2(x[j], 1., x_tr) - c_2(x[j], 1., 1.))/(x_tr - 1.) )));
        }

    }
    return  x_av * sum ;
}

double F_gammap(struct parameters* param, double x[], double x_tr){
    
    double x_av = (1./6.)*(5. + x_tr);
    double sum = 0;

    for (int i = 1; i < 3; i++) {
        sum += param->charg_Vmix[i][1]*conj(param->charg_Vmix[i][2])*x_tr*((x[i]/x_tr)*f_7(x[i]/x_tr) - x[i]*f_7(x[i]))/( -x[i]+x[i]/x_tr );
    }
    return (1./9.)*x_av*sum;
}
double F_box(struct parameters* param, double x[], double x_tr, double x_snu){
    double x_av = (1./6.)*(5.+x_tr);
    double sum = 0 ;
    
    if( fabs(x_tr-1)<EPS_PREC) return F_box(param, x, x_tr-1.5*EPS_PREC,x_snu);
    for (int i = 1; i < 3; i++) {
        for (int j = 1; j < 3; j++) {
            sum += param->charg_Vmix[i][1] * conj(param->charg_Vmix[i][2]) * pow(cabs(param->charg_Vmix[j][1]),2.) * ( (d_2(x[i], x[j], x_tr, x_snu) - d_2(x[i], x[j], 1., x_snu))/(x_tr - 1.) );
        }
    }
    return 4.*x_av*sum;
}
/*******************************Wilson Coeffs**********************************/
double C7_gluino_gamma(Input_var * var, struct parameters *param){
    double factor1, factor2;
    double x = pow(var->mass_gluino/var->M_sq,2.);
    factor1 = ( sqrt(2.)/(pow(var->M_sq,2.)*param->Gfermi) * (1./3.)* (4./3.) * pi*param->alphas_MZ/(conj(param->Vts)*param->Vtb));
    factor2 = ((2.*param->delta_d_23_LL + (param->mass_s/param->mass_b)*param->delta_d_23_RR) * (1./4.)*P_loop(1,3,2,x,x)) - param->delta_d_23_RL*P_loop(1,2,2,x,x)*(var->mass_gluino/param->mass_b) ;
    
    return factor1*factor2;    
}

double C9_gluino_gamma(Input_var * var, struct parameters * param){

    double factor1, factor2;
    double x = pow(var->mass_gluino/var->M_sq,2.);
    factor1 = -( sqrt(2.)/(pow(var->M_sq,2.)*param->Gfermi) * (1./3.)* (4./3.) * pi*param->alphas_MZ/(conj(param->Vts)*param->Vtb));
    
    factor2 = (1./3.)*P_loop(0,4,2,x,x)*param->delta_d_23_LL;
    return factor1*factor2;
}

double C9_gluino_Z_2MI(Input_var * var, struct parameters * param){
    
    double x = pow(var->mass_gluino/var->M_sq,2.);
    double res = (4.*SW2-1.)*((param->delta_d_33_LR*param->delta_d_23_RL)/(param->Vtb*conj(param->Vts)))*(4./3.)*(param->alphas_MZ*inv_alpha_MZ/12.)*P_loop(0,3,2,x,x);
    return res;
}

double C10_gluino_2_2MI(Input_var * var, struct parameters * param){
    
    double x = pow(var->mass_gluino/var->M_sq,2.);
    
    double res= ((param->delta_d_33_LR*param->delta_d_23_RL)/(param->Vtb*conj(param->Vts)))*(4./3.)*(param->alphas_MZ/(12.*inv_alpha_MZ))*P_loop(0,3,2,x,x);
    
    return res; 
}

double C7_WWgamma(Input_var * var, struct parameters * param){
    
    double x[] = {0, pow(param->mass_cha1/var->M_sq,2.), pow(param->mass_cha2/var->M_sq,2.)};
    double prefactor = 0;
    double sum = 0;

    prefactor = -param->delta_u_23_LL *pow(param->mass_W/var->M_sq,2.)*(1./3.)*(conj(param->Vcs)/conj(param->Vts));
    for (int i = 1; i < 3; i++) {
    
        sum += (param->charg_Vmix[i][1]*conj(param->charg_Vmix[i][1]))*((3./2.)*P_loop(2,2,2,x[i],x[i]) + P_loop(1,3,2,x[i], x[i]));
    }
    return prefactor*sum;
}


double C9_WWgamma(Input_var * var, struct parameters * param){
    double x[] = {0, pow(param->mass_cha1/var->M_sq,2.), pow(param->mass_cha2/var->M_sq,2.)};
    double prefactor = 0;
    double sum = 0;

    prefactor = -param->delta_u_23_LL *pow(param->mass_W/var->M_sq,2.)*(2./3.)*(conj(param->Vcs)/conj(param->Vts));
    
    for (int i = 1; i < 3; i++) {
        sum += (param->charg_Vmix[i][1]*conj(param->charg_Vmix[i][1]))*(P_loop(3,1,2,x[i], x[i]) - (1./3.)*P_loop(0,4,2,x[i],x[i]) + x[i]*P_loop(3,1,3,x[i],x[i]));
    }
    return prefactor*sum;
}

double C7_HWgamma(Input_var * var, struct parameters * param){
    double x[] = {0, pow(param->mass_cha1/var->M_sq,2.), pow(param->mass_cha2/var->M_sq,2.)};
    double m_cha[] = {0, param->mass_cha1, param->mass_cha2};
    double sum = 0;
    
    for (int i = 1; i < 3; i++) {
        sum += ( ((conj(param->charg_Vmix[i][2])*param->charg_Vmix[i][1]) * (param->yut[3]/param->g2) * ( 0.5*P_loop(2,2,2,x[i],x[i]) + (1./3.)*P_loop(1,3,2,x[i],x[i])) * param->delta_u_23_LR) ) + ( (conj(param->charg_Umix[i][2])*param->charg_Vmix[i][1] * (m_cha[i]/param->mass_b) * (param->yub[3]/param->g2))*(P_loop(2,1,2, x[i], x[i]) + (2./3.)*P_loop(1,2,2,x[i], x[i]))* param->delta_u_23_LL);
        
    }
    return pow(param->mass_W/var->M_sq,2.)*(conj(param->Vcs)/conj(param->Vts))*sum ;
}

double C9_HWgamma(Input_var * var, struct parameters * param){
    double x[] = {0, pow(param->mass_cha1/var->M_sq,2.), pow(param->mass_cha2/var->M_sq,2.)};
    double sum = 0;

    double prefactor = param->delta_u_23_LR * pow(param->mass_W/var->M_sq,2.)*(2./3)*(param->yut[3]/param->g2)*(conj(param->Vcs)/conj(param->Vts));

    for (int i = 1; i < 3; i++) {
        sum += ( conj(param->charg_Vmix[i][2])*param->charg_Vmix[i][1] ) * ( P_loop(3,1,2,x[i],x[i]) - (1./3.)*P_loop(0,4,2,x[i], x[i]) + x[i]*P_loop(3,1,3,x[i],x[i]) );
    }

    return prefactor*sum;
}


double C9_HWZ(Input_var * var, struct parameters * param){
    double prefactor;
    double x[] = { 0, pow(param->mass_cha1/var->M_sq,2.), pow(param->mass_cha2/var->M_sq,2.) };
    double sum = 0;

    prefactor= (4.*SW2 - 1.)*param->delta_u_23_LR*(param->yut[3]/param->g2)*(conj(param->Vcs)/conj(param->Vts))*(1./(4.*SW2));
    for (int i = 1; i < 3; i++) {
        for( int j = 1; j < 3; j++){
            sum += param->charg_Vmix[i][1]*conj(param->charg_Vmix[j][2]) * ( (conj(param->charg_Umix[i][1])*param->charg_Umix[j][1])*sqrt(x[i]*x[j])*P_loop(1,1,2,x[i],x[j]) + (conj(param->charg_Vmix[i][1])*param->charg_Vmix[j][1]*P_loop(1,1,1,x[i],x[j])) - 0.5*Krodelta(i,j)*P_loop(0,2,1,x[i],x[j]));
        }
    }
    return prefactor*sum;
}


double C9_WWZ(Input_var * var, struct parameters * param){
    double prefactor;
    double x[] = { 0, pow(param->mass_cha1/var->M_sq,2.), pow(param->mass_cha2/var->M_sq,2.) };
    double sum = 0;
    prefactor = (1.-4.*SW2)* param->delta_u_23_LL * (conj(param->Vcs)/conj(param->Vts))*(1./4.*SW2);
    for (int i = 1; i < 3; i++) {
        for(int j = 1 ; j < 3 ; j++){

            sum += param->charg_Vmix[i][1]*conj(param->charg_Vmix[j][1]) * ( (conj(param->charg_Umix[i][1])*param->charg_Umix[j][1])*sqrt(x[i]*x[j])*P_loop(1,1,2,x[i],x[j]) + (conj(param->charg_Vmix[i][1])*param->charg_Vmix[j][1]*P_loop(1,1,1,x[i],x[j])) - Krodelta(i,j)*P_loop(0,2,1,x[i],x[j]));
        }
    }
    
    return prefactor * sum ;
}

double C9_boxW(Input_var * var, struct parameters * param){
    double prefactor;
    double x[] = { 0, pow(param->mass_cha1/var->M_sq,2.), pow(param->mass_cha2/var->M_sq,2.) };
    double x_nu = pow(var->mass_sneu/var->M_sq,2.);
    double sum = 0;
    prefactor = param->delta_u_23_LL * (conj(param->Vcs)/conj(param->Vts)) * pow(param->mass_W/var->M_sq,2.)*(1./SW2);
    for (int i = 1; i < 3; i++) {
         for (int j = 1; j < 3; j++) {
             sum += (conj(param->charg_Vmix[i][1]) * param->charg_Vmix[j][1] * param->charg_Vmix[i][1] * conj(param->charg_Vmix[j][1]) ) * f_loop(x[i], x[j], x_nu);
         }
        
    }
    return prefactor * sum ;

}

double C9_boxHW(Input_var * var, struct parameters * param){
    double prefactor;
    double x[] = { 0, pow(param->mass_cha1/var->M_sq,2.), pow(param->mass_cha2/var->M_sq,2.) };
    double x_nu = pow(var->mass_sneu/var->M_sq,2.);
    double sum = 0;
    
    prefactor = - param->delta_u_23_LR * (conj(param->Vcs)/conj(param->Vts)) * (param->yut[3]/(param->g2*4.*SW2))* pow(param->mass_W/var->M_sq,2.);
    for (int i = 1; i < 3; i++) {
         for (int j = 1; j < 3; j++) {
             sum += (conj(param->charg_Vmix[i][1]) * param->charg_Vmix[j][1] * param->charg_Vmix[i][1] * conj(param->charg_Vmix[j][2]) ) * f_loop(x[i], x[j], x_nu);
         }
    }
    printf("In %s, prefactor = %f\t sum=%f\n",__func__, prefactor, sum) ;
    return prefactor*sum;
}

double C9_WWZ_2MI(Input_var * var, struct parameters * param){
    double prefactor;
    double x[] = { 0, pow(param->mass_cha1/var->M_sq,2.), pow(param->mass_cha2/var->M_sq,2.) };
    double sum = 0;
    prefactor = (1-4.*SW2)* param->delta_u_23_LR * param->delta_u_33_LR * (conj(param->Vcs)/conj(param->Vts))*(1./4.*SW2);
    for (int i = 1; i < 3; i++) {
         for (int j = 1; j <3; j++) {
             sum += param->charg_Vmix[i][1]*conj(param->charg_Vmix[j][1]) * ( (conj(param->charg_Umix[i][1])*param->charg_Umix[j][1])*sqrt(x[i]*x[j])*P_loop(1,2,3,x[i],x[j]) + (0.5*conj(param->charg_Vmix[i][1])*param->charg_Vmix[j][1]*P_loop(1,2,2,x[i],x[j])) - (1./3.)*Krodelta(i,j)*P_loop(0,3,2,x[i],x[j]));
            }
         }
    
    return prefactor*sum;
}





//******************************************************* Coeffs 1205.1500


double C7_delta_u23LR(Input_var * var, struct parameters * param){

    double x_tr = pow(param->MtR_Q/var->M_sq,2.);
    double inv_x_tr = 1./x_tr ;
    double x_av = (1./6.) * (5. + x_tr);
    double inv_x[]  = {0, pow(var->M_sq/param->mass_cha1,2.), pow(var->M_sq/param->mass_cha2,2.)};
    double sum = 0 ;
    double prefactor = (conj(param->Vcs)/conj(param->Vts)) * (param->yut[3]/param->g2) * pow(param->mass_W/var->M_sq,2.) * param->delta_u_23_LR * (1./6.) * x_av;
    for (int i = 1; i < 3; i++) {
        double f_1_term = pow(inv_x[i],2.)* (f_1(inv_x[i]/inv_x_tr)-f_1(inv_x[i]))/((inv_x[i]/inv_x_tr)-inv_x[i]) ;
        sum  += (param->charg_Vmix[i][1]*(param->charg_Vmix[i][2])) * f_1_term;
    }
    return prefactor*sum;
}

double C7_HWgamma_wo_delta_u_23LR(Input_var * var, struct parameters * param){
    double x[] = {0, pow(param->mass_cha1/var->M_sq,2.), pow(param->mass_cha2/var->M_sq,2.)};
    double m_cha[] = {0, param->mass_cha1, param->mass_cha2};
    double sum = 0;
    for (int i = 1; i < 3; i++) {
        sum += ( (param->charg_Umix[i][2]*param->charg_Vmix[i][1] * (m_cha[i]/param->mass_b) * (param->yub[3]/param->g2))*(P_loop(2,1,2, x[i], x[i]) + (2./3.)*P_loop(1,2,2,x[i], x[i]))*param->delta_u_23_LL);
        
    }
    return pow(param->mass_W/var->M_sq,2.)*(conj(param->Vcs)/conj(param->Vts))*sum ;
}

double C9_delta_u23LR_HWZ(Input_var * var, struct parameters * param){
    double prefactor = (conj(param->Vcs)/conj(param->Vts))*(param->yut[3]/param->g2)* (1./(4.*SW2)) * (4.*SW2 - 1.)*param->delta_u_23_LR  ; 
    double x[]= {
        0,
       pow(param->mass_cha1/var->M_sq,2.),
       pow(param->mass_cha2/var->M_sq,2.)
    };
    double x_tr = pow(param->MtR_Q/var->M_sq,2.);
    return prefactor*F_Zp(param, x, x_tr);
}

double C9_delta_u23LR_HWgamma(Input_var * var , struct parameters * param){

    double prefactor = (conj(param->Vcs)/conj(param->Vts))*(param->yut[3]/param->g2) * pow(param->mass_W/var->M_sq,2.)*param->delta_u_23_LR ; 
    double x[]= {
        0,
       pow((param->mass_cha1/var->M_sq),2.),
       pow((param->mass_cha2/var->M_sq),2.)
    };
    double x_tr = pow(param->MtR_Q/var->M_sq,2.);
    return prefactor * F_gammap(param, x, x_tr);

}

double C9_delta_u23LR_boxHW(Input_var * var, struct parameters * param){

    double x[]= {
        0,
       pow((param->mass_cha1/var->M_sq),2.),
       pow((param->mass_cha2/var->M_sq),2.)
    };

    double x_tr = pow(param->MtR_Q/var->M_sq,2.);
    double x_snu = pow(var->mass_sneu/var->M_sq,2.);
    double prefactor = -1.* (conj(param->Vcs)/conj(param->Vts))*(1./(4.*SW2))*(param->yut[3]/param->g2) * pow(param->mass_W/var->M_sq,2.) * param->delta_u_23_LR ; 
    double res = prefactor * F_box(param, x, x_tr, x_snu);
    return res;

}

double C10_delta_u23LR_HWZ(Input_var * var, struct parameters *param ){
        
    double prefactor; 
    double x[] = {
        0,
        pow((param->mass_cha1/var->M_sq),2.),    
        pow((param->mass_cha2/var->M_sq),2.)
    };
    double x_tr = pow(param->MtR_Q/var->M_sq,2.);
    prefactor = conj(param->Vcs)/conj(param->Vts) * (1./(4.*SW2)) * (param->yut[3]/param->g2) * param->delta_u_23_LR;
       
    return F_Zp(param, x, x_tr) * prefactor ;
}

double C10_delta_u23LR_boxHW(Input_var * var, struct parameters * param){
   
    double x[] = {
        0,
        pow((param->mass_cha1/var->M_sq),2.),    
        pow((param->mass_cha2/var->M_sq),2.)
    };
    double x_snu = pow(var->mass_sneu/var->M_sq,2.);
    double x_tr = pow(param->MtR_Q/var->M_sq,2.);
    double prefactor =  conj(param->Vcs)/conj(param->Vts) * (1./(4.*SW2)) * (param->yut[3]/param->g2) * pow(param->mass_W/var->M_sq,2.) * param->delta_u_23_LR;
    return prefactor * F_box(param, x, x_tr, x_snu);
}
/**************** Coefficient calculator****************************/
void CW_calculator_NMFV_1205_1500(Input_var * var , struct parameters * param, double * CW_NMFV){
//    printf("Calculating C7...\n");
    double C7_sum = 0;
    double C7_NMFV[]={
        C7_gluino_gamma(var, param),
        C7_WWgamma(var, param),
        C7_HWgamma_wo_delta_u_23LR(var, param),
        C7_delta_u23LR(var, param)
    };
    char* C7_names[]={
        "gl-gamma",
        "WWgamma",
        "HWgamma_no_delta",
        "delta_u23LR"
    };
      
//    printf("Calculating C9...\n");
    double C9_sum = 0;
    double C9_NMFV[] = {
        C9_gluino_gamma(var, param),
        C9_WWgamma(var, param),
        C9_delta_u23LR_HWgamma(var, param),
        C9_WWZ_2MI(var, param),
        C9_gluino_Z_2MI(var, param),
        C9_delta_u23LR_HWZ(var, param),
        C9_WWZ(var, param),
        C9_boxW(var, param),
        C9_delta_u23LR_boxHW(var, param)
    };
    char* C9_names[] = {
        "gl-gamma",
        "WWgamma",
        "HWgamma_delta",
        "WWZ-2MI",
        "gluino-Z 2MI",
        "HWZ_delta",
        "WWZ-1MI",
        "boxW",
        "boxHW_delta"
    };
    double C10_sum = 0;
    double C10_NMFV[] = {
    C10_delta_u23LR_HWZ(var, param),
    C10_delta_u23LR_boxHW(var, param),
    -1.*C9_WWZ(var, param)/(1.-4.*SW2),
    -1.*C9_boxW(var, param),
    -1.*C9_WWZ_2MI(var, param)/(1.-4.*SW2),
    C10_gluino_2_2MI(var,param)
    };
    char* C10_names[] = {
        "Z-peng",
        "boxHW",
        "WWZ-1MI",
        "boxW",
        "WWZ-2MI",
        "gl-Z-2MI"
    };
    for (int i = 0; i < 9 ; i++) {
        if(i < 4){   C7_sum += C7_NMFV[i];
//        printf("C7:%s = %f \n",C7_names[i], C7_NMFV[i] );
        }
        C9_sum += C9_NMFV[i];
//        printf("C9:%s = %f \n",C9_names[i], C9_NMFV[i] );
    }
    for (int i = 0; i <6 ; i++) {
        C10_sum += C10_NMFV[i];
//        printf("C10:%s = %f \n",C10_names[i], C10_NMFV[i] );
    }
   CW_NMFV[7] += C7_sum;
   CW_NMFV[9] += C9_sum;
   CW_NMFV[10] += C10_sum;
//   
//    printf("C9 = %f\n",CW_NMFV[9]);
}
void CW_calculator_NMFV(Input_var * var , struct parameters * param, double * CW_NMFV){
//    printf("Calculating C7...\n");
    double C7_sum = 0;
    double C7_NMFV[]={
        C7_gluino_gamma(var, param),
        C7_WWgamma(var, param),
        C7_HWgamma_wo_delta_u_23LR(var, param),
        C7_delta_u23LR(var, param)
    };
    char* C7_names[]={
        "gl-gamma",
        "WWgamma",
        "HWgamma_no_delta",
        "delta_u23LR"
    };
      
//    printf("Calculating C9...\n");
    double C9_sum = 0;
    double C9_NMFV[] = {
        C9_gluino_gamma(var, param),
        C9_WWgamma(var, param),
        C9_delta_u23LR_HWgamma(var, param),
        C9_WWZ_2MI(var, param),
        C9_gluino_Z_2MI(var, param),
        C9_delta_u23LR_HWZ(var, param),
        C9_WWZ(var, param),
        C9_boxW(var, param),
        C9_delta_u23LR_boxHW(var, param)
    };
    char* C9_names[] = {
        "gl-gamma",
        "WWgamma",
        "HWgamma_delta",
        "WWZ-2MI",
        "gluino-Z 2MI",
        "HWZ_delta",
        "WWZ-1MI",
        "boxW",
        "boxHW_delta"
    };
    double C10_sum = 0;
    double C10_NMFV[] = {
    C10_delta_u23LR_HWZ(var, param),
    C10_delta_u23LR_boxHW(var, param),
    C10_gluino_2_2MI(var, param)
    };
    char* C10_names[] = {
        "Z-peng",
        "boxHW",
        "g-2MI"
    };
    printf("Done.\n");
    printf("====\t\t Coeff[Graph] = Value\t\t====\n\n");
    for (int i = 0; i < 9 ; i++) {
        if(i < 2){
            printf("C10[%s] = \t %f\n", C10_names[i], C10_NMFV[i]);
            C10_sum += C10_NMFV[i];
        }
        if(i < 3) {
            printf("C7[%s] =\t %f\t C9[%s] = %f\n", C7_names[i], C7_NMFV[i], C9_names[i], C9_NMFV[i]);
            C7_sum += C7_NMFV[i];
        }
        else if( i == 3 ){
            printf("C7[%s] = \t %f\nC9[%s] = %f\n", C7_names[i], C7_NMFV[i], C9_names[i], C9_NMFV[i]);
            C7_sum += C7_NMFV[i];
        }
        else printf("C9[%s] =\t %f\n", C9_names[i], C9_NMFV[i]);
        C9_sum += C9_NMFV[i];
    }
    printf("\n Contribution sum : \n");
    printf("C10_NMFV = %f\n", C10_sum);
    printf("C9_NMFV = %f\n", C9_sum);
    CW_NMFV[7] += C7_sum;
    printf("C7_NMFV = %f\n", C7_sum);
    CW_NMFV[9] += C9_sum;
    CW_NMFV[10] += C10_sum;
}
