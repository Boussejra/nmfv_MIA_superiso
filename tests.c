#include "src/include.h"

int main(){
  // Test files for functions : 
  
  // abs(z) > 1
  int tests = 0; int nb_tests = 0;
  double a[] = {2,2,4,5,8,3,6};
  double b[] = {3,2,5,1,2,9,6};
  double c[] =  {5,4,8,5,9,12,65};
  double z[] = {-14,-3,-1.01,-2,-2,-3, -0.8990};
  double res_hypgeotests[] = {0.0175084685, 0.2069937346,
    0.1159282299, 0.333333333, 0.1319881048, 0.0318181955,0.6307879488}; 
  printf("Tests for z < 1\n");
  for (int i = 0; i <= 6; i++){ 
    nb_tests++;
//    printf("2F1( %.3f, %.3f ; %.3f ; %.3f) =\t %.10f \n", a[i],b[i],c[i],z[i], hypgeo2F1(a[i],b[i],c[i],z[i]));
    if(fabs(hypgeo2F1(a[i],b[i],c[i],z[i]) - res_hypgeotests[i])<1e-9){
      printf("Hypgeo test %d ok.\n",i); tests++ ;
    }
  }
  printf("2F1(-1/2, 1/3, 4/3, -1) = %.10f\n", hypgeo2F1(-1./2., 1./3., 4./3., -1.));
  printf("2F1(-2,3,5,-1.1) = %.10f\n", hypgeo2F1(-2.,3.,5.,-1.1));
  
  printf("Tests passed : %d/%d \n", tests, nb_tests);
  // Tests wilson coeffs Silvestrini : 

  // Should load params from benchmark points given in paper and compare  output. 


  Input_var nmfv_var;
  struct parameters param;

  Init_param(&param);
  slha_adjust(&param);

  Les_Houches_Reader("1205_1500.lha", &param);
  sprintf(nmfv_var.filename,"1205_1500.lha") ;
  set_yukawa(&nmfv_var, &param);
//  nmfv_var.M_sq=(param.mass_dnl+param.mass_dnr+param.mass_upl+param.mass_upr+param.mass_stl+param.mass_str+param.mass_chl+param.mass_chr+param.mass_b1+param.mass_b2+param.mass_t1+param.mass_t2)/12;
  nmfv_var.M_sq = 1000;
  nmfv_var.mass_gluino = param.M3_Q;
  nmfv_var.mass_sneu = param.mass_nuel;

  double CW_NMFV[11] = {0};
  CW_calculator_NMFV_1205_1500(&nmfv_var, &param, CW_NMFV);
  printf("C9 = %f\n", CW_NMFV[9]);
  printf("C7 = %f\n", CW_NMFV[7]);
  printf("C10 = %f\n", CW_NMFV[10]);
    printf("Doing RGE to m_b scale...\n");
    double mu_W,mu_b,mu_spec;	
    double complex C0b[11],C0spec[11],C1b[11],C1spec[11],C0w[11],C1w[11],C2w[11],C2b[11],Cpb[11],CQpb[3],CQ0b[3],CQ1b[3] = {0};
    double complex C0eb[11],C1eb[11],C0ew[11],C1ew[11],C2ew[11],C2eb[11],Cpeb[11],CQpeb[3],CQ0eb[3],CQ1eb[3];              
    mu_W = 120; 
    mu_b = param.mass_b_pole;
    
    // Need to RGE only the NMFV coeffs 
    CW_calculator(2,C0w,C1w,C2w,mu_W,&param); /* 2 = muon */ 
    // C0w[7] = 0;
    // C0w[9] = 0;
    // C0w[10] = 0;
    // C1w[7] = 0;
    // C1w[9] = 0;
    // C1w[10] = 0;
    // C2w[7] = 0;
    // C2w[9] = 0;
    // C2w[10] = 0;
    for (int i = 0; i < 11; i++)  C0w[i] += CW_NMFV[i];
    C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,5.3,&param);
    printf("C_9(mb) = %f \n", creal(C0b[9]));
    printf("C_10(mb) = %f\n", creal(C0b[10]));
    printf("C_7(mb) = %f\n", creal(C0b[7]));




  return 0;
}
