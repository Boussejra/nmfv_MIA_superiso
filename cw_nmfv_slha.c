// This file is an example of script to calculate the Wilson coefficients
// C7,C9,C10 in the NMFV-MSSM using the Mass Insertion Approximation.
#include "src/include.h"
#include <errno.h>

// Functions defined at the end.
void calc_wcf(struct parameters * param, Input_var * nmfv_param, char * name);
int lha_reader_init(struct parameters * param, Input_var * nmfv_var);
int test_deltas(struct parameters * param);

int main(int argc, char* argv[]){
  errno = 0;
  if(argc !=2){
    printf("This program takes exactly 1 parameter\n"
        "SLHA1_file \t\t SLHA1 file with and additional BLOCK NMFVPAR:\n"
        "BLOCK NMFVPAR\t   #NMFV params\n"
         "1      val     # delta_u_23_LR \n"
         "2      val     # delta_u_23_LL \n"
         "3      val     # delta_u_33_LR \n"
         "4      val     # delta_d_23_LL \n"
         "5      val     # delta_d_23_RR \n"
         "6      val     # delta_d_23_RL \n"
         "7      val     # delta_d_23_LR \n"
         "8      val     # delta_d_33_RL \n"
         "9      val     # delta_d_33_LR \n"
         "with val a floating point value between [-1,1]\n");
    exit(EXIT_FAILURE);
  }
  char filename[200];
  sscanf(argv[1],"%s",filename);
  if(errno!=0){
    printf("%s",strerror(errno));
    exit(EXIT_FAILURE);
  }
  printf("Treating file %s...\n",filename);

  struct parameters param; Init_param(&param);
  Input_var nmfv_param; init_Inputvar(&nmfv_param);
  
  strcpy(nmfv_param.filename,filename);

  if(lha_reader_init(&param,&nmfv_param)){
    printf("Valid File ! Calculating...\n");
    char * outfile = "out.txt";
    calc_wcf(&param, &nmfv_param, outfile);
    printf("Output written in out.txt\n");
  }
  else{
    printf("The input file you gave was rejected for some reason. Check LSP & excluded_masses conditions.\n");
    exit(EXIT_FAILURE);
  }
  
  return 0;
}

void calc_wcf(struct parameters * param, Input_var *nmfv_param, char *name){
  printf("Calculating wcf and saving...\n");
 	double complex C0w[11],C1w[11],C2w[11];
  double complex C0b[11],C1b[11],C2b[11];
 	double mu_W, mu_b;

  // Calculate complete contribution and run to m_b 
 	mu_W=2*param->mass_W;
  mu_b = param->mass_b_pole ;
  double C0w_NMFV[11] = {0};
 	
  CW_calculator(2,C0w,C1w,C2w,mu_W,param); /* 2 = muon */
  CW_calculator_NMFV_1205_1500(nmfv_param, param,C0w_NMFV);
  
  for (int i = 0; i < 11; i++) C0w[i] += C0w_NMFV[i] + 0.*I;
  
  printf("C9w_p = %.10f \n", creal(C0w[9]));
  printf("C9w_MIA = %.10f \n", C0w_NMFV[9]);
  
  C_calculator_base1(C0w, C1w, C2w,mu_W, C0b, C1b ,C2b, mu_b,param );

	double alphas_mb=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
  double cw[11] = {0}; double cw_b[11] = {0};
  
  for(int i =0; i < 11; i++){
      cw[i] = creal(C0w[i]) + (param->alphas_MZ/(4.*pi))*creal(C1w[i]) + pow(param->alphas_MZ/(4.*pi),2.)*creal(C2w[i]);
      cw_b[i] = creal(C0b[i]) + (alphas_mb/(4.*pi))*creal(C1b[i]) + pow(alphas_mb/(4.*pi),2.)*creal(C2b[i]);
  }

  FILE * fname = fopen(name, "a");
  fprintf(fname,"   C_7(mu_W)    C_9(mu_W)     C_10(mu_W)    C_7(mu_b)     C_9(mu_b)     C_10(mu_b)\n");
  fprintf(fname,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
      cw[7], cw[9], cw[10], cw_b[7], cw_b[9], cw_b[10]
      );
  fclose(fname);
}

int lha_reader_init(struct parameters * param, Input_var * nmfv_param){
  char  msg[200];
  if(Les_Houches_Reader(nmfv_param->filename, param)
     && !(charged_LSP(param)) 
     && !excluded_masses(param) 
     && !(test_deltas(param))){

      printf("LesHouches tests passed. Initializing parameters\n");
       nmfv_param->mu = param->mu_Min;
       nmfv_param->mass_sneu = param->mass_nuel;
       nmfv_param->tan_beta = param->tan_beta;
       nmfv_param->M_sq=(param->mass_dnl 
                      +param->mass_dnr 
                      +param->mass_upl 
                      +param->mass_upr 
                      +param->mass_stl 
                      +param->mass_str 
                      +param->mass_chl 
                      +param->mass_chr)/8.;
       nmfv_param->mass_gluino = param->M3_Min;
       nmfv_param->M2 = param->M2_Min;
    return 1;
  }
  else{
    sprintf(msg,"%s is not a valid SLHA1 file!\n", nmfv_param->filename);
    WARNING(msg);
    return 0;
  }
 return 0;
}

int test_deltas(struct parameters * param){
    if(
        (param->delta_u_23_LR == 1)
    &&  (param->delta_u_23_LL == 1)
    &&  (param->delta_u_33_LR == 1)    
    &&  (param->delta_d_23_LL == 1)
    &&  (param->delta_d_23_RR == 1)
    &&  (param->delta_d_23_RL == 1)
    &&  (param->delta_d_23_LR == 1)
    &&  (param->delta_d_33_RL == 1)
    &&  (param->delta_d_33_LR == 1)){
      printf("You either either initialized all MI parameters to 1 or forgot the NMFVPAR block in the input file\n");
      return 0;
    }
    else return 1;
  }
