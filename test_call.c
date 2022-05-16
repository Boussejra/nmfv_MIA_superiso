#include "src/include.h"
#define slha_chi_path "./slha_chi2.x"
int main(int argc, char* argv[]){
    char name[500];
    int test; 

    if(argc<2)  { 
        printf(" This program needs 1 parameter:\n"
                "   name    name of the SLHA file\n");
        exit(1); 
    } 
    else{
    	sscanf(argv[1],"%s",name);
    }
                                                   
//        sprintf(tmp_char, "%s leshouches < %s > %s", SOFTSUSY_PATH, fileLesHouches, name);
    char tmp_char[5000];
    sprintf(tmp_char, "%s %s",slha_chi_path, name);
    double chi  = system(tmp_char);
    printf("chi2  = %f\n", chi);
    return 0;
}
// Test inconclusive. Must investigate nature of return type of main.
