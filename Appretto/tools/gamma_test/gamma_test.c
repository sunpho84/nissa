#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

int main(int narg,char **arg)
{

  init_appretto();

  if(rank==0)
    {
      //check of the first 0-5 gamma
      for(int i=0;i<=5;i++)
	{
	  printf("Gamma%d = %s\n",i,gtag[i]);
	  print_dirac(&(base_gamma[i]));
	  printf("\n");
	}

      //check of gamma6
      dirac_matr gamma6;
      dirac_prod(&gamma6,&(base_gamma[1]),&(base_gamma[5]));
      printf("Gamma6 = Gamma1Gamma5  = %s as:\n",gtag[6]);
      printf("(a) defined by Silvano\n");
      print_dirac(&(base_gamma[6]));
      printf("(b) obtained by making explicite products\n");
      print_dirac(&gamma6);
      printf("\n");

      //check of gamma7
      dirac_matr gamma7;
      dirac_prod(&gamma7,&(base_gamma[2]),&(base_gamma[5]));
      printf("Gamma7 = Gamma2Gamma5  = %s as:\n",gtag[7]);
      printf("(a) defined by Silvano\n");
      print_dirac(&base_gamma[7]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(&gamma7);
      printf("\n");

      //check of gamma8
      dirac_matr gamma8;
      dirac_prod(&gamma8,&(base_gamma[3]),&(base_gamma[5]));
      printf("Gamma8 = Gamma3Gamma5  = %s as:\n",gtag[8]);
      printf("(a) defined by Silvano\n");
      print_dirac(&(base_gamma[8]));
      printf("(b) obtained by making explicite products\n");
      print_dirac(&gamma8);
      printf("\n");

      //check of gamma9
      dirac_matr gamma9;
      dirac_prod(&gamma9,&(base_gamma[4]),&(base_gamma[5]));
      printf("Gamma9 = Gamma4Gamma5  = %s as:\n",gtag[9]);
      printf("(a) defined by Silvano\n");
      print_dirac(&(base_gamma[9]));
      printf("(b) obtained by making explicite products\n");
      print_dirac(&gamma9);
      printf("\n");
      
      //check of gamma10
      dirac_matr gamma10;
      dirac_prod(&gamma10,&(base_gamma[4]),&(base_gamma[1]));
      printf("Gamma10 = Gamma4Gamma1  = %s as:\n",gtag[10]);
      printf("(a) defined by Silvano\n");
      print_dirac(&(base_gamma[10]));
      printf("(b) obtained by making explicite products\n");
      print_dirac(&gamma10);
      printf("\n");

      //check of gamma11
      dirac_matr gamma11;
      dirac_prod(&gamma11,&(base_gamma[4]),&(base_gamma[2]));
      printf("Gamma11 = Gamma4Gamma2  = %s as:\n",gtag[11]);
      printf("(a) defined by Silvano\n");
      print_dirac(&(base_gamma[11]));
      printf("(b) obtained by making explicite products\n");
      print_dirac(&gamma11);
      printf("\n");

      //check of gamma12
      dirac_matr gamma12;
      dirac_prod(&gamma12,&(base_gamma[4]),&(base_gamma[3]));
      printf("Gamma12 = Gamma4Gamma3  = %s as:\n",gtag[12]);
      printf("(a) defined by Silvano\n");
      print_dirac(&(base_gamma[12]));
      printf("(b) obtained by making explicite products\n");
      print_dirac(&gamma12);
      printf("\n");

      //check of gamma13
      dirac_matr gamma13;
      dirac_prod(&gamma13,&(base_gamma[2]),&(base_gamma[3]));
      printf("Gamma13 = Gamma2Gamma3  = %s as:\n",gtag[13]);
      printf("(a) defined by Silvano\n");
      print_dirac(&(base_gamma[13]));
      printf("(b) obtained by making explicite products\n");
      print_dirac(&gamma13);
      printf("\n");

      //check of gamma14
      dirac_matr gamma14;
      dirac_prod(&gamma14,&(base_gamma[3]),&(base_gamma[1]));
      printf("Gamma14 = Gamma3Gamma1  = %s as:\n",gtag[14]);
      printf("(a) defined by Silvano\n");
      print_dirac(&(base_gamma[14]));
      printf("(b) obtained by making explicite products\n");
      print_dirac(&gamma14);
      printf("\n");

      //check of gamma15
      dirac_matr gamma15;
      dirac_prod(&gamma15,&(base_gamma[1]),&(base_gamma[2]));
      printf("Gamma15 = Gamma1Gamma2  = %s as:\n",gtag[15]);
      printf("(a) defined by Silvano\n");
      print_dirac(&(base_gamma[15]));
      printf("(b) obtained by making explicite products\n");
      print_dirac(&gamma15);
      printf("\n");


      complex one_over_rad2={1/sqrt(2),0};

      //check of Pplus
      dirac_matr gamma_Pplus;
      dirac_compl_prod(&gamma_Pplus,&(base_gamma[5]),I);
      dirac_summ(&gamma_Pplus,&gamma_Pplus,&(base_gamma[0]));
      dirac_compl_prod(&gamma_Pplus,&gamma_Pplus,one_over_rad2);
      printf("Pminus=(1+iGamma5)/sqrt(2) as:\n");
      printf("(a) defined by Silvano\n");
      print_dirac(&Pplus);
      printf("(b) obtained by making explicite computations\n");
      print_dirac(&gamma_Pplus);
      printf("\n");

      //check of Pminus
      dirac_matr gamma_Pminus;
      complex mI={0,-1};
      dirac_compl_prod(&gamma_Pminus,&(base_gamma[5]),mI);
      dirac_summ(&gamma_Pminus,&gamma_Pminus,&(base_gamma[0]));
      dirac_compl_prod(&gamma_Pminus,&gamma_Pminus,one_over_rad2);
      printf("Pminus=(1-iGamma5)/sqrt(2) as:\n");
      printf("(a) defined by Silvano\n");
      print_dirac(&Pminus);
      printf("(b) obtained by making explicite computations\n");
      print_dirac(&gamma_Pminus);
      printf("\n");

    }
  
  return 0;
}
