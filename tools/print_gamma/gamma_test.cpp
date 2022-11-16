#include <math.h>

#include <nissa.hpp>

using namespace nissa;

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  if(rank==0)
    {
      //check of the first 0-5 gamma
      for(int i=0;i<=5;i++)
	{
	  printf("Gamma%d = %s\n",i,gtag[i]);
	  print_dirac(base_gamma[i]);
	  printf("\n");
	}
      
      //check of gamma5
      dirac_matr gamma5;
      gamma5=base_gamma[4]*base_gamma[1]*base_gamma[2]*base_gamma[3];
      printf("Gamma5 = Gamma4Gamma1Gamma2Gamma3  = %s as:\n",gtag[5]);
      print_dirac(gamma5);
      printf("\n");
      
      //check of gamma6
      dirac_matr gamma6=base_gamma[1]*base_gamma[5];
      printf("Gamma6 = Gamma1Gamma5  = %s as:\n",gtag[6]);
      printf("(a) defined by Silvano\n");
      print_dirac(base_gamma[6]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(gamma6);
      printf("\n");
      
      //check of gamma7
      dirac_matr gamma7=base_gamma[2]*base_gamma[5];
      printf("Gamma7 = Gamma2Gamma5  = %s as:\n",gtag[7]);
      printf("(a) defined by Silvano\n");
      print_dirac(base_gamma[7]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(gamma7);
      printf("\n");
      
      //check of gamma8
      dirac_matr gamma8=base_gamma[3]*base_gamma[5];
      printf("Gamma8 = Gamma3Gamma5  = %s as:\n",gtag[8]);
      printf("(a) defined by Silvano\n");
      print_dirac(base_gamma[8]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(gamma8);
      printf("\n");
      
      //check of gamma9
      dirac_matr gamma9=base_gamma[4]*base_gamma[5];
      printf("Gamma9 = Gamma4Gamma5  = %s as:\n",gtag[9]);
      printf("(a) defined by Silvano\n");
      print_dirac(base_gamma[9]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(gamma9);
      printf("\n");
      
      //check of gamma10
      dirac_matr gamma10=base_gamma[4]*base_gamma[1];
      printf("Gamma10 = Gamma4Gamma1  = %s as:\n",gtag[10]);
      printf("(a) defined by Silvano\n");
      print_dirac(base_gamma[10]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(gamma10);
      printf("\n");
      
      //check of gamma11
      dirac_matr gamma11=base_gamma[4]*base_gamma[2];
      printf("Gamma11 = Gamma4Gamma2  = %s as:\n",gtag[11]);
      printf("(a) defined by Silvano\n");
      print_dirac(base_gamma[11]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(gamma11);
      printf("\n");
      
      //check of gamma12
      dirac_matr gamma12=base_gamma[4]*base_gamma[3];
      printf("Gamma12 = Gamma4Gamma3  = %s as:\n",gtag[12]);
      printf("(a) defined by Silvano\n");
      print_dirac(base_gamma[12]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(gamma12);
      printf("\n");
      
      //check of gamma13
      dirac_matr gamma13=base_gamma[2]*base_gamma[3];
      printf("Gamma13 = Gamma2Gamma3  = %s as:\n",gtag[13]);
      printf("(a) defined by Silvano\n");
      print_dirac(base_gamma[13]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(gamma13);
      printf("\n");
      
      //check of gamma14
      dirac_matr gamma14=base_gamma[3]*base_gamma[1];
      printf("Gamma14 = Gamma3Gamma1  = %s as:\n",gtag[14]);
      printf("(a) defined by Silvano\n");
      print_dirac(base_gamma[14]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(gamma14);
      printf("\n");
      
      //check of gamma15
      dirac_matr gamma15=base_gamma[1]*base_gamma[2];
      printf("Gamma15 = Gamma1Gamma2  = %s as:\n",gtag[15]);
      printf("(a) defined by Silvano\n");
      print_dirac(base_gamma[15]);
      printf("(b) obtained by making explicite products\n");
      print_dirac(gamma15);
      printf("\n");
      
      complex one_over_rad2={1/sqrt(2),0};
      
      //check of Pplus
      complex c={0,1};
      dirac_matr gamma_Pplus=dirac_summ(dirac_compl_prod(base_gamma[5],c),base_gamma[0]);
      gamma_Pplus=dirac_compl_prod(gamma_Pplus,one_over_rad2);
      printf("Pminus=(1+iGamma5)/sqrt(2) as:\n");
      printf("(a) defined by Silvano\n");
      print_dirac(Pplus);
      printf("(b) obtained by making explicite computations\n");
      print_dirac(gamma_Pplus);
      printf("\n");
      
      //check of Pminus
      dirac_matr gamma_Pminus;
      complex mI={0,-1};
      gamma_Pminus=dirac_compl_prod(base_gamma[5],mI);
      gamma_Pminus=dirac_summ(gamma_Pminus,base_gamma[0]);
      gamma_Pminus=dirac_compl_prod(gamma_Pminus,one_over_rad2);
      printf("Pminus=(1-iGamma5)/sqrt(2) as:\n");
      printf("(a) defined by Silvano\n");
      print_dirac(Pminus);
      printf("(b) obtained by making explicite computations\n");
      print_dirac(gamma_Pminus);
      printf("\n");
      
      printf("Pmunu\n");
      char P[4][4][2][50];
      for(int id=0;id<4;id++)
	for(int jd=0;jd<4;jd++)
	  for(int ri=0;ri<2;ri++)
	    strcpy(P[id][jd][ri],"");
      for(int imunu=0;imunu<6;imunu++)
	for(int id=0;id<4;id++)
	  {
	    int jd=smunu_pos[id][imunu];
	    for(int ri=0;ri<2;ri++)
	      if(smunu_entr[id][imunu][ri])
		{
		  char temp[50];
		  snprintf(temp,50,"%s%sF%d",P[id][jd][ri],(((int)smunu_entr[id][imunu][ri]==+1)?(strlen(P[id][jd][ri])?"+":""):"-"),imunu);
		  strcpy(P[id][jd][ri],temp);
		}
	  }
      
      for(int id=0;id<4;id++)
	{
	  for(int jd=0;jd<4;jd++) for(int ri=0;ri<2;ri++) if(strlen(P[id][jd][ri])==0) strcpy(P[id][jd][ri],"0");
	  printf("(%s,%s)",P[id][0][0],P[id][0][1]);
	  for(int jd=1;jd<4;jd++)
	    printf(",\t(%s,%s)",P[id][jd][0],P[id][jd][1]);
	  printf("\n");
	}
    }
    
  close_nissa();

  return 0;
}
