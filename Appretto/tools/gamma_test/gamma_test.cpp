#include <mpi.h>
#include <fstream>
#include <lemon.h>

#include "appretto.h"

using namespace std;

int main(int narg,char **arg)
{

  init_appretto();

  if(rank==0)
    {
      //check of the first 0-5 gamma
      for(int i=0;i<=5;i++)
	{
	  cout<<"Gamma"<<i<<" = "<<gtag[i]<<endl;
	  print_dirac(base_gamma[i]);
	  cout<<endl;
	}

      //check of gamma6
      dirac_matr gamma6;
      dirac_prod(gamma6,base_gamma[1],base_gamma[5]);
      cout<<"Gamma6 = Gamma1Gamma5  = "<<gtag[6]<<" as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(base_gamma[6]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_dirac(gamma6);
      cout<<endl;

      //check of gamma7
      dirac_matr gamma7;
      dirac_prod(gamma7,base_gamma[2],base_gamma[5]);
      cout<<"Gamma7 = Gamma2Gamma5 = "<<gtag[7]<<" as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(base_gamma[7]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_dirac(gamma7);
      cout<<endl;

      //check of gamma8
      dirac_matr gamma8;
      dirac_prod(gamma8,base_gamma[3],base_gamma[5]);
      cout<<"Gamma8 = Gamma3Gamma5 = "<<gtag[8]<<" as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(base_gamma[8]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_dirac(gamma8);
      cout<<endl;

      //check of gamma9
      dirac_matr gamma9;
      dirac_prod(gamma9,base_gamma[4],base_gamma[5]);
      cout<<"Gamma9 = Gamma4Gamma5 = "<<gtag[9]<<" as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(base_gamma[9]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_dirac(gamma9);
      cout<<endl;
      
      //check of gamma10
      dirac_matr gamma10;
      dirac_prod(gamma10,base_gamma[4],base_gamma[1]);
      cout<<"Gamma10 = Gamma4Gamma1 = "<<gtag[10]<<" as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(base_gamma[10]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_dirac(gamma10);
      cout<<endl;

      //check of gamma11
      dirac_matr gamma11;
      dirac_prod(gamma11,base_gamma[4],base_gamma[2]);
      cout<<"Gamma11 = Gamma4Gamma2 = "<<gtag[11]<<" as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(base_gamma[11]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_dirac(gamma11);
      cout<<endl;

      //check of gamma12
      dirac_matr gamma12;
      dirac_prod(gamma12,base_gamma[4],base_gamma[3]);
      cout<<"Gamma12 = Gamma4Gamma3 = "<<gtag[12]<<" as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(base_gamma[12]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_dirac(gamma12);
      cout<<endl;

      //check of gamma13
      dirac_matr gamma13;
      dirac_prod(gamma13,base_gamma[2],base_gamma[3]);
      cout<<"Gamma13 = Gamma2Gamma3 = "<<gtag[13]<<" as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(base_gamma[13]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_dirac(gamma13);
      cout<<endl;

      //check of gamma14
      dirac_matr gamma14;
      dirac_prod(gamma14,base_gamma[3],base_gamma[1]);
      cout<<"Gamma14 = Gamma3Gamma1 = "<<gtag[14]<<" as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(base_gamma[14]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_dirac(gamma14);
      cout<<endl;

      //check of gamma15
      dirac_matr gamma15;
      dirac_prod(gamma15,base_gamma[1],base_gamma[2]);
      cout<<"Gamma15 = Gamma1Gamma2 = "<<gtag[15]<<" as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(base_gamma[15]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_dirac(gamma15);
      cout<<endl;


      complex one_over_rad2={1/sqrt(2),0};

      //check of Pplus
      dirac_matr gamma_Pplus;
      dirac_compl_prod(gamma_Pplus,base_gamma[5],I);
      dirac_summ(gamma_Pplus,gamma_Pplus,base_gamma[0]);
      dirac_compl_prod(gamma_Pplus,gamma_Pplus,one_over_rad2);
      cout<<"Pplus=(1+iGamma5)/sqrt(2) as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(Pplus);
      cout<<"(b) obtained by making explicite computations"<<endl;
      print_dirac(gamma_Pplus);
      cout<<endl;

      //check of Pminus
      dirac_matr gamma_Pminus;
      complex mI={0,-1};
      dirac_compl_prod(gamma_Pminus,base_gamma[5],mI);
      dirac_summ(gamma_Pminus,gamma_Pminus,base_gamma[0]);
      dirac_compl_prod(gamma_Pminus,gamma_Pminus,one_over_rad2);
      cout<<"Pminus=(1-iGamma5)/sqrt(2) as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_dirac(Pminus);
      cout<<"(b) obtained by making explicite computations"<<endl;
      print_dirac(gamma_Pminus);
      cout<<endl;

    }
  
  return 0;
}
