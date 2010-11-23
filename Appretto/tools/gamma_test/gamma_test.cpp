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
	  cout<<"Gamma"<<i<<endl;
	  print_gamma(base_gamma[i]);
	  cout<<endl;
	}

      //check of gamma6
      gamma gamma6;
      gamma_prod(gamma6,base_gamma[1],base_gamma[5]);
      cout<<"Gamma6=Gamma1Gamma5 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[6]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma6);
      cout<<endl;

      //check of gamma7
      gamma gamma7;
      gamma_prod(gamma7,base_gamma[2],base_gamma[5]);
      cout<<"Gamma7=Gamma2Gamma5 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[7]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma7);
      cout<<endl;

      //check of gamma8
      gamma gamma8;
      gamma_prod(gamma8,base_gamma[3],base_gamma[5]);
      cout<<"Gamma8=Gamma3Gamma5 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[8]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma8);
      cout<<endl;

      //check of gamma9
      gamma gamma9;
      gamma_prod(gamma9,base_gamma[4],base_gamma[5]);
      cout<<"Gamma9=Gamma4Gamma5 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[9]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma9);
      cout<<endl;

      //check of gamma10
      gamma gamma10;
      gamma_prod(gamma10,base_gamma[4],base_gamma[1]);
      cout<<"Gamma10=Gamma4Gamma1 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[10]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma10);
      cout<<endl;

      //check of gamma11
      gamma gamma11;
      gamma_prod(gamma11,base_gamma[4],base_gamma[2]);
      cout<<"Gamma11=Gamma4Gamma2 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[11]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma11);
      cout<<endl;

      //check of gamma12
      gamma gamma12;
      gamma_prod(gamma12,base_gamma[4],base_gamma[3]);
      cout<<"Gamma12=Gamma4Gamma3 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[12]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma12);
      cout<<endl;

      //check of gamma13
      gamma gamma13;
      gamma_prod(gamma13,base_gamma[2],base_gamma[3]);
      cout<<"Gamma13=Gamma2Gamma3 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[13]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma13);
      cout<<endl;

      //check of gamma14
      gamma gamma14;
      gamma_prod(gamma14,base_gamma[3],base_gamma[1]);
      cout<<"Gamma14=Gamma3Gamma1 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[14]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma14);
      cout<<endl;

      //check of gamma15
      gamma gamma15;
      gamma_prod(gamma15,base_gamma[1],base_gamma[2]);
      cout<<"Gamma15=Gamma1Gamma2 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[15]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma15);
      cout<<endl;


    }
  
  return 0;
}
