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
      for(int i=0;i<=5;i++)
	{
	  cout<<"Gamma"<<i<<endl;
	  print_gamma(base_gamma[i]);
	  cout<<endl;
	}

      gamma gamma6;
      gamma_prod(gamma6,base_gamma[1],base_gamma[5]);
      cout<<"Gamma6=Gamma1Gamma5 as:"<<endl;
      cout<<"(a) defined by Silvano"<<endl;
      print_gamma(base_gamma[6]);
      cout<<"(b) obtained by making explicite products"<<endl;
      print_gamma(gamma6);
      cout<<endl;
    }
  
  return 0;
}
