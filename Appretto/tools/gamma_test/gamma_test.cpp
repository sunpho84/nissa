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
      print_gamma(base_gamma[0]);
      cout<<endl;
    }

  return 0;
}
