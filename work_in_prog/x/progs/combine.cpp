#include <stdio.h>
#include <stdlib.h>

#include "nissa.hpp"
using namespace std;

#include "../src/types/types.hpp"
#include "../src/routines/read_and_write.hpp"

corr16 *loaded_corr,*combined_corr;
int N;

//initialize the program
void init_calc(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<6||narg%2) crash("use %s T L outfile file1 weight1 file2...",arg[0]);
  
  int T=atoi(arg[1]);
  int L=atoi(arg[2]);
  
  //init the grid
  init_grid(T,L);
  
  //allocate correlator
  loaded_corr=nissa_malloc("loaded_corr",loc_vol,corr16);
  combined_corr=nissa_malloc("combined_corr",loc_vol,corr16);
}

//close the program
void close_calc()
{
  nissa_free(combined_corr);
  nissa_free(loaded_corr);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_calc(narg,arg);
  
  memset(combined_corr,0,sizeof(loc_vol)*sizeof(corr16));
  
  int N=narg-4;
  for(int i=4;i<narg;i+=2)
    {
      double w;
      sscanf(arg[i+1],"%lg",&w);
      master_printf("Loading file %s and combining with weight %lg\n",arg[i],w);
      
      read_corr16(loaded_corr,arg[i]);
      
      NISSA_LOC_VOL_LOOP(ivol)
	for(int ig=0;ig<16;ig++)
	  complex_summ_the_prod_double(combined_corr[ivol][ig],loaded_corr[ivol][ig],w);
    }
  
  write_corr16(arg[3],combined_corr,64);
  
  close_calc();
  
  return 0;
}
