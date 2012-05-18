#include <stdio.h>

#include "nissa.h"

#include "../src/types.h"
#include "../src/read_and_write.h"

corr16 *unav_corr;

//initialize the program
void init_calc()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(48,24);
  
  //allocate correlator
  unav_corr=nissa_malloc("corr",loc_vol,corr16);
}

//close the program
void close_calc()
{
  nissa_free(unav_corr);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_calc();
  if(rank_tot>1) crash("only available in scalar");
  
  read_corr16(unav_corr,"corr");
  
  FILE *f00=open_file("corr00_tau32-0_L24_T48_free.dat","w");
  
  coords x;
  for(x[0]=0;x[0]<glb_size[1]/3;x[0]++)
    for(x[1]=0;x[1]<glb_size[1]/3;x[1]++)
      for(x[2]=0;x[2]<glb_size[2]/3;x[2]++)
	for(x[3]=0;x[3]<glb_size[3]/3;x[3]++)
	  {
	    int i=glblx_of_coord(x);
	    
	    fprintf(f00,"%d %d %d %d %16.16le\n",x[1],x[2],x[3],x[0],unav_corr[i][5][0]);
	  }
  
  fclose(f00);
  
  close_calc();
  
  return 0;
}
