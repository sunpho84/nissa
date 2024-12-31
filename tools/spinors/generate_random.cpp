#include "nissa.h"

#include <stdlib.h>

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();
  
  if(nissa_nranks>1) CRASH("cannot run in parallel");
  
  if(narg<4) CRASH("use: %s L T seed file_out",arg[0]);

  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  int seed=atoi(arg[3]);

  //Init the MPI grid 
  init_grid(T,L);

  ///////////////////////////////////////////

  quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);
  
  start_loc_rnd_gen(seed);
  
  rnd_fill_unif_loc_vector((double*)conf,72,-1,1);
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      su3_unitarize_explicitly_inverting(conf[ivol][mu],conf[ivol][mu]);
  
  write_ildg_gauge_conf(arg[4],conf,64);
  
  nissa_free(conf);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
