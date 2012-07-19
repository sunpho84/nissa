#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/types/types_routines.h"
#include "../src/routines/read_and_write.h"
#include "../src/routines/correlations.h"

spinspin *prop;
corr16 *corr;

//initialize the program
void init_calc(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();
  
  if(narg<3) crash("use %s T L",arg[0]);
  int T=atoi(arg[1]);
  int L=atoi(arg[2]);

  //init the grid
  init_grid(T,L);
  
  //allocatepropagators
  prop=nissa_malloc("prop",loc_vol,spinspin);
  corr=nissa_malloc("corr",loc_vol,corr16);
}

//close the program
void close_calc()
{
  nissa_free(prop);
  nissa_free(corr);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_calc(narg,arg);
  
  //anti-periodic boundary condition in time
  double theta[4]={1,0,0,0};
  
  //kappa and mass
  double kappa=1.0/8;
  double mass=0.00;
  quark_info qu=create_twisted_quark_info(kappa,mass,theta);
  
  compute_x_space_twisted_propagator_by_fft(prop,qu);  
  compute_all_2pts_qdagq_correlations(corr,prop,prop);
  
  write_corr16("tree",corr,64);
  
  close_calc();
  
  return 0;
}
