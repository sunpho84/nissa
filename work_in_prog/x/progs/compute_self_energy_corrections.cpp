#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/diagrams/propagator_self_energy.h"
#include "../src/types/types_routines.h"
#include "../src/routines/read_and_write.h"
#include "../src/routines/correlations.h"

spinspin *prop,*self_prop;
corr16 *corr;

//initialize the program
void init_calc()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(48,24);
  
  //allocate propagators
  prop=nissa_malloc("prop",loc_vol,spinspin);
  self_prop=nissa_malloc("self_prop",loc_vol,spinspin);
  corr=nissa_malloc("corr",loc_vol,corr16);
}

//close the program
void close_calc()
{
  nissa_free(prop);
  nissa_free(self_prop);
  nissa_free(corr);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_calc();
  
  //kappa and mass
  double kappa=1.0/8;
  double mass=0.00;
  double quark_theta[4]={1,0,0,0};
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double alpha=0;
  double gluon_theta[4]={0,0,0,0};
  gluon_info gl=create_tlSym_gluon_info(alpha,gluon_theta);
  
  /////////////////////////////// propagator and pion computed analytically //////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(prop,qu);
  compute_self_energy_twisted_propagator_in_x_space(self_prop,qu,gl);
  
  //////////////////////////////// compute correlation and write them on disk ////////////////////////

  compute_all_2pts_qdagq_correlations(corr,prop,self_prop);
  
  write_corr16("self_energy_corr",corr,64);
  
  close_calc();
  
  return 0;
}
