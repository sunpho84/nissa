#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/diagrams/propagator_self_energy.h"
#include "../src/types/types_routines.h"
#include "../src/routines/shift.h"
#include "../src/routines/read_and_write.h"
#include "../src/routines/correlations.h"
#include "../src/stochastic/stochastic_twisted_propagator.h"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.h"

spinspin *prop,*self_prop;
corr16 *corr;

int L=4;

//initialize the program
void init_calc()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(2*L,L);
  
  //allocate propagators
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
  init_calc();
  
  //kappa and mass
  double kappa=1.0/8;
  double mass=0.00;
  double quark_theta[4]={0.3,0.3,0.3,0.3};
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double alpha=0;
  double gluon_theta[4]={0,0,0,0};
  gluon_info gl=create_tlSym_gluon_info(alpha,gluon_theta);
  
  ////////////////////////////////////// propagators computed analytically ////////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(prop,qu);
  
  //////////////////////////////// compute correlation and write them on disk ////////////////////////

  compute_all_2pts_qdagq_correlations(corr,prop,prop);
  
  //compute a particular point in a different way
  nissa_loc_vol_loop(P)
    {
      spinspin t;
      compute_x_space_propagator_to_sink_from_source(t,prop,qu.bc,glb_coord_of_loclx[0],glb_coord_of_loclx[P]);
      spinspin pr;
      unsafe_spinspin_spinspin_prod(pr,prop[P],t);
      complex c;
      trace_spinspin(c,pr);
      master_printf("%d %lg %lg\n",P,corr[P][0][0],c[0]);
    }
  
  write_corr16("tree_level_corr",corr,64);
  
  close_calc();
  
  return 0;
}
