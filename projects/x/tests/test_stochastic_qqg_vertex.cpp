#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.hpp"
#include "../src/types/types.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/vertex/x_space_stochastic_qqg_vertex.hpp"

spinspin *S0_prop,*Sphi_prop;
spin1field *phi,*eta;

//initialize the program
void init_test(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //init the grid
  init_grid(8,4);
  
  //start loc rnd gen
  start_loc_rnd_gen(1);
  
  //allocate propagators
  S0_prop=nissa_malloc("S0_prop",loc_vol,spinspin);
  Sphi_prop=nissa_malloc("Sphi_prop",loc_vol,spinspin);
  
  //allocate stoch field propagator
  phi=nissa_malloc("phi",loc_vol,spin1field);
  eta=nissa_malloc("eta",loc_vol,spin1field);
}

//close the program
void close_test()
{
  nissa_free(S0_prop);
  nissa_free(Sphi_prop);
  nissa_free(eta);
  nissa_free(phi);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_test(narg,arg);
  
  //non-trivial boundary condition in one space direction
  double theta_gluon[4]={0,1,0,0};
  //anti-periodic boundary condition in time
  double theta_quark[4]={0,1,0,0};
  
  //covariant gauge fixing constant
  double alpha=0.3;
  //kappa
  double kappa=1.0/8;
  //mass
  double mass=0;
  
  //gluon information
  gluon_info gl=create_tlSym_gluon_info(alpha,theta_gluon);
  quark_info qu=create_twisted_quark_info(kappa,mass,theta_quark);
  
  //compute the point to be printed
  coords ix={1,2,1,2};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  //////////////////////////////// S0 quark propagator //////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(S0_prop,qu);  
  
  //////////////////////////// stochastic gluon propagator //////////////////////
  
  generate_stochastic_source_and_tlSym_gluon_propagator(phi,eta,gl);

  //////////////////////////////////// vertex ///////////////////////////////////
  
  stochastic_x_space_qqg_vertex(Sphi_prop,S0_prop,qu,phi,gl);
  
  if(rx==rank)
    {
      spinspin_print(S0_prop[lx]);
      spinspin_print(Sphi_prop[lx]);
    }
  
  close_test();
  
  return 0;
}
