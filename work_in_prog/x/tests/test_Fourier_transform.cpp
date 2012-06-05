#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/stochastic/stochastic_source.h"
#include "../src/routines/fourier.h"
#include "../src/types/types_routines.h"

spin1field *eta;
spin1field *etab;

//initialize the program
void init_test()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(8,4);
  
  start_loc_rnd_gen(100);
  
  //allocatepropagators
  eta=nissa_malloc("eta",loc_vol,spin1field);
  etab=nissa_malloc("etab",loc_vol,spin1field);
}

//close the program
void close_test()
{
  nissa_free(eta);
  nissa_free(etab);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_test();
  
  //covariant gauge fixing constant
  double alpha=0.3;
  double theta_gluon[4]={0.1,0.7,0.3,0.8};
  gluon_info gl=create_tlSym_gluon_info(alpha,theta_gluon);
  
  generate_stochastic_source_eta(eta);
  pass_spin1field_from_x_to_mom_space(etab,eta,gl.bc);
  pass_spin1field_from_mom_to_x_space(etab,etab,gl.bc);
  
  //compute the point to be printed
  coords ix={1,2,1,2};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);

  if(rank==rx)
    {
      printf("\n\nComparing the vector on site of coordinates: (%d,%d,%d,%d), rank: %d\n",ix[0],ix[1],ix[2],ix[3],rx);
      print_spin(eta[lx]);
      printf("\n");
      print_spin(etab[lx]);
    }
  
  close_test();
    
  return 0;
}
