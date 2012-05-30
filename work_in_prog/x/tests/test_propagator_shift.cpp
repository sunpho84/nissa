#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/types/types_routines.h"
#include "../src/routines/fourier.h"
#include "../src/routines/shift.h"

spinspin  *q_prop,*q_prop_sh;

//initialize the program
void init_test()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(8,4);
  
  //allocate propagators
  q_prop_sh=nissa_malloc("q_prop_sh",loc_vol,spinspin);
  q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
}

//close the program
void close_test()
{
  nissa_free(q_prop_sh);
  nissa_free(q_prop);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_test();
  
  //quark
  double quark_theta[4]={0.3,0.1,0.3,0.34};
  double kappa=0.177;
  double mass=0.5;
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  /////////////////////////////////////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  shift_spinspin_source_up(q_prop_sh,q_prop,quark_theta,0);
  shift_spinspin_source_up(q_prop_sh,q_prop_sh,quark_theta,1);
  shift_spinspin_source_up(q_prop_sh,q_prop_sh,quark_theta,2);
  shift_spinspin_source_up(q_prop_sh,q_prop_sh,quark_theta,3);
  shift_spinspin_source_up(q_prop_sh,q_prop_sh,quark_theta,2);
  shift_spinspin_sink_up(q_prop_sh,q_prop_sh,quark_theta,3);
  shift_spinspin_sink_up(q_prop_sh,q_prop_sh,quark_theta,3);
  
  /////////////////////////////////////////////////////////
  
  coords r={1,1,2,-1};
  shift_spinspin_source(q_prop,q_prop,quark_theta,r);

  double d=0;
  nissa_loc_vol_loop(ivol)
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	for(int ri=0;ri<2;ri++)
	  {
	    double t=q_prop[ivol][id1][id2][ri]-q_prop_sh[ivol][id1][id2][ri];
	    d+=t*t;
	  }
  d=sqrt(glb_reduce_double(d)/glb_vol);
  master_printf("Difference between shift_up and shift: %lg\n",d);
  
  //////////////////////////////////// output ////////////////////////////
  
  close_test();
  
  return 0;
}
