#include <math.h>

#include "nissa.h"

#include "../src/twisted_propagator.h"
#include "../src/types_routines.h"
#include "../src/read_and_write.h"

spinspin *prop;
corr16 *corr;

//initialize the program
void init_calc()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(48,24);
  
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
  init_calc();
  
  //anti-periodic boundary condition in time
  double theta[4]={1,0,0,0};
  
  //kappa and mass
  double kappa=1.0/8;
  double mass=0.00;
  quark_info qu=create_twisted_quark_info(kappa,mass,theta);
  
  /////////////////////////////// propagator and pion computed analytically //////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(prop,qu);  
  
  dirac_matr t1[16],t2[16];
  for(int igamma=0;igamma<16;igamma++)
    {
      dirac_prod(&(t1[igamma]),&(base_gamma[igamma]),&(base_gamma[5]));
      dirac_prod(&(t2[igamma]),&(base_gamma[5]),&(base_gamma[igamma]));
    }
      
  nissa_loc_vol_loop(ivol)
    for(int igamma=0;igamma<16;igamma++)
      {
	site_trace_g_sdag_g_s(corr[ivol][igamma],&(t1[igamma]),prop[ivol],&(t1[igamma]),prop[ivol]); 
	complex_prodassign_double(corr[ivol][igamma],3);
      }
  
  write_corr16("corr",corr,64);
  
  close_calc();
  
  return 0;
}
