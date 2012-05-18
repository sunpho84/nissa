#include <math.h>

#include "nissa.h"

#include "../src/twisted_propagator.h"
#include "../src/types_routines.h"

spinspin *prop_fft;
spinspin *prop_inv;

//initialize the program
void init_test()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(8,4);
  
  //allocatepropagators
  prop_fft=nissa_malloc("prop_fft",loc_vol,spinspin);
  prop_inv=nissa_malloc("prop_fft",loc_vol,spinspin);
}

//close the program
void close_test()
{
  nissa_free(prop_fft);
  nissa_free(prop_inv);
  
  close_nissa();
}

void compute_pion_correlator(complex *corr,spinspin *p)
{
  //stupid propagators
  su3spinspin *stp=nissa_malloc("stp",loc_vol,su3spinspin);
  
  //stupidly promote the propagator to su3
  memset(stp,0,sizeof(su3spinspin)*loc_vol);
  nissa_loc_vol_loop(ivol)
    for(int ic1=0;ic1<1;ic1++)
      for(int ic2=0;ic2<1;ic2++)
	memcpy(stp[ivol][ic1][ic2],p[ivol],sizeof(spinspin));
  
  //gamma5=gamma0 because of revert
  trace_g_ccss_dag_g_ccss(corr,&(base_gamma[0]),stp,&(base_gamma[0]),stp,1); 
  
  nissa_free(stp);
}

int main(int narg,char **arg)
{
  init_test();
  
  //anti-periodic boundary condition in time
  double theta[4]={1,0,0,0};
  
  //kappa and mass
  double kappa=0.177;
  double mass=0.01;
  quark_info qu=create_twisted_quark_info(kappa,mass,theta);
  
  //pion correlator
  complex pion_corr_fft[glb_size[0]];
  complex pion_corr_inv[glb_size[0]];
  int five[1]={5};
  
  //compute the point to be printed
  coords ix={1,2,1,2};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  /////////////////////////////// propagator and pion computed analytically //////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(prop_fft,qu);  
  compute_pion_correlator(pion_corr_fft,prop_fft);

  /////////////////////////////// propagator and pion computed numerically //////////////////////////
  
  compute_x_space_twisted_propagator_by_inverting(prop_inv,qu);
  compute_pion_correlator(pion_corr_inv,prop_inv);

  /////////////////////////////////////////// output ////////////////////////////////////////////////
  
  master_printf("\nComparing the pion correlator\n");
  
  print_contractions_to_file(stdout,1,five,five,pion_corr_fft,0,"",1.0);
  print_contractions_to_file(stdout,1,five,five,pion_corr_inv,0,"",1.0);

  if(rank==rx)
    {
      printf("\n\nComparing the propagator on site of coordinates: (%d,%d,%d,%d), rank: %d\n",ix[0],ix[1],ix[2],ix[3],rx);
      print_spinspin(prop_fft[lx]);
      printf("\n");
      print_spinspin(prop_inv[lx]);
    }

  //take the squared norm of the differnce between the two computed propagators
  double loc_d=0;
  nissa_loc_vol_loop(ivol)
    for(int id=0;id<32;id++)
      {
	double t=((double*)prop_fft[ivol])[id]-((double*)prop_inv[ivol])[id];
	loc_d+=t*t;
      }
  double glb_d=sqrt(glb_reduce_double(loc_d)/glb_vol);
  master_printf("\n\nAverage norm2 difference between fft and inv computed propagators: %lg\n\n",glb_d);
      
  close_test();
  
  return 0;
}
