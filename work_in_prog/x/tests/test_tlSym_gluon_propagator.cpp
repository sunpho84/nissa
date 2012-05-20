#include <math.h>

#include "nissa.h"

#include "../src/types/types.h"
#include "../src/types/types_routines.h"
#include "../src/propagators/tlSym_gluon_propagator.h"


spin1prop *prop_fft;
spin1prop *prop_inv;

//initialize the program
void init_test()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(8,4);
  
  //allocatepropagators
  prop_fft=nissa_malloc("prop_fft",loc_vol,spin1prop);
  prop_inv=nissa_malloc("prop_inv",loc_vol,spin1prop);
}

//close the program
void close_test()
{
  nissa_free(prop_fft);
  nissa_free(prop_inv);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_test();
  
  //anti-periodic boundary condition in one space direction
  double theta[4]={0,1,0,0};
  
  //covariant gauge fixing constant
  double alpha=0.3;
  
  //gluon info
  gluon_info gl=create_tlSym_gluon_info(alpha,theta);
  
  //compute the point to be printed
  coords ix={1,2,1,2};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  /////////////////////////////// propagator and pion computed analytically //////////////////////////
  
  compute_x_space_tlSym_gluon_propagator_by_fft(prop_fft,gl);
  
  /////////////////////////////// propagator and pion computed numerically //////////////////////////
  
  //compute_x_space_tlSym_gluon_propagator_by_inv(prop_inv,gl);
  
  /////////////////////////////////////////// output ////////////////////////////////////////////////
  
  if(rank==rx)
    {
      printf("\n\nComparing the propagator on site of coordinates: (%d,%d,%d,%d), rank: %d\n",ix[0],ix[1],ix[2],ix[3],rx);
      print_spinspin(prop_fft[lx]);
      //printf("\n");
      //print_spinspin(prop_inv[lx]);
    }

  //take the squared norm of the differnce between the two computed propagators
  
  /*
  double loc_d=0;
  nissa_loc_vol_loop(ivol)
    for(int id=0;id<32;id++)
      {
	double t=((double*)prop_inv[ivol])[id]-((double*)prop_fft[ivol])[id];
	loc_d+=t*t;
      }
  double glb_d=sqrt(glb_reduce_double(loc_d)/glb_vol);
  master_printf("\n\nAverage norm2 difference between fft and inv computed tlSym propagators: %lg\n\n",glb_d);
  */
  
  close_test();
  
  return 0;
}
