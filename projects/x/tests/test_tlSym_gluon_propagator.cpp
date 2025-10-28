#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/types/types.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/propagators/tlSym_gluon_propagator.hpp"


spin1prop *prop_fft;
spin1prop *prop_inv;

//initialize the program
void init_test(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
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
  init_test(narg,arg);
  
  //anti-periodic boundary condition in one space direction
  double theta[4]={0,0,0,0};
  
  //covariant gauge fixing constant
  double alpha=0.3;
  
  //gluon info
  gluon_info gl=create_tlSym_gluon_info(alpha,theta);
  
  //compute the point to be printed
  coords ix={3,3,2,1};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  /////////////////////////////// propagator and pion computed analytically //////////////////////////
  
  compute_x_space_tlSym_gluon_propagator_by_fft(prop_fft,gl);
  pass_spin1prop_from_x_to_mom_space(prop_fft,prop_fft,gl.bc);

  /////////////////////////////// propagator and pion computed numerically //////////////////////////
  
  //compute_x_space_tlSym_gluon_propagator_by_inv(prop_inv,gl);
  
  /////////////////////////////////////////// output ////////////////////////////////////////////////
  
  if(rank==rx)
    {
      printf("\n\nComparing the propagator on site of coordinates: (%d,%d,%d,%d), rank: %d\n",ix[0],ix[1],ix[2],ix[3],rx);
      spinspin_print(prop_fft[lx]);
      //printf("\n");
      //spinspin_print(prop_inv[lx]);
    }

  //take the squared norm of the differnce between the two computed propagators
  
  /*
  double loc_d=0;
  NISSA_LOC_VOL_LOOP(ivol)
    for(int id=0;id<32;id++)
      {
	double t=((double*)prop_inv[ivol])[id]-((double*)prop_fft[ivol])[id];
	loc_d+=t*t;
      }
  double glb_d=sqrt(glb_reduce_double(loc_d)/glb_vol);
  MASTER_PRINTF("\n\nAverage norm2 difference between fft and inv computed tlSym propagators: %lg\n\n",glb_d);
  */
  
  close_test();
  
  return 0;
}
