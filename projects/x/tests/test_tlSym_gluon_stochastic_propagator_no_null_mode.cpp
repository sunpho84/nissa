#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/types/types.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/propagators/Wilson_gluon_propagator.hpp"
#include "../src/operators/Wilson_gluon_Klein_Gordon_operator.hpp"
#include "../src/propagators/tlSym_gluon_propagator.hpp"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.hpp"
#include "../src/inverters/cg_Wilson_gluon_operator.hpp"

spin1field *eta,*phi_fft,*phi_inv;

//initialize the program
void init_test(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //init the grid
  init_grid(24,32);
  
  //start loc rnd gen
  start_loc_rnd_gen(1);
  
  //allocate stoch field propagator
  eta=nissa_malloc("eta",loc_vol,spin1field);
  phi_fft=nissa_malloc("phi_fft",loc_vol,spin1field);
  phi_inv=nissa_malloc("phi_inv",loc_vol+bord_vol,spin1field);
}

//close the program
void close_test()
{
  nissa_free(eta);
  nissa_free(phi_fft);
  nissa_free(phi_inv);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_test(narg,arg);
  
  //anti-periodic boundary condition in one space direction
  double theta[4]={0,0,0,0};
  
  //covariant gauge fixing constant
  double alpha=0.3;
  
  //gluon information
  gluon_info gl=create_tlSym_gluon_info(alpha,theta,0);
  
  //compute the point to be printed
  coords ix={1,2,1,2};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  /////////////////////////////// stochastic propagator //////////////////////////
  
  //generate eta and phi
  generate_stochastic_source_and_tlSym_gluon_propagator(phi_fft,eta,gl);
  
  //compute the value of the null mode
  spin1field loc_ave;
  spin_put_to_zero(loc_ave);
  NISSA_LOC_VOL_LOOP(ivol)
    spin_summassign(loc_ave,eta[ivol]);
  spin1field glb_ave;
  MPI_Allreduce(loc_ave,glb_ave,sizeof(spin1field)/sizeof(double),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  spin_prodassign_double(glb_ave,1.0/glb_vol);
  
  //subtract the null mode
  NISSA_LOC_VOL_LOOP(ivol)
    spin_subtassign(eta[ivol],glb_ave);
  set_borders_invalid(eta);
  
  //invert
  inv_Wilson_gluon_Klein_Gordon_operator(phi_inv,NULL,gl,1000000,5,1.e-26,eta);
  
  //compute value of the null mode of the invertion
  spin_put_to_zero(loc_ave);
  NISSA_LOC_VOL_LOOP(ivol)
    spin_summassign(loc_ave,eta[ivol]);
  MPI_Allreduce(loc_ave,glb_ave,sizeof(spin1field)/sizeof(double),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MASTER_PRINTF("null mode of the out: %lg\n",glb_ave[0][0]);
  
  
  //take the squared norm of the difference between the two computed propagators  
  double loc_d=0;
  NISSA_LOC_VOL_LOOP(ivol)
    for(int id=0;id<32;id++)
      {
	double a=((double*)phi_fft[ivol])[id];
	double b=((double*)phi_inv[ivol])[id];
	double t=a-b;
	loc_d+=t*t;
      }
  double glb_d=sqrt(glb_reduce_double(loc_d)/glb_vol);
  MASTER_PRINTF("Difference between fft with explicit null mode and inversion: %lg\n",glb_d);
    
  close_test();
  
  return 0;
}
