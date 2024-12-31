#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/types/types.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/propagators/Wilson_gluon_propagator.hpp"
#include "../src/operators/Wilson_gluon_Klein_Gordon_operator.hpp"
#include "../src/propagators/tlSym_gluon_propagator.hpp"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.hpp"

spin1prop *prop_wi;
spin1prop *prop_wi_stoch;
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
  prop_wi=nissa_malloc("prop_wi",loc_vol,spin1prop);
  prop_wi_stoch=nissa_malloc("prop_wi_stoch",loc_vol,spin1prop);
  
  //allocate stoch field propagator
  phi=nissa_malloc("phi",loc_vol,spin1field);
  eta=nissa_malloc("eta",loc_vol,spin1field);
}

//close the program
void close_test()
{
  nissa_free(prop_wi);
  nissa_free(prop_wi_stoch);
  nissa_free(eta);
  nissa_free(phi);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_test(narg,arg);
  
  //anti-periodic boundary condition in one space direction
  double theta[4]={0,1,0,0};
  
  //covariant gauge fixing constant
  double alpha=0.3;
  
  //gluon information
  gluon_info gl=create_tlSym_gluon_info(alpha,theta,0);
  
  //compute the point to be printed
  coords ix={1,2,1,2};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  /////////////////////////////// analytic propagator //////////////////////////
  
  compute_x_space_tlSym_gluon_propagator_by_fft(prop_wi,gl);
  
  /////////////////////////////// stochastic propagator //////////////////////////
  
  memset(prop_wi_stoch,0,sizeof(spin1prop)*loc_vol);
  
  int nsource=1<<18;
  int how=1;
  for(int isource=0;isource<nsource;isource++)
    {
      //generate eta and phi
      generate_stochastic_source_and_tlSym_gluon_propagator(phi,eta,gl);
      
      //broadcast eta[0]
      spin1field eta_0;
      if(rank==0) memcpy(eta_0,eta[0],sizeof(spin1field));
      MPI_Bcast(eta_0,sizeof(spin1field),MPI_CHAR,0,MPI_COMM_WORLD);
      
      //reconstruct prop[ivol][mu][nu]+=phi[ivol][mu]*eta[0][nu]
      NISSA_LOC_VOL_LOOP(ivol)
	for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++)
	    complex_summ_the_conj2_prod(prop_wi_stoch[ivol][mu][nu],phi[ivol][mu],eta_0[nu]);
  
      if((isource+1)%how==0)
	{
	  //take the squared norm of the difference between the two computed propagators  
	  double loc_d=0;
	  NISSA_LOC_VOL_LOOP(ivol)
	    for(int id=0;id<32;id++)
	      {
		double a=((double*)prop_wi[ivol])[id];
		double b=((double*)prop_wi_stoch[ivol])[id]/(isource+1);
		double t=a-b;
		loc_d+=t*t;
	      }
	  how*=2;
	  double glb_d=sqrt(glb_reduce_double(loc_d)/glb_vol);
	  MASTER_PRINTF("Average norm2 difference between analytic and stochastic gluon propagators with %d sources: %lg\n",isource+1,glb_d);
	}
    }
  
  //normalize
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      for(int nu=0;nu<4;nu++)
	complex_prodassign_double(prop_wi_stoch[ivol][mu][nu],1.0/nsource);
  
  if(rx==rank)
    {
      spinspin_print(prop_wi_stoch[lx]);
      spinspin_print(prop_wi[lx]);
    }
  
  close_test();
  
  return 0;
}
