#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/propagators/twisted_propagator_g2_corr.h"
#include "../src/stochastic/stochastic_twisted_propagator.h"
#include "../src/types/types_routines.h"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.h"

spinspin *prop_S0;
spinspin *d2_corr,*d2_corr_ave;
spin1field *phi,*eta;

//initialize the program
void init_test()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(8,4);
  
  //allocate propagators
  prop_S0=nissa_malloc("prop_S0",loc_vol,spinspin);
  d2_corr=nissa_malloc("d2_corr",loc_vol,spinspin);
  d2_corr_ave=nissa_malloc("d2_corr_ave",loc_vol,spinspin);
  //allocate stoch field propagator
  phi=nissa_malloc("phi",loc_vol,spin1field);
  eta=nissa_malloc("eta",loc_vol,spin1field);
}

//close the program
void close_test()
{
  nissa_free(prop_S0);
  nissa_free(d2_corr);
  nissa_free(d2_corr_ave);
  nissa_free(eta);
  nissa_free(phi);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  int nsource=1;
  
  init_test();
  
  //quark
  double quark_theta[4]={1,0,0,0};
  double kappa=1.0/8;
  double mass=0;
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double gluon_theta[4]={0,1,0,0};
  double alpha=0.3;
  gluon_info gl=create_Wilson_gluon_info(alpha,gluon_theta);
  
  start_loc_rnd_gen(100);
  
  /////////////////////////////// S0 propagator=id //////////////////////////
  
  nissa_loc_vol_loop(ivol)
    {
      spinspin_put_to_id(prop_S0[ivol]);
      spinspin_put_to_zero(d2_corr_ave[ivol]);
    }
  
  /////////////////////////////// D2 correction //////////////////////////
  
  for(int isource=0;isource<nsource;isource++)
    {
      generate_stochastic_source_and_tlSym_gluon_propagator(phi,eta,gl);
      generate_stochastic_A_B_dag_twisted_propagator_source_by_inv(d2_corr,prop_S0,qu,phi,eta,gl);
      
      nissa_loc_vol_loop(ivol)
	for(int id=0;id<4;id++)
	  spin_summ(d2_corr_ave[ivol][id],d2_corr_ave[ivol][id],d2_corr[ivol][id]);
    }      

  //////////////////////////////////// output ////////////////////////////
  
  nissa_loc_vol_loop(ivol)
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	complex_prodassign_double(d2_corr[ivol][id1][id2],nsource);

  pass_spinspin_from_mom_to_x_space(d2_corr_ave,d2_corr_ave,quark_theta);
  
  print_spinspin(d2_corr_ave[32]);
  
  spinspin temp;
  mom_space_twisted_propagator_g2_d2_corr_of_imom(temp,qu,gl,32);
  print_spinspin(temp);

  close_test();
  
  return 0;
}
