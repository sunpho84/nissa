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
  init_grid(4,4);
  
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
  int nsource=64*64*64*64;
  
  init_test();
  
  //quark
  double quark_theta[4]={0,0,0,0};
  double kappa=1.0/8;
  double mass=0;
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double gluon_theta[4]={0,0,0,0};
  double alpha=0.3;
  gluon_info gl=create_Wilson_gluon_info(alpha,gluon_theta);
  
  start_loc_rnd_gen(100);
  
  /////////////////////////////// S0 propagator=id //////////////////////////
  
  memset(d2_corr_ave,0,sizeof(spinspin)*loc_vol);
  memset(prop_S0,0,sizeof(spinspin)*loc_vol);
  spinspin_put_to_id(prop_S0[0]);
  
  /////////////////////////////// D2 correction //////////////////////////
  
  int log2=0;
  for(int isource=0;isource<nsource;isource++)
    {
      generate_stochastic_source_and_tlSym_gluon_propagator(phi,eta,gl);
      
      generate_stochastic_A_B_dag_twisted_propagator_source(d2_corr,prop_S0,qu,phi,eta,gl);
      
      nissa_loc_vol_loop(ivol)
	spinspin_summ(d2_corr_ave[ivol],d2_corr_ave[ivol],d2_corr[ivol]);
      
      //////////////////////////////////// output ////////////////////////////
      
      if(isource==1<<log2)
	{
	  log2++;
	  nissa_loc_vol_loop(ivol)
	    for(int id1=0;id1<4;id1++)
	      for(int id2=0;id2<4;id2++)
		complex_prod_double(d2_corr[ivol][id1][id2],d2_corr_ave[ivol][id1][id2],-4.0/isource);
	  
	  pass_spinspin_from_x_to_mom_space(d2_corr,d2_corr,quark_theta);
	  
	  coords ix={0,0,0,1};
	  int lx,rx;
	  get_loclx_and_rank_of_coord(&lx,&rx,ix);
	  master_printf("%d\n",isource);
	  print_spinspin(d2_corr[lx]);
	  master_printf("\n");
	}
    }
  
  close_test();
  
  return 0;
}
