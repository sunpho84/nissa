#include <math.h>
#include <time.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/diagrams/propagator_self_energy.h"
#include "../src/types/types_routines.h"
#include "../src/routines/read_and_write.h"
#include "../src/routines/correlations.h"
#include "../src/stochastic/stochastic_twisted_propagator.h"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.h"

spin1field *phi,*eta;
spinspin *prop,*self_prop;//,*id;
corr16 *corr,*summ_corr,*temp_corr;

int L=4;
int T=8;

//initialize the program
void init_calc()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(T,L);
  
  //allocate propagators
  //id=nissa_malloc("id",loc_vol,spinspin);
  prop=nissa_malloc("prop",loc_vol,spinspin);
  self_prop=nissa_malloc("self_prop",loc_vol,spinspin);
  summ_corr=nissa_malloc("summ_corr",loc_vol,corr16);
  temp_corr=nissa_malloc("temp_corr",loc_vol,corr16);
  corr=nissa_malloc("corr",loc_vol,corr16);
  phi=nissa_malloc("phi",loc_vol+bord_vol,spin1field);
  eta=nissa_malloc("eta",loc_vol+bord_vol,spin1field);
}

//close the program
void close_calc()
{
  //nissa_free(id);
  nissa_free(prop);
  nissa_free(self_prop);
  nissa_free(corr);
  nissa_free(temp_corr);
  nissa_free(summ_corr);
  nissa_free(eta);
  nissa_free(phi);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_calc();
  
  //kappa and mass
  double kappa=1.0/8;
  double mass=0.00;
  double quark_theta[4]={1,0,0,0};
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double alpha=0;
  double gluon_theta[4]={0,0,0,0};
  gluon_info gl=create_tlSym_gluon_info(alpha,gluon_theta);
  
  /////////////////////////////////// propagators computed analytically //////////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(prop,qu);
  compute_self_energy_twisted_propagator_in_x_space(self_prop,qu,gl);
  
  ////////////////////////////// compute correlation and write them on disk //////////////////////////

  compute_all_2pts_qdagq_correlations(corr,prop,self_prop);
  write_corr16("self_energy_corr",corr,64);
  
  master_printf("%lg\n",corr[9][0][0]);
    
  /////////////////////////////////////////////// tough way /////////////////////////////////////////
  
  compute_self_energy_twisted_propagator_in_x_space_tough_way(self_prop,qu,gl);

  compute_all_2pts_qdagq_correlations(corr,prop,self_prop);
  write_corr16("self_energy_corr_tough",corr,64);
  
  master_printf("%lg\n",corr[9][0][0]);
  
  ////////////////////////////////// propagators computed stochastically ////////////////////////////
  
  //memset(id,0,sizeof(spinspin)*loc_vol);
  //spinspin_put_to_id(id[0]);
  int nsources=1;
  
  start_loc_rnd_gen(100);
  memset(summ_corr,0,sizeof(corr16)*loc_vol);
  for(int isource=0;isource<nsources;isource++)
    {  
      generate_stochastic_source_and_tlSym_gluon_propagator(phi,eta,gl);
      generate_stochastic_A_B_dag_twisted_propagator(self_prop,prop,qu,phi,eta,gl);
      
      compute_all_2pts_qdagq_correlations(temp_corr,prop,self_prop);
      
      double d2=0,t2=0;
      nissa_loc_vol_loop(ivol)
        {
	  for(int ig=0;ig<16;ig++)
	    {
	      complex_summassign(summ_corr[ivol][ig],temp_corr[ivol][ig]);
	      complex_prod_double(temp_corr[ivol][ig],summ_corr[ivol][ig],1.0/(isource+1));
	    }
	  double codi=temp_corr[ivol][5][0]-corr[ivol][5][0];
	  double cosu=corr[ivol][5][0];
	  d2+=codi*codi;
	  t2+=cosu*cosu;
	}
      d2=glb_reduce_double(d2);
      t2=glb_reduce_double(t2);

      master_printf("time %d, isource=%d, diff=(%lg/%lg)=%lg\n",time(0),isource+1,d2,t2,d2/t2);
    }
  
  //////////////////////////////// compute correlation and write them on disk ////////////////////////

  write_corr16("self_energy_stoch_corr",temp_corr,64);
  
  close_calc();
  
  return 0;
}
