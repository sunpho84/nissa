#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/diagrams/meson_exchange.h"
#include "../src/types/types_routines.h"
#include "../src/routines/read_and_write.h"
#include "../src/routines/correlations.h"
#include "../src/stochastic/stochastic_twisted_propagator.h"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.h"

spin1field *phi,*eta;
spinspin *prop,*prop_phi,*prop_eta;
corr16 *corrA,*corrB,*corrC,*summ_corr,*temp_corr;

int L=2;

//initialize the program
void init_calc()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(L,L);
  
  //allocate propagators
  prop=nissa_malloc("prop",loc_vol,spinspin);
  prop_phi=nissa_malloc("prop_phi",loc_vol,spinspin);
  prop_eta=nissa_malloc("prop_eta",loc_vol,spinspin);
  summ_corr=nissa_malloc("summ_corr",loc_vol,corr16);
  temp_corr=nissa_malloc("temp_corr",loc_vol,corr16);
  corrA=nissa_malloc("corrA",loc_vol,corr16);
  corrB=nissa_malloc("corrB",loc_vol,corr16);
  corrC=nissa_malloc("corrC",loc_vol,corr16);
  phi=nissa_malloc("phi",loc_vol+bord_vol,spin1field);
  eta=nissa_malloc("eta",loc_vol+bord_vol,spin1field);
}

//close the program
void close_calc()
{
  nissa_free(prop);
  nissa_free(prop_phi);
  nissa_free(prop_eta);
  nissa_free(corrA);
  nissa_free(corrB);
  nissa_free(corrC);
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
  //double quark_theta[4]={-0.3,-0.3,-0.3,-0.3};
  //double quark_theta[4]={0.3,0.3,0.3,0.3};
  double quark_theta[4]={0,0,0,0};
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double alpha=0;
  double gluon_theta[4]={0,0,0,0};
  //double gluon_theta[4]={0.3,0.3,0.3,0.3};
  gluon_info gl=create_tlSym_gluon_info(alpha,gluon_theta);
  
  ////////////////////////////////////// compute correlation analytically ////////////////////////////
  
  compute_meson_exchange_correction_analyticallyA(corrA,qu,gl);
  write_corr16("exchange_corr",corrA,64);
  
  compute_meson_exchange_correction_analyticallyB(corrB,qu,gl);

  compute_meson_exchange_correction_analyticallyC(corrC,qu,gl);
  
  //print the various determinations of P5P5
  master_printf("%s\t%s\t%s\n","XA","XB","P");
  nissa_loc_vol_loop(ivol)
    {
      int ig=5;
      
      master_printf("%lg\t%lg\t%lg\n",corrA[ivol][ig][0],corrB[ivol][ig][0],corrC[ivol][ig][0]);
    }

  ////////////////////////////////////// compute correlation stochastically ////////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(prop,qu);
  
  int nsources=1;
  const int corr_id=0;
  start_loc_rnd_gen(100);
  memset(summ_corr,0,sizeof(corr16)*loc_vol);
  for(int isource=0;isource<nsources;isource++)
    {  
      generate_stochastic_source_and_tlSym_gluon_propagator(phi,eta,gl);
      generate_stochastic_A_twisted_propagator(prop_phi,prop,qu,phi,gl);
      generate_stochastic_A_twisted_propagator(prop_eta,prop,qu,eta,gl);
      
      compute_all_2pts_qdagq_correlations(temp_corr,prop_phi,prop_eta);
      
      double d2=0,t2=0;
      nissa_loc_vol_loop(ivol)
        {
	  for(int ig=0;ig<16;ig++)
	    {
	      complex_summassign(summ_corr[ivol][ig],temp_corr[ivol][ig]);
	      complex_prod_double(temp_corr[ivol][ig],summ_corr[ivol][ig],1.0/(isource+1));
	    }
	  double codi=temp_corr[ivol][corr_id][0]-corrA[ivol][corr_id][0];
	  if(isource==nsources-1) master_printf("%d %lg\t%lg\t%lg\t%lg\n",ivol,temp_corr[ivol][corr_id][0],corrA[ivol][corr_id][0],corrB[ivol][corr_id][0],corrC[ivol][corr_id][0]);
	  double cosu=corrA[ivol][corr_id][0];
	  d2+=codi*codi;
	  t2+=cosu*cosu;
	}
      
      d2=glb_reduce_double(d2);
      t2=glb_reduce_double(t2);
      
      if(isource==nsources-1)
	{
	  master_printf("isource=%d, diff=(%lg/%lg)=%lg\n",isource+1,d2,t2,d2/t2);
	  if(d2/t2<1.5e-5) master_printf("OK!\n");
	  else master_printf("ERROR!!!!\n");
	}
    }
  
  //////////////////////////////// compute correlation and write them on disk ////////////////////////
  
  write_corr16("exchange_stoch_corr",temp_corr,64);
  
  close_calc();
  
  return 0;
}
