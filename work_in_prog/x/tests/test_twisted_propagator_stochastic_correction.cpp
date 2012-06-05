#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/propagators/Wilson_gluon_propagator.h"
#include "../src/propagators/twisted_propagator_g2_corr.h"
#include "../src/stochastic/stochastic_twisted_propagator.h"
#include "../src/types/types_routines.h"
#include "../src/routines/shift.h"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.h"

spinspin *d2_stoch_corr,*d2_stoch_corr_ave;
spinspin *d2_stoch_corr_new,*d2_stoch_corr_new_ave;
spinspin *d2_corr,*temp_corr;
spin1field *phi,*eta,*phi_sh,*eta_sh;
spinspin *q_prop,*q_prop_sh,*id;
spin1prop *g_prop,*g_prop_sh;

int map_mu[4]={4,1,2,3};

//initialize the program
void init_test()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(4,4);
  
  //allocate propagators
  id=nissa_malloc("1",loc_vol,spinspin);
  q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  g_prop_sh=nissa_malloc("g_prop_sh",loc_vol,spin1prop);
  q_prop_sh=nissa_malloc("q_prop_sh",loc_vol,spinspin);
  g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  d2_corr=nissa_malloc("d2_corr",loc_vol,spinspin);
  temp_corr=nissa_malloc("temp_corr",loc_vol,spinspin);
  d2_stoch_corr=nissa_malloc("d2_stoch_corr",loc_vol,spinspin);
  d2_stoch_corr_new=nissa_malloc("d2_stoch_corr_new",loc_vol,spinspin);
  d2_stoch_corr_new_ave=nissa_malloc("d2_stoch_corr_new_ave",loc_vol,spinspin);
  d2_stoch_corr_ave=nissa_malloc("d2_stoch_corr_ave",loc_vol,spinspin);
  //allocate stoch field propagator
  phi=nissa_malloc("phi",loc_vol,spin1field);
  eta=nissa_malloc("eta",loc_vol,spin1field);
  phi_sh=nissa_malloc("phi_sh",loc_vol,spin1field);
  eta_sh=nissa_malloc("eta_sh",loc_vol,spin1field);
}

//close the program
void close_test()
{
  nissa_free(id);
  nissa_free(q_prop);
  nissa_free(q_prop_sh);
  nissa_free(g_prop);
  nissa_free(g_prop_sh);
  nissa_free(d2_corr);
  nissa_free(temp_corr);
  nissa_free(d2_stoch_corr);
  nissa_free(d2_stoch_corr_new);
  nissa_free(d2_stoch_corr_new_ave);
  nissa_free(d2_stoch_corr_ave);
  nissa_free(eta);
  nissa_free(phi);
  nissa_free(eta_sh);
  nissa_free(phi_sh);
  
  close_nissa();
}

void summ_the_stoch_contribution(spinspin *osi,spinspin *q,spin1field *phi,spin1field *eta,int nu,int mu,spinspin *oso,double sign)
{
  nissa_loc_vol_loop(ivol)
  {
    complex g;
    unsafe_complex_conj2_prod(g,phi[ivol][nu],eta[0][mu]);
    
    spinspin t;
    unsafe_spinspin_spinspin_prod(t,q[ivol],oso[mu]);
    safe_spinspin_spinspin_prod(t,osi[nu],t);
    spinspin_prod_double(t,t,-sign/4);
    spinspin_summ_the_complex_prod(d2_stoch_corr[ivol],t,g);
  }
}

void summ_the_contribution(spinspin *osi,spinspin *q,spin1prop *g,int nu,int mu,spinspin *oso,double sign)
{
  nissa_loc_vol_loop(ivol)
  {
    spinspin t;
    unsafe_spinspin_spinspin_prod(t,q[ivol],oso[mu]);
    safe_spinspin_spinspin_prod(t,osi[nu],t);
    spinspin_prod_double(t,t,-sign/4);
    spinspin_summ_the_complex_prod(d2_corr[ivol],t,g[ivol][nu][mu]);
  }
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
  
  ///////////////////////////////// compute (1+-g_\mu) /////////////////////////////////
  
  spinspin opg[4],omg[4];
  
  for(int mu=0;mu<4;mu++)
    {
      spinspin_put_to_id(opg[mu]);
      spinspin_put_to_id(omg[mu]);
      
      spinspin_dirac_summ_the_prod_double(opg[mu],&(base_gamma[map_mu[mu]]),+1);
      spinspin_dirac_summ_the_prod_double(omg[mu],&(base_gamma[map_mu[mu]]),-1);
    }
  
  coords ix={0,0,0,1};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  /////////////////////////////// S0 propagator=id //////////////////////////
  
  memset(id,0,sizeof(spinspin)*loc_vol);
  spinspin_put_to_id(id[0]);
  
  memset(d2_stoch_corr,0,sizeof(spinspin)*loc_vol);
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_x_space_Wilson_gluon_propagator_by_fft(g_prop,gl);
  
  /////////////////////////////// D2 correction //////////////////////////
  
  //reset the corr
  memset(d2_corr,0,sizeof(spinspin)*loc_vol);
  for(int nu=0;nu<4;nu++)
    for(int mu=0;mu<4;mu++)
      {
        //  term 1: +(1-gnu) G(x-nu,0)   S(x-nu,mu) (1-gmu)
        shift_spinspin_sink_dw(g_prop_sh,g_prop,gl.bc,nu);
        //memcpy(g_prop_sh,g_prop_sh,sizeof(spin1prop)*loc_vol);
        shift_spinspin_sink_dw(q_prop_sh,q_prop,qu.bc,nu);
        shift_spinspin_source_up(q_prop_sh,q_prop_sh,qu.bc,mu);
        summ_the_contribution(omg,q_prop_sh,g_prop_sh,nu,mu,omg,+1.0);
        
        //  term 2: -(1+gnu) G(x,0)      S(x+nu,mu) (1-gmu)
        //memcpy(g_prop_sh,g_prop,sizeof(spin1prop)*loc_vol);
        //memcpy(g_prop_sh,g_prop_sh,sizeof(spin1prop)*loc_vol);
        shift_spinspin_sink_up(q_prop_sh,q_prop,qu.bc,nu);
        shift_spinspin_source_up(q_prop_sh,q_prop_sh,qu.bc,mu);
        summ_the_contribution(opg,q_prop_sh,g_prop,nu,mu,omg,-1.0);
        
        //  term 3: -(1-gnu) G(x-nu,-mu) S(x-nu,-mu) (1+gmu)
        shift_spinspin_sink_dw(g_prop_sh,g_prop,gl.bc,nu);
        shift_spinspin_source_dw(g_prop_sh,g_prop_sh,gl.bc,mu);
        shift_spinspin_sink_dw(q_prop_sh,q_prop,qu.bc,nu);
        shift_spinspin_source_dw(q_prop_sh,q_prop_sh,qu.bc,mu);
        summ_the_contribution(omg,q_prop_sh,g_prop_sh,nu,mu,opg,-1.0);
        
        //  term 4: +(1+gnu) G(x,-mu) S(x+nu,-mu) (1+gmu)
        //memcpy(g_prop_sh,g_prop,sizeof(spin1prop)*loc_vol);
        shift_spinspin_source_dw(g_prop_sh,g_prop,gl.bc,mu);
        shift_spinspin_sink_up(q_prop_sh,q_prop,qu.bc,nu);
        shift_spinspin_source_dw(q_prop_sh,q_prop_sh,qu.bc,mu);
        summ_the_contribution(opg,q_prop_sh,g_prop_sh,nu,mu,opg,+1.0);
      }
  
  pass_spinspin_from_x_to_mom_space(d2_corr,d2_corr,quark_theta);
  
  print_spinspin(d2_corr[lx]);
  master_printf("\n");
  
  /////////////////////////////// D2 stoch correction //////////////////////////
  
  int log2=0;
  for(int isource=0;isource<nsource;isource++)
    {
      generate_stochastic_source_and_tlSym_gluon_propagator(phi,eta,gl);
      
      generate_stochastic_A_B_dag_twisted_propagator_source(temp_corr,id,qu,phi,eta,gl);
      nissa_loc_vol_loop(ivol)
	spinspin_summ(d2_stoch_corr_new[ivol],d2_stoch_corr_new[ivol],temp_corr[ivol]);
            
      for(int nu=0;nu<4;nu++)
	for(int mu=0;mu<4;mu++)
	  {
	    //  term 1: +(1-gnu) G(x-nu,0)   S(x-nu,mu) (1-gmu)
	    shift_spin_up(phi_sh,phi,gl.bc,nu);
	    //memcpy(eta_sh,eta,sizeof(spin)*loc_vol);
	    shift_spinspin_sink_dw(q_prop_sh,q_prop,qu.bc,nu);
	    shift_spinspin_source_up(q_prop_sh,q_prop_sh,qu.bc,mu);
	    summ_the_stoch_contribution(omg,q_prop_sh,phi_sh,eta,nu,mu,omg,+1.0);
	    
	    //  term 2: -(1+gnu) G(x,0)      S(x+nu,mu) (1-gmu)
	    //memcpy(phi_sh,phi,sizeof(spin)*loc_vol);
	    //memcpy(eta_sh,eta,sizeof(spin)*loc_vol);
	    shift_spinspin_sink_up(q_prop_sh,q_prop,qu.bc,nu);
	    shift_spinspin_source_up(q_prop_sh,q_prop_sh,qu.bc,mu);
	    summ_the_stoch_contribution(opg,q_prop_sh,phi,eta,nu,mu,omg,-1.0);
	    
	    //  term 3: -(1-gnu) G(x-nu,-mu) S(x-nu,-mu) (1+gmu)
	    shift_spin_up(phi_sh,phi,gl.bc,nu);
	    shift_spin_up(eta_sh,eta,gl.bc,mu);
	    shift_spinspin_sink_dw(q_prop_sh,q_prop,qu.bc,nu);
	    shift_spinspin_source_dw(q_prop_sh,q_prop_sh,qu.bc,mu);
	    summ_the_stoch_contribution(omg,q_prop_sh,phi_sh,eta_sh,nu,mu,opg,-1.0);
	    
	    //  term 4: +(1+gnu) G(x,-mu) S(x+nu,-mu) (1+gmu)
	    //memcpy(phi_sh,phi,sizeof(spin)*loc_vol);
	    shift_spin_up(eta_sh,eta,gl.bc,mu);
	    shift_spinspin_sink_up(q_prop_sh,q_prop,qu.bc,nu);
	    shift_spinspin_source_dw(q_prop_sh,q_prop_sh,qu.bc,mu);
	    summ_the_stoch_contribution(opg,q_prop_sh,phi,eta_sh,nu,mu,opg,+1.0);
	  }
      
      //////////////////////////////////// output ////////////////////////////
      
      if((isource+1)==1<<log2)
	{
	  log2++;
	  nissa_loc_vol_loop(ivol)
	    {
	      spinspin_prod_double(d2_stoch_corr_ave[ivol],d2_stoch_corr[ivol],1.0/(isource+1));
	      spinspin_prod_double(d2_stoch_corr_new_ave[ivol],d2_stoch_corr_new[ivol],1.0/(isource+1));
	    }
	  
	  pass_spinspin_from_x_to_mom_space(d2_stoch_corr_ave,d2_stoch_corr_ave,quark_theta);
	  pass_spinspin_from_x_to_mom_space(d2_stoch_corr_new_ave,d2_stoch_corr_new_ave,quark_theta);
	 
	  master_printf("%d\n",isource+1);
	  master_printf("exa\n");
	  print_spinspin(d2_corr[lx]);
	  master_printf("new\n");
	  print_spinspin(d2_stoch_corr_new_ave[lx]);
	  master_printf("old\n");
	  print_spinspin(d2_stoch_corr_ave[lx]);
	  master_printf("\n");
	  
	  double rold=0,rnew=0;
	  nissa_loc_vol_loop(ivol)
	  {
	    spinspin t;
	    spinspin_subt(t,d2_stoch_corr_ave[ivol],d2_corr[ivol]);
	    rold+=real_part_of_trace_spinspin_prod_spinspin_dag(t,t);

	    spinspin_subt(t,d2_stoch_corr_new_ave[ivol],d2_corr[ivol]);
	    rnew+=real_part_of_trace_spinspin_prod_spinspin_dag(t,t);
	  }
	  master_printf("after %d sources: old %lg  new %lg\n",isource+1,sqrt(rold/glb_vol),sqrt(rnew/glb_vol));
	  fflush(stdout);
	  
	}
    }
  
  close_test();
  
  return 0;
}
