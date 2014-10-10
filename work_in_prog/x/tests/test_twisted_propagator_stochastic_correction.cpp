#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/propagators/Wilson_gluon_propagator.hpp"
#include "../src/propagators/twisted_propagator_g2_corr.hpp"
#include "../src/stochastic/stochastic_twisted_propagator.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/routines/shift.hpp"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.hpp"
#include "../src/diagrams/propagator_self_energy.hpp"

spinspin *d2_stoch_corr,*d2_stoch_corr_ave;
spinspin *d2_corr,*temp_corr;
spin1field *phi,*eta;
spinspin *q_prop,*id;
spin1prop *g_prop;

//initialize the program
void init_test(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //init the grid
  init_grid(12,12);
  
  //allocate propagators
  id=nissa_malloc("1",loc_vol+bord_vol,spinspin);
  q_prop=nissa_malloc("q_prop",loc_vol+bord_vol,spinspin);
  g_prop=nissa_malloc("g_prop",loc_vol+bord_vol,spin1prop);
  d2_corr=nissa_malloc("d2_corr",loc_vol,spinspin);
  temp_corr=nissa_malloc("temp_corr",loc_vol,spinspin);
  d2_stoch_corr=nissa_malloc("d2_stoch_corr",loc_vol,spinspin);
  d2_stoch_corr_ave=nissa_malloc("d2_stoch_corr_ave",loc_vol,spinspin);
  //allocate stoch field propagator
  phi=nissa_malloc("phi",loc_vol+bord_vol,spin1field);
  eta=nissa_malloc("eta",loc_vol+bord_vol,spin1field);
}

//close the program
void close_test()
{
  nissa_free(id);
  nissa_free(q_prop);
  nissa_free(g_prop);
  nissa_free(d2_corr);
  nissa_free(temp_corr);
  nissa_free(d2_stoch_corr);
  nissa_free(d2_stoch_corr_ave);
  nissa_free(eta);
  nissa_free(phi);
  
  close_nissa();
}

void summ_the_stoch_contribution(spinspin *osi,spinspin *q,spin1field *phi,spin1field *eta,int nu,int mu,spinspin *oso,double sign)
{
  NISSA_LOC_VOL_LOOP(ivol)
  {
    complex g;
    unsafe_complex_conj2_prod(g,phi[ivol][nu],eta[0][mu]);
    
    spinspin t;
    unsafe_spinspin_prod_spinspin(t,q[ivol],oso[mu]);
    safe_spinspin_prod_spinspin(t,osi[nu],t);
    spinspin_prod_double(t,t,-sign/4);
    spinspin_summ_the_complex_prod(d2_stoch_corr[ivol],t,g);
  }
}

void set_spinspin_to_delta(spinspin *in)
{
  memset(in,0,sizeof(spinspin)*loc_vol);
  if(rank==0) spinspin_put_to_id(in[0]);
  set_borders_invalid(in);
}

int main(int narg,char **arg)
{
  int nsource=64*64*64*64;
  
  init_test(narg,arg);
  
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
  
  coords ix={0,0,0,1};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  compute_self_energy_twisted_diagram_in_x_space(d2_corr,qu,gl);
  pass_spinspin_from_x_to_mom_space(d2_corr,d2_corr,quark_theta);
  
  /////////////////////////////// S0 propagator=id //////////////////////////
  
  set_spinspin_to_delta(id);
  
  memset(d2_stoch_corr,0,sizeof(spinspin)*loc_vol);
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_x_space_Wilson_gluon_propagator_by_fft(g_prop,gl);
  
  /////////////////////////////// D2 stoch correction //////////////////////////
  
  int log2=0;
  for(int isource=0;isource<nsource;isource++)
    {
      generate_stochastic_source_and_tlSym_gluon_propagator(phi,eta,gl);
      generate_stochastic_A_B_dag_twisted_propagator_source(temp_corr,id,qu,phi,eta,gl);

      NISSA_LOC_VOL_LOOP(ivol)
	spinspin_summassign(d2_stoch_corr[ivol],temp_corr[ivol]);
      
      //////////////////////////////////// output ////////////////////////////
      
      if((isource+1)==1<<log2)
	{
	  log2++;
	  NISSA_LOC_VOL_LOOP(ivol)
	    spinspin_prod_double(d2_stoch_corr_ave[ivol],d2_stoch_corr[ivol],1.0/(isource+1));
	  
	  pass_spinspin_from_x_to_mom_space(d2_stoch_corr_ave,d2_stoch_corr_ave,quark_theta);
	 
	  if(rank==rx)
	    {
	      master_printf("%d\n",isource+1);
	      master_printf("exa\n");
	      spinspin_print(d2_corr[lx]);
	      master_printf("sto\n");
	      spinspin_print(d2_stoch_corr_ave[lx]);
	      master_printf("\n");
	    }
	  
	  double rs=0,rd=0;
	  NISSA_LOC_VOL_LOOP(ivol)
	  {
	    spinspin td;
	    spinspin_subt(td,d2_stoch_corr_ave[ivol],d2_corr[ivol]);
	    rd+=real_part_of_trace_spinspin_prod_spinspin_dag(td,td);
	    rs+=real_part_of_trace_spinspin_prod_spinspin_dag(d2_corr[ivol],d2_corr[ivol]);
	  }
	  rs=glb_reduce_double(rs)/glb_vol;
	  rd=glb_reduce_double(rd)/glb_vol;
	  
	  if(rank==rx)
	    {
	      printf("Relative difference after %d sources: %lg\n",isource+1,sqrt(rd/rs));
	      fflush(stdout);
	    }
	}
    }
  
  close_test();
  
  return 0;
}
