#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/propagators/tlSym_gluon_propagator.h"
#include "../src/types/types_routines.h"
#include "../src/routines/fourier.h"
#include "../src/routines/shift.h"

spinspin  *q_prop,*q_prop_sh;
spin1prop *g_prop,*g_prop_sh;
spinspin *corr_x,*corr_p;

//initialize the program
void init_test()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(8,4);
  
  //allocate propagators
  q_prop_sh=nissa_malloc("q_prop_sh",loc_vol,spinspin);
  g_prop_sh=nissa_malloc("g_prop_sh",loc_vol,spinspin);
  q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  corr_x=nissa_malloc("corr_x",loc_vol,spinspin);
  corr_p=nissa_malloc("corr_p",loc_vol,spinspin);
}

//close the program
void close_test()
{
  nissa_free(q_prop_sh);
  nissa_free(g_prop_sh);
  nissa_free(q_prop);
  nissa_free(g_prop);
  nissa_free(corr_x);
  nissa_free(corr_p);
  
  close_nissa();
}

void summ_the_contribution(spinspin a,spinspin *q,spin1prop *g,int nu,int mu,spinspin b,double sign)
{
  nissa_loc_vol_loop(ivol)
    {
      spinspin t;
      unsafe_spinspin_spinspin_prod(t,a,q[ivol]);
      safe_spinspin_spinspin_prod(t,t,b);
      spinspin_prod_double(t,t,sign);
      spinspin_summ_the_complex_prod(corr_x[ivol],t,g[ivol][nu][mu]);
    }
}

int main(int narg,char **arg)
{
  init_test();
  
  //quark
  double quark_theta[4]={1,1,1,1};
  double kappa=1.0/8;
  double mass=0;
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double gluon_theta[4]={1,1,1,1};
  double alpha=0.3;
  gluon_info gl=create_Wilson_gluon_info(alpha,gluon_theta);
  
  /////////////////////////////// quark and gluon propagator //////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);
  
  ///////////////////////////////// compute (1+-g_\mu) /////////////////////////////////
  
  spinspin opg[4],omg[4];
  
  for(int mu=0;mu<4;mu++)
    {
      int nu=(mu==0)?4:mu;
      
      spinspin_put_to_id(opg[mu]);
      spinspin_put_to_id(omg[mu]);
      
      spinspin_dirac_summ_the_prod_double(opg[mu],&(base_gamma[nu]),+1);
      spinspin_dirac_summ_the_prod_double(omg[mu],&(base_gamma[nu]),-1);
    }
  
  /////////////////////////////////// correction D2 ///////////////////////////////////
  
  //reset the corr
  memset(corr_x,0,sizeof(spinspin)*loc_vol);
  
  for(int nu=0;nu<4;nu++)
    for(int mu=0;mu<4;mu++)
      {
	//  term 1: +(1-gnu) G(x,mu)    S(x-nu,mu) (1-gmu)
	//memcpy(g_prop_sh,g_prop,sizeof(spin1prop)*loc_vol);
	shift_spinspin_source_up(g_prop_sh,g_prop,gluon_theta,mu);
	shift_spinspin_sink_dw(q_prop_sh,q_prop,quark_theta,nu);
	shift_spinspin_source_up(q_prop_sh,q_prop_sh,quark_theta,mu);
	summ_the_contribution(omg[nu],q_prop_sh,g_prop_sh,nu,mu,omg[mu],+1.0);
	
	//  term 2: -(1+gnu) G(x+nu,mu) S(x+nu,mu) (1-gmu)
	shift_spinspin_sink_up(g_prop_sh,g_prop,gluon_theta,nu);
	shift_spinspin_source_up(g_prop_sh,g_prop_sh,gluon_theta,mu);
	shift_spinspin_sink_up(q_prop_sh,q_prop,quark_theta,nu);
	shift_spinspin_source_up(q_prop_sh,q_prop_sh,quark_theta,mu);
	summ_the_contribution(opg[nu],q_prop_sh,g_prop_sh,nu,mu,omg[mu],-1.0);
	
	//  term 3: -(1-gnu) G(x,0)     S(x-nu,-mu) (1+gmu)
	//memcpy(g_prop_sh,g_prop,sizeof(spin1prop)*loc_vol);
	//memcpy(g_prop_sh,g_prop_sh,sizeof(spin1prop)*loc_vol);
	shift_spinspin_sink_dw(q_prop_sh,q_prop,quark_theta,nu);
	shift_spinspin_source_dw(q_prop_sh,q_prop_sh,quark_theta,mu);
	summ_the_contribution(omg[nu],q_prop_sh,g_prop,nu,mu,opg[mu],-1.0);
	
	//  term 4: +(1+gnu) G(x+nu,0)  S(x+nu,-mu) (1+gmu)
	shift_spinspin_sink_up(g_prop_sh,g_prop,gluon_theta,nu);
	//memcpy(g_prop_sh,g_prop_sh,sizeof(spin1prop)*loc_vol);
	shift_spinspin_sink_up(q_prop_sh,q_prop,quark_theta,nu);
	shift_spinspin_source_dw(q_prop_sh,q_prop_sh,quark_theta,mu);
	summ_the_contribution(opg[nu],q_prop_sh,g_prop_sh,nu,mu,opg[mu],+1.0);
      }

  ///////////////////////////////// correction in P ////////////////////////////  
  
  if(rank_tot==1)
    {
      memset(corr_p,0,sizeof(spinspin)*loc_vol);
      
      compute_mom_space_twisted_propagator(q_prop,qu);
      compute_mom_space_tlSym_gluon_propagator(g_prop,gl);
      
      nissa_loc_vol_loop(ip)
	nissa_loc_vol_loop(iq)
	{
	  spinspin w[4];
	  
	  coords ipp_mu;
	  //compute the weights and find pp
	  for(int mu=0;mu<4;mu++)
	    {
	      //pp=p-q
	      ipp_mu[mu]=(glb_coord_of_loclx[ip][mu]-glb_coord_of_loclx[iq][mu]+glb_size[mu])%glb_size[mu];

	      double p_mu=M_PI*(2*glb_coord_of_loclx[ip][mu]+qu.bc[mu])/glb_size[mu];
	      double q_mu=M_PI*(2*glb_coord_of_loclx[iq][mu]+gl.bc[mu])/glb_size[mu];
	      
	      double c=2*p_mu-q_mu;
	      double a=4*cos(c/2.0)/glb_vol,b=4*sin(c/2.0)/glb_vol;

	      int nu=(mu==0)?4:mu;
	      spinspin_put_to_zero(w[mu]);
	      for(int ig=0;ig<4;ig++)
		{
		  complex_summ_the_prod_double( w[mu][ig][base_gamma[nu].pos[ig]],base_gamma[nu].entr[ig],a);
		  complex_subt_the_prod_idouble(w[mu][ig][base_gamma[0].pos[ig]],base_gamma[0].entr[ig],b);
		}
	    }
	  
	  int ipp=glblx_of_coord(ipp_mu);
	  for(int mu=0;mu<4;mu++)
	    for(int nu=0;nu<4;nu++)
	      {
		spinspin temp;
		unsafe_spinspin_spinspin_prod(temp,w[mu],q_prop[ipp]);
		safe_spinspin_spinspin_prod(temp,temp,w[nu]);
		spinspin_summ_the_complex_prod(corr_p[ip],temp,g_prop[iq][mu][nu]);
	      }
	}
      
      pass_spinspin_from_mom_to_x_space(corr_p,corr_p,quark_theta);
    }
      
  //////////////////////////////////// output ////////////////////////////
  
  //compute the point to be printed
  coords ix={5,2,1,2};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  lx=0;
  if(rank_tot==1)
    {
      print_spinspin(corr_p[lx]);
      printf("\n");
    }
  
  if(rank==rx) print_spinspin(corr_x[lx]);
  
  close_test();
  
  return 0;
}
