#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/propagators/tlSym_gluon_propagator.h"
#include "../src/propagators/Wilson_gluon_propagator.h"
#include "../src/types/types_routines.h"
#include "../src/routines/fourier.h"
#include "../src/routines/shift.h"

spinspin  *q_prop,*q_prop_sh;
spin1prop *g_prop,*g_prop_sh;
spinspin *corr_x,*corr_p;

int comp=01;
int map_mu[4]={4,1,2,3};

//initialize the program
void init_test()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(4,4);
  
  //allocate propagators
  q_prop_sh=nissa_malloc("q_prop_sh",loc_vol,spinspin);
  g_prop_sh=nissa_malloc("g_prop_sh",loc_vol,spin1prop);
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

void summ_the_contribution(spinspin *osi,spinspin *q,spin1prop *g,int nu,int mu,spinspin *oso,double sign)
{
  nissa_loc_vol_loop(ivol)
    {
      spinspin t;
      unsafe_spinspin_spinspin_prod(t,q[ivol],oso[mu]);
      safe_spinspin_spinspin_prod(t,osi[nu],t);
      spinspin_prod_double(t,t,sign);
      spinspin_summ_the_complex_prod(corr_x[ivol],t,g[ivol][nu][mu]);
    }
}

int main(int narg,char **arg)
{
  init_test();
  
  double null_theta[4]={0,0,0,0};
  double small=1.e-6;
  double small_theta[4]={small,small,small,small};
  double rand_theta[4]={0.1,0.3,0.6,0.4};
  
  //quark
  double quark_theta[4];memcpy(quark_theta,null_theta,sizeof(double)*4);
  double kappa=1.0/8;
  double mass=0;
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double gluon_theta[4];memcpy(gluon_theta,null_theta,sizeof(double)*4);
  double alpha=0.3;
  gluon_info gl=create_Wilson_gluon_info(alpha,gluon_theta);
  
  /////////////////////////////// quark and gluon propagator //////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_x_space_Wilson_gluon_propagator_by_fft(g_prop,gl);
  
  ///////////////////////////////// compute (1+-g_\mu) /////////////////////////////////
  
  spinspin opg[4],omg[4];
  
  for(int mu=0;mu<4;mu++)
    {
      spinspin_put_to_id(opg[mu]);
      spinspin_put_to_id(omg[mu]);
      
      spinspin_dirac_summ_the_prod_double(opg[mu],&(base_gamma[map_mu[mu]]),+1);
      spinspin_dirac_summ_the_prod_double(omg[mu],&(base_gamma[map_mu[mu]]),-1);
    }
  
  /////////////////////////////////// correction D2 ///////////////////////////////////
  
  //reset the corr
  memset(corr_x,0,sizeof(spinspin)*loc_vol);
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
  
  nissa_loc_vol_loop(ivol)
    spinspin_prod_double(corr_x[ivol],corr_x[ivol],0.25);
  
  pass_spinspin_from_x_to_mom_space(corr_x,corr_x,qu.bc);
  
  ///////////////////////////////// correction in P ////////////////////////////  
  
  if(rank_tot==1 && comp)
    {
      memset(corr_p,0,sizeof(spinspin)*loc_vol);
      
      compute_mom_space_twisted_propagator(q_prop,qu);
      compute_mom_space_Wilson_gluon_propagator(g_prop,gl);
      
      bool zmp=1;
      for(int mu=0;mu<4;mu++) zmp&=(gl.bc[mu]==0);
      
      if(zmp)
	nissa_loc_vol_loop(ip)
	  nissa_loc_vol_loop(iq)
	    {
	      spinspin w[4];
	      
	      coords ir_mu;
	      //compute the weights and find r
	      for(int mu=0;mu<4;mu++)
		{
		  //r=p-q
		  ir_mu[mu]=(glb_coord_of_loclx[ip][mu]-glb_coord_of_loclx[iq][mu]+glb_size[mu])%glb_size[mu];
		  
		  double p_mu=M_PI*(2*glb_coord_of_loclx[ip][mu]+qu.bc[mu])/glb_size[mu];
		  double q_mu=M_PI*(2*glb_coord_of_loclx[iq][mu]+gl.bc[mu])/glb_size[mu];
		  
		  double s_mu=(2*p_mu-q_mu)/2;
		  double co=cos(s_mu),si=sin(s_mu);
		  
		  spinspin_dirac_prod_double(w[mu],&(base_gamma[map_mu[mu]]),co);
		  spinspin_dirac_summ_the_prod_idouble(w[mu],&(base_gamma[0]),si);
		}
	      
	      int ir=glblx_of_coord(ir_mu);
	      for(int mu=0;mu<4;mu++)
		for(int nu=0;nu<4;nu++)
		  {
		    spinspin temp;
		    unsafe_spinspin_spinspin_prod(temp,w[mu],q_prop[ir]);
		    safe_spinspin_spinspin_prod(temp,temp,w[nu]);
		    spinspin_summ_the_complex_prod(corr_p[ip],temp,g_prop[iq][mu][nu]);
		  }
	    }
      else
	{
	  crash("Non periodic boundary conditions!");
	  /*
	  nissa_loc_vol_loop(ip)
	    nissa_loc_vol_loop(iq)
	      nissa_loc_vol_loop(ir)
	        {
		  spinspin w[4];
		  
		  //compute the weights and find pp
		  for(int mu=0;mu<4;mu++)
		    {
		  
		  double p_mu=M_PI*(2*glb_coord_of_loclx[ip][mu]+qu.bc[mu])/glb_size[mu];
		  double r_mu=M_PI*(2*glb_coord_of_loclx[ir][mu]+qu.bc[mu])/glb_size[mu];
		  
		  double s_mu=(p_mu+r_mu)/2;
		  double co=cos(s_mu),si=sin(s_mu);
		  
		  spinspin_dirac_prod_double(w[mu],&(base_gamma[map_mu[mu]]),co);
		  spinspin_dirac_summ_the_prod_idouble(w[mu],&(base_gamma[0]),si);
		}
	      
	      int ipp=glblx_of_coord(ipp_mu);
	      for(int mu=0;mu<4;mu++)
		for(int nu=0;nu<4;nu++)
		  {
		    spinspin temp;
		    unsafe_spinspin_spinspin_prod(temp,w[mu],q_prop[ipp]);
		    safe_spinspin_spinspin_prod(temp,temp,w[nu]);
		    spinspin_prod_double(temp,temp,4.0/glb_vol);
		    spinspin_summ_the_complex_prod(corr_p[ip],temp,g_prop[iq][mu][nu]);
		  }
	  */
	}
      //pass_spinspin_from_mom_to_x_space(corr_p,corr_p,qu.bc);
    }
  
  if(comp)
  nissa_loc_vol_loop(ivol)
    {
      spinspin temp;
      coords opp;
      for(int mu=0;mu<4;mu++) opp[mu]=(glb_size[mu]-glb_coord_of_loclx[ivol][mu])%glb_size[mu];
      int iopp=glblx_of_coord(opp);
      spinspin_subt(temp,corr_x[ivol],corr_p[ivol]);
      double t2=real_part_of_trace_spinspin_prod_spinspin_dag(temp,temp);
      if(fabs(t2)>1.e-10) 
	{
	  master_printf("%d %lg\n",ivol,t2);
	  print_spinspin(temp);
	  master_printf("x\n");
	  print_spinspin(corr_x[ivol]);
	  master_printf("p\n");
	  print_spinspin(corr_p[ivol]);
	}
    }

  //////////////////////////////////// output ////////////////////////////
  
  //compute the point to be printed
  
  coords ix={0,0,0,1};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  double a2p2=0,eps=8.239;
  for(int mu=0;mu<4;mu++)
    {
      double p=M_PI*(2*ix[mu]+qu.bc[mu])/glb_size[mu];
      a2p2+=p*p;
    }
  printf("a2p2: %lg\n",a2p2);
  printf("att: %lg\n",-4*a2p2*(eps-5.39*alpha-0.5*(3-2*alpha)*log(a2p2))/(16*M_PI*M_PI));
  
  if(rank_tot==1 && comp)
    {
      print_spinspin(corr_p[lx]);
      printf("\n");
    }
  
  if(rank==rx) print_spinspin(corr_x[lx]);
  
  close_test();
  
  return 0;
}
