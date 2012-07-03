#include <math.h>

#include "nissa.h"

#include "../propagators/twisted_propagator.h"
#include "../propagators/tlSym_gluon_propagator.h"
#include "../routines/fourier.h"
#include "../routines/shift.h"

#include "propagator_self_energy.h"

/*
     ___/____A____/___
    /   \   {_    \   \
   /         _}        \
  X         {_          O
   \          }        /
    \___\____{_____\__/
        /    B     /
 */

//#define AMBA

//#ifdef AMBA

void summ_the_exchange_contributionA(corr16 corr,spinspin OB,spinspin pB,spinspin BX,spinspin XA,spinspin pA,spinspin AO,complex AB,double w)
{
  spinspin pO,XO; //upper line
  unsafe_spinspin_spinspin_prod(pO, pA,AO);
  unsafe_spinspin_spinspin_prod(XO, XA,pO);
  
  spinspin pX,OX; //lower line
  unsafe_spinspin_spinspin_prod(pX, pB,BX);
  unsafe_spinspin_spinspin_prod(OX, OB,pX);
  
  //no debug: commented product by AB
  complex wAB;
  complex_prod_double(wAB,AB,w);
  //complex wAB={w,0};

  spinspin XO_AB; //upper line with gluon and weight
  unsafe_spinspin_complex_prod(XO_AB,XO,wAB);
  
  //trace
  for(int ig=0;ig<16;ig++)
    {
      spinspin GOX,GXO_AB;
      spinspin_dirac_spinspin_prod(GOX,base_gamma+ig,OX);
      spinspin_dirac_spinspin_prod(GXO_AB,base_gamma+ig,XO_AB);
      summ_the_trace_prod_spinspins(corr[ig],GOX,GXO_AB);
    }
}

void compute_meson_exchange_correction_analyticallyA(corr16 *corr,quark_info qu,gluon_info gl)
{
  if(rank_tot>1) crash("works only on scalar");
  
  //freeing memory
  memset(corr,0,sizeof(corr16)*loc_vol);
  
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spinspin *temp=nissa_malloc("temp",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);
  
  int O=0;
  nissa_loc_vol_loop(A)
    for(int mu_A=0;mu_A<4;mu_A++)
      nissa_loc_vol_loop(X)
	nissa_loc_vol_loop(B)
	  for(int mu_B=0;mu_B<4;mu_B++)
	    {
	      int Aup=loclx_neighup[A][mu_A];
	      int Adw=loclx_neighdw[A][mu_A];
	      int Bup=loclx_neighup[B][mu_B];
	      int Bdw=loclx_neighdw[B][mu_B];
	      
	      spinspin OB,Bup_X,Bdw_X,XA,Aup_O,Adw_O;
	      spin1prop AB; //this is different for each contrib
	      
	      compute_x_space_propagator_to_sink_from_source(Adw_O,q_prop,qu.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[O]);
	      compute_x_space_propagator_to_sink_from_source(Aup_O,q_prop,qu.bc,glb_coord_of_loclx[Aup],glb_coord_of_loclx[O]);
	      compute_x_space_propagator_to_sink_from_source(XA,   q_prop,qu.bc,glb_coord_of_loclx[X],  glb_coord_of_loclx[A]);
	      compute_x_space_propagator_to_sink_from_source(Bdw_X,q_prop,qu.bc,glb_coord_of_loclx[Bdw],glb_coord_of_loclx[X]);
	      compute_x_space_propagator_to_sink_from_source(Bup_X,q_prop,qu.bc,glb_coord_of_loclx[Bup],glb_coord_of_loclx[X]);
	      compute_x_space_propagator_to_sink_from_source(OB,   q_prop,qu.bc,glb_coord_of_loclx[O],  glb_coord_of_loclx[B]);
	      
	      // - G S(0,B) (1-g_mu_B) S(Bup,X) G S(X,A) (1-g_mu_A) S(Aup,0) G(A,B)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[A],glb_coord_of_loclx[B]);
	      summ_the_exchange_contributionA(corr[X],OB,nissa_omg[mu_B],Bup_X,XA,nissa_omg[mu_A],Aup_O,AB[mu_A][mu_B],-0.25);
	      // + G S(0,B) (1-g_mu_B) S(Bup,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,B)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[B]);
	      summ_the_exchange_contributionA(corr[X],OB,nissa_omg[mu_B],Bup_X,XA,nissa_opg[mu_A],Adw_O,AB[mu_A][mu_B],+0.25);
	      // + G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1-g_mu_A) S(Aup,0) G(A,Bdw)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[A],glb_coord_of_loclx[Bdw]);
	      summ_the_exchange_contributionA(corr[X],OB,nissa_opg[mu_B],Bdw_X,XA,nissa_omg[mu_A],Aup_O,AB[mu_A][mu_B],+0.25);
	      // - G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,Bdw)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[Bdw]);
	      summ_the_exchange_contributionA(corr[X],OB,nissa_opg[mu_B],Bdw_X,XA,nissa_opg[mu_A],Adw_O,AB[mu_A][mu_B],-0.25);
	    }
  
  nissa_free(q_prop);
  nissa_free(temp);
  nissa_free(g_prop);
}

//#else

void summ_the_exchange_contributionB(corr16 corr,spinspin XB,spinspin pB,spinspin BO,spinspin XA,spinspin pA,spinspin AO,complex AB,double w)
{
  spinspin pO_up,XO_up; //upper line
  unsafe_spinspin_spinspin_prod(pO_up, pA,AO);
  unsafe_spinspin_spinspin_prod(XO_up, XA,pO_up);
  
  spinspin pO_dw,XO_dw; //lower line
  unsafe_spinspin_spinspin_prod(pO_dw, pB,BO);
  unsafe_spinspin_spinspin_prod(XO_dw, XB,pO_dw);
  
  spinspin OX_dw; //revert lower line
  unsafe_spinspin_hermitian(OX_dw,XO_dw);
  //safe_spinspin_prod_dirac(OX_dw,OX_dw,base_gamma+5);
  //safe_dirac_prod_spinspin(OX_dw,base_gamma+5,OX_dw);
  
  complex wAB;
  complex_prod_double(wAB,AB,w);
  
  //trace
  for(int ig=0;ig<16;ig++)
    {
      //no debug: comment prodcut by wAB
      spinspin GOX_dw;
      spinspin_dirac_spinspin_prod(GOX_dw,base_gamma+ig,OX_dw);
      
      spinspin GXO_up;
      spinspin_dirac_spinspin_prod(GXO_up,base_gamma+ig,XO_up);
      
      spinspin t;
      unsafe_spinspin_spinspin_prod(t,GXO_up,GOX_dw);
      complex c;
      trace_spinspin(c,t);
      complex_summ_the_prod(corr[ig],c,wAB);
      //complex_summ_the_prod_double(corr[ig],c,w);
    }
}

void compute_meson_exchange_correction_analyticallyB(corr16 *corr,quark_info qu,gluon_info gl)
{
  if(rank_tot>1) crash("works only on scalar");
  
  memset(corr,0,sizeof(corr16)*loc_vol);
  
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spinspin *temp=nissa_malloc("temp",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);

  int O=0;
  nissa_loc_vol_loop(A)
    for(int mu_A=0;mu_A<4;mu_A++)
      nissa_loc_vol_loop(X)
	nissa_loc_vol_loop(B)
	  for(int mu_B=0;mu_B<4;mu_B++)
	    {
	      int Aup=loclx_neighup[A][mu_A];
	      int Adw=loclx_neighdw[A][mu_A];
	      int Bup=loclx_neighup[B][mu_B];
	      int Bdw=loclx_neighdw[B][mu_B];
	      
	      spinspin XB,Bup_O,Bdw_O,XA,Aup_O,Adw_O;
	      spin1prop AB; //this is different for each contrib
	      
	      compute_x_space_propagator_to_sink_from_source(Adw_O,q_prop,qu.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[O]);
	      compute_x_space_propagator_to_sink_from_source(Aup_O,q_prop,qu.bc,glb_coord_of_loclx[Aup],glb_coord_of_loclx[O]);
	      compute_x_space_propagator_to_sink_from_source(XA,   q_prop,qu.bc,glb_coord_of_loclx[X],  glb_coord_of_loclx[A]);
	      compute_x_space_propagator_to_sink_from_source(Bdw_O,q_prop,qu.bc,glb_coord_of_loclx[Bdw],glb_coord_of_loclx[O]);
	      compute_x_space_propagator_to_sink_from_source(Bup_O,q_prop,qu.bc,glb_coord_of_loclx[Bup],glb_coord_of_loclx[O]);
	      compute_x_space_propagator_to_sink_from_source(XB,   q_prop,qu.bc,glb_coord_of_loclx[X],  glb_coord_of_loclx[B]);
	      
	      // - G S(0,B) (1-g_mu_B) S(Bup,X) G S(X,A) (1-g_mu_A) S(Aup,0) G(A,B)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[A],glb_coord_of_loclx[Bup]);
	      summ_the_exchange_contributionB(corr[X],XB,nissa_omg[mu_B],Bup_O,XA,nissa_omg[mu_A],Aup_O,AB[mu_A][mu_B],-0.25);
	      // + G S(0,B) (1-g_mu_B) S(Bup,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,B)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[Bup]);
	      summ_the_exchange_contributionB(corr[X],XB,nissa_omg[mu_B],Bup_O,XA,nissa_opg[mu_A],Adw_O,AB[mu_A][mu_B],+0.25);
	      // + G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1-g_mu_A) S(Aup,0) G(A,Bdw)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[A],glb_coord_of_loclx[B]);
	      summ_the_exchange_contributionB(corr[X],XB,nissa_opg[mu_B],Bdw_O,XA,nissa_omg[mu_A],Aup_O,AB[mu_A][mu_B],+0.25);
	      // - G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,Bdw)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[B]);
	      summ_the_exchange_contributionB(corr[X],XB,nissa_opg[mu_B],Bdw_O,XA,nissa_opg[mu_A],Adw_O,AB[mu_A][mu_B],-0.25);
	    }
  
  nissa_free(q_prop);
  nissa_free(temp);
  nissa_free(g_prop);
}

//#endif

int site_comb(int b,int wb,int c,int wc)
{
  coords co;
  for(int mu=0;mu<4;mu++)
    {
      co[mu]=glb_coord_of_loclx[b][mu]*wb+glb_coord_of_loclx[c][mu]*wc;
      while(co[mu]<0) co[mu]+=glb_size[mu];
      co[mu]%=glb_size[mu];
    }
  
  return loclx_of_coord(co);
}

void mom_space_qq_vertex_function(spinspin v,int imom_sum,quark_info qu,int mu)
{
  if(mu<0||mu>3) crash("mu=%d",mu);

  double p=M_PI*(2*glb_coord_of_loclx[imom_sum][mu]+qu.bc[mu])/glb_size[mu];
  double ph=p/2;

  spinspin_put_to_id(v);
  spinspin_prodassign_double(v,-sin(ph));
  spinspin_dirac_summ_the_prod_idouble(v,base_gamma+nissa_map_mu[mu],-cos(ph));
  spinspin_prodassign_double(v,glb_vol);
}

void compute_meson_exchange_correction_analyticallyC(corr16 *corr,quark_info qu,gluon_info gl)
{
  if(rank_tot>1) crash("works only on scalar");
  
  //resetting memory
  memset(corr,0,sizeof(corr16)*loc_vol);
  
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  
  compute_mom_space_twisted_propagator(q_prop,qu);
  compute_mom_space_tlSym_gluon_propagator(g_prop,gl);
  
  nissa_loc_vol_loop(p)
    nissa_loc_vol_loop(q)
      nissa_loc_vol_loop(r)
        for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++)
	    {
	      int ppq=site_comb(p,+1,q,+1);
	      int ppr=site_comb(p,+1,r,+1);
	      int rmq=site_comb(r,+1,q,-1);
	      
	      spinspin vmu,vnu;
	      mom_space_qq_vertex_function(vmu,site_comb(ppq,+1,ppr,+1),qu,mu);
	      mom_space_qq_vertex_function(vnu,site_comb(q,+1,r,+1),qu,nu);
	      
	      spinspin up;
	      unsafe_spinspin_spinspin_prod(up,vmu,q_prop[ppq]);
	      safe_spinspin_spinspin_prod(up,q_prop[ppr],up);
	      
	      spinspin dw;
	      unsafe_spinspin_spinspin_prod(dw,vnu,q_prop[r]);
	      safe_spinspin_spinspin_prod(dw,q_prop[q],dw);
	      
	      spinspin up_w;
	      unsafe_spinspin_complex_prod(up_w,up,g_prop[rmq][mu][nu]);
	      
	      for(int ig=0;ig<16;ig++)
		{
		  spinspin gup_w,gupg_w;
		  spinspin_dirac_spinspin_prod(gup_w,base_gamma+ig,up_w);
		  safe_spinspin_prod_dirac(gupg_w,gup_w,base_gamma+ig);
		  summ_the_trace_prod_spinspins(corr[p][ig],gupg_w,dw);
		}
	    }
  
  pass_spinspin_from_mom_to_x_space((spinspin*)corr,(spinspin*)corr,qu.bc);

  nissa_free(q_prop);
  nissa_free(g_prop);
}
