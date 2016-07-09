#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../propagators/twisted_propagator.hpp"
#include "../propagators/tlSym_gluon_propagator.hpp"
#include "../routines/fourier.hpp"
#include "../routines/correlations.hpp"
#include "../routines/shift.hpp"
#include "../vertex/vertex.hpp"
#include "../stochastic/stochastic_twisted_propagator.hpp"
#include "../stochastic/stochastic_tlSym_gluon_propagator.hpp"

#include "propagator_self_energy.hpp"

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
  unsafe_spinspin_prod_spinspin(pO, pA,AO);
  unsafe_spinspin_prod_spinspin(XO, XA,pO);
  
  spinspin pX,OX; //lower line
  unsafe_spinspin_prod_spinspin(pX, pB,BX);
  unsafe_spinspin_prod_spinspin(OX, OB,pX);
  
  complex wAB;
  complex_prod_double(wAB,AB,w);

  spinspin XO_AB; //upper line with gluon and weight
  unsafe_spinspin_complex_prod(XO_AB,XO,wAB);
  
  //trace
  for(int ig=0;ig<16;ig++)
    {
      spinspin GOX,GXO_AB;
      unsafe_dirac_prod_spinspin(GOX,base_gamma+ig,OX);
      unsafe_dirac_prod_spinspin(GXO_AB,base_gamma+ig,XO_AB);
      summ_the_trace_prod_spinspins(corr[ig],GOX,GXO_AB);
    }
}

void compute_meson_exchange_correction_analyticallyA(corr16 *corr,quark_info qu,gluon_info gl)
{
  if(nranks>1) CRASH("works only on scalar");
  
  //freeing memory
  memset(corr,0,sizeof(corr16)*loc_vol);
  
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spinspin *temp=nissa_malloc("temp",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);
  
  int O=0;
  NISSA_LOC_VOL_LOOP(A)
    for(int mu_A=0;mu_A<4;mu_A++)
      NISSA_LOC_VOL_LOOP(X)
	NISSA_LOC_VOL_LOOP(B)
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
	      summ_the_exchange_contributionA(corr[X],OB,omg[mu_B],Bup_X,XA,omg[mu_A],Aup_O,AB[mu_A][mu_B],-0.25);
	      // + G S(0,B) (1-g_mu_B) S(Bup,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,B)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[B]);
	      summ_the_exchange_contributionA(corr[X],OB,omg[mu_B],Bup_X,XA,opg[mu_A],Adw_O,AB[mu_A][mu_B],+0.25);
	      // + G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1-g_mu_A) S(Aup,0) G(A,Bdw)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[A],glb_coord_of_loclx[Bdw]);
	      summ_the_exchange_contributionA(corr[X],OB,opg[mu_B],Bdw_X,XA,omg[mu_A],Aup_O,AB[mu_A][mu_B],+0.25);
	      // - G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,Bdw)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[Bdw]);
	      summ_the_exchange_contributionA(corr[X],OB,opg[mu_B],Bdw_X,XA,opg[mu_A],Adw_O,AB[mu_A][mu_B],-0.25);
	    }
  
  nissa_free(q_prop);
  nissa_free(temp);
  nissa_free(g_prop);
}

//#else

void summ_the_exchange_contributionB(corr16 corr,spinspin XB,spinspin pB,spinspin BO,spinspin XA,spinspin pA,spinspin AO,complex AB,double w)
{
  spinspin pO_up,XO_up; //upper line
  unsafe_spinspin_prod_spinspin(pO_up, pA,AO);
  unsafe_spinspin_prod_spinspin(XO_up, XA,pO_up);
  
  spinspin pO_dw,XO_dw; //lower line
  unsafe_spinspin_prod_spinspin(pO_dw, pB,BO);
  unsafe_spinspin_prod_spinspin(XO_dw, XB,pO_dw);
  
  spinspin OX_dw; //revert lower line
  unsafe_spinspin_hermitian(OX_dw,XO_dw);
  safe_spinspin_prod_dirac(OX_dw,OX_dw,base_gamma+5);
  safe_dirac_prod_spinspin(OX_dw,base_gamma+5,OX_dw);
  
  complex wAB;
  complex_prod_double(wAB,AB,w);
  
  //trace
  for(int ig=0;ig<16;ig++)
    {
      //no debug: comment product by wAB
      spinspin GOX_dw;
      unsafe_dirac_prod_spinspin(GOX_dw,base_gamma+ig,OX_dw);
      
      spinspin GXO_up;
      unsafe_dirac_prod_spinspin(GXO_up,base_gamma+ig,XO_up);
      
      //close the diagram and trace
      spinspin t;
      unsafe_spinspin_prod_spinspin(t,GXO_up,GOX_dw);
      complex c;
      trace_spinspin(c,t);
      
      //multiply by gluon
      complex_summ_the_prod(corr[ig],c,wAB);
      //complex_summ_the_prod_double(corr[ig],c,w);
    }
}

void compute_meson_exchange_correction_analyticallyB(corr16 *corr,quark_info qu,gluon_info gl)
{
  if(nranks>1) CRASH("works only on scalar");
  
  memset(corr,0,sizeof(corr16)*loc_vol);
  
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spinspin *temp=nissa_malloc("temp",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);

  int O=0;
  NISSA_LOC_VOL_LOOP(A)
    for(int mu_A=0;mu_A<4;mu_A++)
      NISSA_LOC_VOL_LOOP(X)
	NISSA_LOC_VOL_LOOP(B)
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
	      summ_the_exchange_contributionB(corr[X],XB,omg[mu_B],Bup_O,XA,omg[mu_A],Aup_O,AB[mu_A][mu_B],-0.25);
	      // + G S(0,B) (1-g_mu_B) S(Bup,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,B)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[Bup]);
	      summ_the_exchange_contributionB(corr[X],XB,omg[mu_B],Bup_O,XA,opg[mu_A],Adw_O,AB[mu_A][mu_B],+0.25);
	      // + G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1-g_mu_A) S(Aup,0) G(A,Bdw)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[A],glb_coord_of_loclx[B]);
	      summ_the_exchange_contributionB(corr[X],XB,opg[mu_B],Bdw_O,XA,omg[mu_A],Aup_O,AB[mu_A][mu_B],+0.25);
	      // - G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,Bdw)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,gl.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[B]);
	      summ_the_exchange_contributionB(corr[X],XB,opg[mu_B],Bdw_O,XA,opg[mu_A],Adw_O,AB[mu_A][mu_B],-0.25);
	    }
  
  nissa_free(q_prop);
  nissa_free(temp);
  nissa_free(g_prop);
}

//#endif

