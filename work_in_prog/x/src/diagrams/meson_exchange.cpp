#include <math.h>

#include "nissa.h"

#include "../propagators/twisted_propagator.h"
#include "../propagators/tlSym_gluon_propagator.h"
#include "../routines/fourier.h"
#include "../routines/shift.h"

#include "propagator_self_energy.h"

/*
    __A__
   /  _} \
  0  {_   X
   \__ }_/
      B
 */
 

void summ_the_exchange_contribution(corr16 corr,spinspin OB,spinspin pB,spinspin BX,spinspin XA,spinspin pA,spinspin AO,complex AB,double w)
{
  complex wAB;
  complex_prod_double(wAB,AB,w);
  
  spinspin pO,XO; //upper line
  unsafe_spinspin_spinspin_prod(pO,pA,AO);
  unsafe_spinspin_spinspin_prod(XO,XA,pO);
  
  spinspin pX,OX; //lower line
  unsafe_spinspin_spinspin_prod(pO,pA,AO);
  unsafe_spinspin_spinspin_prod(XO,XA,pO);
  
  spinspin XO_AB; //upper line with gluon and weight
  unsafe_spinspin_complex_prod(XO_AB,XO,wAB);
  
  
  master_printf("OB\n");
  print_spinspin(OB);
  master_printf("pB\n");
  print_spinspin(pB);
  master_printf("BX\n");
  print_spinspin(BX);
  master_printf("XA\n");
  print_spinspin(XA);
  master_printf("pA\n");
  print_spinspin(pA);
  master_printf("AO\n");
  print_spinspin(AO);
  master_printf("AB\n");
  master_printf("%lg %lg\n",AB[0],AB[1]);
  master_printf("w\n");
  master_printf("%lg",w);
  master_printf("XO\n");
  print_spinspin(XO);
  master_printf("XO_AB\n");
  print_spinspin(XO_AB);

  //trace
  for(int ig=0;ig<16;ig++)
    {
      spinspin GXO,GXO_AB;
      spinspin_dirac_spinspin_prod(GXO,base_gamma+ig,XO);
      spinspin_dirac_spinspin_prod(GXO_AB,base_gamma+ig,XO_AB);
      summ_the_trace_prod_spinspins(corr[ig],GXO,GXO_AB);
      if(ig==5)
	{
	  master_printf("GXO\n");
	  print_spinspin(GXO);
	  master_printf("GXO_AB\n");
	  print_spinspin(GXO_AB);
	  master_printf("GXO*GXO_AB\n");
	  spinspin test;
	  unsafe_spinspin_spinspin_prod(test,GXO,GXO_AB);
	  print_spinspin(test);
	}
    }
  master_printf("corr[5]=%lg\n",corr[5][0]);
}

void compute_meson_exchange_correction_analytically(corr16 *corr,quark_info qu,gluon_info gl)
{
  if(rank_tot>1) crash("works only on scalar");
  
  memset(corr,0,sizeof(corr16)*loc_vol);
  
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
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
	      compute_x_space_propagator_to_sink_from_source(XA,q_prop,qu.bc,glb_coord_of_loclx[X],glb_coord_of_loclx[A]);
	      compute_x_space_propagator_to_sink_from_source(Bdw_X,q_prop,qu.bc,glb_coord_of_loclx[Bdw],glb_coord_of_loclx[X]);
	      compute_x_space_propagator_to_sink_from_source(Bup_X,q_prop,qu.bc,glb_coord_of_loclx[Bup],glb_coord_of_loclx[X]);
	      compute_x_space_propagator_to_sink_from_source(OB,q_prop,qu.bc,glb_coord_of_loclx[O],glb_coord_of_loclx[B]);
	      
	      // - G S(0,B) (1-g_mu_B) S(Bup,X) G S(X,A) (1-g_mu_A) S(Aup,0) G(A,B)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,qu.bc,glb_coord_of_loclx[A],glb_coord_of_loclx[B]);
	      summ_the_exchange_contribution(corr[X],OB,nissa_omg[mu_B],Bup_X,XA,nissa_omg[mu_A],Aup_O,AB[mu_A][mu_B],-0.25);
	      // + G S(0,B) (1-g_mu_B) S(Bup,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,B)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,qu.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[B]);
	      summ_the_exchange_contribution(corr[X],OB,nissa_omg[mu_B],Bup_X,XA,nissa_opg[mu_A],Adw_O,AB[mu_A][mu_B],+0.25);
	      // + G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1-g_mu_A) S(Aup,0) G(A,B)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,qu.bc,glb_coord_of_loclx[A],glb_coord_of_loclx[Bdw]);
	      summ_the_exchange_contribution(corr[X],OB,nissa_opg[mu_B],Bdw_X,XA,nissa_omg[mu_A],Aup_O,AB[mu_A][mu_B],+0.25);
	      // - G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,B)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,qu.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[Bdw]);
	      summ_the_exchange_contribution(corr[X],OB,nissa_opg[mu_B],Bdw_X,XA,nissa_opg[mu_A],Adw_O,AB[mu_A][mu_B],-0.25);
	    }
  
  nissa_free(q_prop);
  nissa_free(g_prop);
}
