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
  int deb=0;
  if(deb) master_printf("pA\n");
  if(deb) print_spinspin(pA);
  if(deb) master_printf("AO\n");
  if(deb) print_spinspin(AO);
  spinspin pO,XO; //upper line
  unsafe_spinspin_spinspin_prod(pO, pA,AO);
  if(deb) master_printf("pO\n");
  if(deb) print_spinspin(pO);
  if(deb) master_printf("XA\n");
  if(deb) print_spinspin(XA);
  unsafe_spinspin_spinspin_prod(XO, XA,pO);
  if(deb) master_printf("XO\n");
  if(deb) print_spinspin(XO);
  
  if(deb) master_printf("pB\n");
  if(deb) print_spinspin(pB);
  if(deb) master_printf("BX\n");
  if(deb) print_spinspin(BX);
  spinspin pX,OX; //lower line
  unsafe_spinspin_spinspin_prod(pX, pB,BX);
  if(deb) master_printf("pX\n");
  if(deb) print_spinspin(pX);
  if(deb) master_printf("OB\n");
  if(deb) print_spinspin(OB);
  unsafe_spinspin_spinspin_prod(OX, OB,pX);
  if(deb) master_printf("OX\n");
  if(deb) print_spinspin(OX);
  
  if(deb) master_printf("AB %lg %lg\n",AB[0],AB[1]);
  if(deb) master_printf("w %lg\n",w);
  complex wAB;
  complex_prod_double(wAB,AB,w);
  if(deb) master_printf("wAB %lg %lg\n",wAB[0],wAB[1]);
  
  spinspin XO_AB; //upper line with gluon and weight
  unsafe_spinspin_complex_prod(XO_AB,XO,wAB);
  if(deb) master_printf("XO_AB\n");
  if(deb) print_spinspin(XO_AB);  
  
  //trace
  for(int ig=0;ig<16;ig++)
    {
      spinspin GOX,GXO_AB;
      spinspin_dirac_spinspin_prod(GOX,base_gamma+ig,OX);
      spinspin_dirac_spinspin_prod(GXO_AB,base_gamma+ig,XO_AB);
      summ_the_trace_prod_spinspins(corr[ig],GOX,GXO_AB);
      if(0)
      if(ig==0)
	{
	  master_printf("GOX\n");
	  print_spinspin(GOX);
	  master_printf("GXO_AB\n");
	  print_spinspin(GXO_AB);
	  master_printf("GOX*GXO_AB\n");
	  spinspin test;
	  unsafe_spinspin_spinspin_prod(test,GOX,GXO_AB);
	  print_spinspin(test);
	  complex tr;
	  trace_spinspin(tr,test); 
	  master_printf("trace: %lg\n",tr[0]);
	}
    }
  //master_printf("corr[5]=%lg\n",corr[5][0]);
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
	      // + G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1-g_mu_A) S(Aup,0) G(A,Bdw)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,qu.bc,glb_coord_of_loclx[A],glb_coord_of_loclx[Bdw]);
	      summ_the_exchange_contribution(corr[X],OB,nissa_opg[mu_B],Bdw_X,XA,nissa_omg[mu_A],Aup_O,AB[mu_A][mu_B],+0.25);
	      // - G S(0,B) (1+g_mu_B) S(Bdw,X) G S(X,A) (1+g_mu_A) S(Adw,0) G(Adw,Bdw)
	      compute_x_space_propagator_to_sink_from_source(AB,g_prop,qu.bc,glb_coord_of_loclx[Adw],glb_coord_of_loclx[Bdw]);
	      summ_the_exchange_contribution(corr[X],OB,nissa_opg[mu_B],Bdw_X,XA,nissa_opg[mu_A],Adw_O,AB[mu_A][mu_B],-0.25);
	      //master_printf("corr term\n");
	    }
  
  nissa_free(q_prop);
  nissa_free(g_prop);
}
