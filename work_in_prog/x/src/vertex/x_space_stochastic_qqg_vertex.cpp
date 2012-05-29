#include <string.h>
#include <math.h>

#include "../../../../src/base/global_variables.h"
#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/new_types/complex.h"
#include "../../../../src/new_types/dirac.h"
#include "../../../../src/new_types/spin.h"
#include "../../../../src/base/debug.h"
#include "../../../../src/operations/fft.h"
#include "../../../../src/base/vectors.h"

#include "../propagators/twisted_propagator.h"
#include "../types/types.h"
#include "../routines/derivatives.h"
#include "../routines/fourier.h"
#include "../routines/shift.h"

void stochastic_x_space_qqg_vertex_source(spin *q_out,spin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false)
{
  if(q_out==q_in) crash("q_in==q_out");
  
  spin *shift_q_in=nissa_malloc("shift_q_in",loc_vol+bord_vol,spin);
  spin1field *shift_g_in=nissa_malloc("shift_g_in",loc_vol+bord_vol,spin1field);
  
  //reset the output
  memset(q_out,0,sizeof(spin)*loc_vol);
  
  for(int mu=0;mu<4;mu++)
    {
      dirac_matr gamma_mu=(mu==0)?base_gamma[4]:base_gamma[mu];
      
      //shift=q_in_up_mu
      shift_spin_dw(shift_q_in,q_in,qu.bc,mu);

      nissa_loc_vol_loop(ivol)
        {
	  spin temp;
	  
	  //temp=(1-gamma_mu)*q_in_up_mu
	  unsafe_dirac_prod_spin(temp,&gamma_mu,shift_q_in[ivol]);
	  spin_subt(temp,shift_q_in[ivol],temp);
	  
	  //q_out+=temp*g_in_mu
	  if(g_dag==false) spin_summ_the_complex_prod(q_out[ivol],temp,g_in[ivol][mu]);
	  else             spin_summ_the_complex_conj2_prod(q_out[ivol],temp,g_in[ivol][mu]);
	}

      //shift=q_in_dw_mu
      shift_spin_up(shift_q_in,q_in,qu.bc,mu);
      //shift=g_in_dw_mu_mu
      shift_spin1field_up(shift_g_in,g_in,gl.bc,mu);

      nissa_loc_vol_loop(ivol)
        {
	  spin temp;
	  
	  //temp=(1+gamma_mu)*q_in_dw_mu
	  unsafe_dirac_prod_spin(temp,&gamma_mu,shift_q_in[ivol]);
	  spin_summ(temp,shift_q_in[ivol],temp);
	  
	  //q_out-=temp*g_in_dw_mu_mu
	  if(g_dag==false) spin_subt_the_complex_prod(q_out[ivol],temp,shift_g_in[ivol][mu]);
	  else             spin_subt_the_complex_conj2_prod(q_out[ivol],temp,shift_g_in[ivol][mu]);
	}
    }

  //put -i/2
  nissa_loc_vol_loop(ivol)
    for(int id=0;id<4;id++)
      complex_prodassign_idouble(q_out[ivol][id],-0.5);
  
  set_borders_invalid(q_out);
  
  nissa_free(shift_q_in);
  nissa_free(shift_g_in);
}

void stochastic_x_space_qqg_vertex_source(spinspin *q_out,spinspin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false)
{
  //source and temp prop
  spin *tsource=nissa_malloc("tsource",loc_vol+bord_vol,spin);
  spin *tprop=nissa_malloc("tprop",loc_vol,spin);
    
  //loop over the source index
  for(int id_so=0;id_so<4;id_so++)
    {
      //prepare the source
      nissa_loc_vol_loop(ivol)
        for(int id_si=0;id_si<4;id_si++)
          complex_copy(tsource[ivol][id_si],q_in[ivol][id_si][id_so]);
      
      //operate on single source index
      stochastic_x_space_qqg_vertex_source(tprop,tsource,qu,g_in,gl,g_dag);
      nissa_loc_vol_loop(ivol)
        for(int id_si=0;id_si<4;id_si++)
          complex_copy(q_out[ivol][id_si][id_so],tprop[ivol][id_si]);
    }
  
  set_borders_invalid(q_out);
  
  nissa_free(tsource);
  nissa_free(tprop);
}

void stochastic_x_space_qqg_vertex(spin *q_out,spin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false)
{
  stochastic_x_space_qqg_vertex_source(q_out,q_in,qu,g_in,gl,g_dag);
  multiply_x_space_twisted_propagator_by_inverting(q_out,q_out,qu);
}

void stochastic_x_space_qqg_vertex(spinspin *q_out,spinspin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false)
{
  //source and temp prop
  spin *tsource=nissa_malloc("tsource",loc_vol+bord_vol,spin);
  spin *tprop=nissa_malloc("tprop",loc_vol,spin);
    
  //loop over the source index
  for(int id_so=0;id_so<4;id_so++)
    {
      //prepare the source
      nissa_loc_vol_loop(ivol)
        for(int id_si=0;id_si<4;id_si++)
          complex_copy(tsource[ivol][id_si],q_in[ivol][id_si][id_so]);
      
      //operate on single source index
      stochastic_x_space_qqg_vertex(tprop,tsource,qu,g_in,gl,g_dag);
      nissa_loc_vol_loop(ivol)
        for(int id_si=0;id_si<4;id_si++)
          complex_copy(q_out[ivol][id_si][id_so],tprop[ivol][id_si]);
    }
  
  set_borders_invalid(q_out);
  
  nissa_free(tsource);
  nissa_free(tprop);
}
