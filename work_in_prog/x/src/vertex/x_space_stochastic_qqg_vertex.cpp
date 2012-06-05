#include <string.h>
#include <math.h>

#include "../../../../src/nissa.h"

#include "../propagators/twisted_propagator.h"
#include "../types/types.h"
#include "../types/types_routines.h"
#include "../routines/derivatives.h"
#include "../routines/fourier.h"
#include "../routines/shift.h"

void stochastic_x_space_qqg_vertex_source(spinspin *q_out,spinspin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false)
{
  if(q_out==q_in) crash("q_in==q_out");
  
  spinspin *shift_q_in=nissa_malloc("shift_q_in",loc_vol+bord_vol,spinspin);
  spin1field *shift_g_in=nissa_malloc("shift_g_in",loc_vol+bord_vol,spin1field);
  
  //reset the output
  memset(q_out,0,sizeof(spinspin)*loc_vol);
  
  for(int mu=0;mu<4;mu++)
    {
      dirac_matr gamma_mu=(mu==0)?base_gamma[4]:base_gamma[mu];
      
      //shift=q_in_dw_mu
      shift_spinspin_sink_dw(shift_q_in,q_in,qu.bc,mu);
      //shift=g_in_dw_mu_mu
      shift_spin1field_up(shift_g_in,g_in,gl.bc,mu);

      nissa_loc_vol_loop(ivol)
        {
	  spinspin temp;
	  
	  //temp=(1-gamma_mu)*q_in_dw_mu
	  unsafe_dirac_prod_spinspin(temp,&gamma_mu,shift_q_in[ivol]);
	  spinspin_subt(temp,shift_q_in[ivol],temp);
	  
	  //q_out+=temp*g_in_dw_mu_mu
	  if(g_dag==false) spinspin_summ_the_complex_prod(q_out[ivol],temp,shift_g_in[ivol][mu]);
	  else             spinspin_summ_the_complex_conj2_prod(q_out[ivol],temp,shift_g_in[ivol][mu]);
	}
      
      //shift=q_in_up_mu
      shift_spinspin_sink_up(shift_q_in,q_in,qu.bc,mu);
      
      nissa_loc_vol_loop(ivol)
        {
	  spinspin temp;
	  
	  //temp=(1+gamma_mu)*q_in_up_mu
	  unsafe_dirac_prod_spinspin(temp,&gamma_mu,shift_q_in[ivol]);
	  spinspin_summ(temp,shift_q_in[ivol],temp);
	  
	  //q_out-=temp*g_in_mu
	  if(g_dag==false) spinspin_subt_the_complex_prod(q_out[ivol],temp,g_in[ivol][mu]);
	  else             spinspin_subt_the_complex_conj2_prod(q_out[ivol],temp,g_in[ivol][mu]);
	}
    }
  
  //put i/2
  nissa_loc_vol_loop(ivol)
    spinspin_prodassign_idouble(q_out[ivol],0.5);
  
  set_borders_invalid(q_out);
  
  nissa_free(shift_q_in);
  nissa_free(shift_g_in);
}

/*void stochastic_x_space_qqg_vertex_source(spinspin *q_out,spinspin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false)
{
  //source and temp prop
  spin *tsource=nissa_malloc("tsource",loc_vol+bord_vol,spin);
  spin *tprop=nissa_malloc("tprop",loc_vol,spin);
    
  //loop over the source index
  for(int id_so=0;id_so<4;id_so++)
    {
      //operate on single source index
      get_spin_from_spinspin(tsource,q_in,id_so);
      stochastic_x_space_qqg_vertex_source(tprop,tsource,qu,g_in,gl,g_dag);
      put_spin_into_spinspin(q_out,tprop,id_so);
    }
  
  set_borders_invalid(q_out);
  
  nissa_free(tsource);
  nissa_free(tprop);
}
void stochastic_x_space_qqg_vertex(spin *q_out,spin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false)
{
  stochastic_x_space_qqg_vertex_source(q_out,q_in,qu,g_in,gl,g_dag);
  multiply_x_space_twisted_propagator_by_fft(q_out,q_out,qu);
}
*/

void stochastic_x_space_qqg_vertex(spinspin *q_out,spinspin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false)
{
  stochastic_x_space_qqg_vertex_source(q_out,q_in,qu,g_in,gl,g_dag);
  multiply_x_space_twisted_propagator_by_fft(q_out,q_out,qu);
}
