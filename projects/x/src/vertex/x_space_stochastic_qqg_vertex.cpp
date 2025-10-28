#include <string.h>
#include <math.h>

#include "../../../../src/nissa.hpp"
using namespace std;

#include "../propagators/twisted_propagator.hpp"
#include "../types/types.hpp"
#include "../types/types_routines.hpp"
#include "../routines/derivatives.hpp"
#include "../routines/shift.hpp"

void stochastic_x_space_qqg_vertex_source(spinspin *q_out,spinspin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false)
{
  if(q_out==q_in) CRASH("q_in==q_out");
  
  spinspin *shift_q_in=nissa_malloc("shift_q_in",loc_vol,spinspin);
  spin1field *shift_g_in=nissa_malloc("shift_g_in",loc_vol,spin1field);
  
  //reset the output
  memset(q_out,0,sizeof(spinspin)*loc_vol);
  
  for(int mu=0;mu<4;mu++)
    {
      //shift=q_in_up_mu
      shift_spinspin_sink_up(shift_q_in,q_in,qu.bc,mu);

      NISSA_LOC_VOL_LOOP(ivol)
        {
	  spinspin temp;
	  
	  //temp=(1-gamma_mu)*q_in_up_mu
	  unsafe_dirac_prod_spinspin(temp,base_gamma+map_mu[mu],shift_q_in[ivol]);
	  spinspin_subt(temp,shift_q_in[ivol],temp);
	  
	  //q_out+=temp*g_in_mu
	  if(g_dag==false) spinspin_summ_the_complex_prod(q_out[ivol],temp,g_in[ivol][mu]);
	  else             spinspin_summ_the_complex_conj2_prod(q_out[ivol],temp,g_in[ivol][mu]);
	}
      
      //shift=q_in_dw_mu
      shift_spinspin_sink_dw(shift_q_in,q_in,qu.bc,mu);
      //shift=g_in_dw_mu_mu
      shift_spin1field_up(shift_g_in,g_in,gl.bc,mu);
      
      NISSA_LOC_VOL_LOOP(ivol)
        {
	  spinspin temp;
	  
	  //temp=(1+gamma_mu)*q_in_dw_mu
	  unsafe_dirac_prod_spinspin(temp,base_gamma+map_mu[mu],shift_q_in[ivol]);
	  spinspin_summ(temp,shift_q_in[ivol],temp);
	  
	  //q_out-=temp*g_in_dw_mu_mu
	  if(g_dag==false) spinspin_subt_the_complex_prod(q_out[ivol],temp,shift_g_in[ivol][mu]);
	  else             spinspin_subt_the_complex_conj2_prod(q_out[ivol],temp,shift_g_in[ivol][mu]);
	}
    }
  
  //put i/2
  NISSA_LOC_VOL_LOOP(ivol)
    spinspin_prodassign_idouble(q_out[ivol],0.5);
  
  set_borders_invalid(q_out);
  
  nissa_free(shift_q_in);
  nissa_free(shift_g_in);
}

void stochastic_x_space_qqg_vertex(spinspin *q_out,spinspin *q_in,quark_info qu,spin1field *g_in,gluon_info gl,bool g_dag=false)
{
  stochastic_x_space_qqg_vertex_source(q_out,q_in,qu,g_in,gl,g_dag);
  multiply_from_left_by_x_space_twisted_propagator_by_fft(q_out,q_out,qu);
}
