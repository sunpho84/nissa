#include <string.h>
#include <math.h>

#include "../../../../src/base/global_variables.hpp"
#include "../../../../src/new_types/new_types_definitions.hpp"
#include "../../../../src/new_types/complex.hpp"
#include "../../../../src/new_types/dirac.hpp"
#include "../../../../src/new_types/spin.hpp"
#include "../../../../src/base/debug.hpp"
#include "../../../../src/operations/fft.hpp"
#include "../../../../src/base/vectors.hpp"

#include "../propagators/twisted_propagator.hpp"
#include "../types/types.hpp"
#include "../routines/derivatives.hpp"
#include "../routines/fourier.hpp"
#include "../routines/shift.hpp"

void stochastic_mom_space_qqg_vertex(spin *q_out,spin *q_in,quark_info qu,spin1field *g_in,gluon_info gl)
{
  if(nranks>1) CRASH("implemented only in scalar");
  if(q_out==q_in) CRASH("q_in==q_out");
  
  spin *shift_q_in=nissa_malloc("shift_q_in",loc_vol+bord_vol,spin);
  spin1field *shift_g_in=nissa_malloc("shift_g_in",loc_vol+bord_vol,spin1field);
  
  //reset the output
  memset(q_out,0,sizeof(spin)*loc_vol);
  
  for(int mu=0;mu<4;mu++)
    {
      dirac_matr gamma_mu=(mu==0)?base_gamma[4]:base_gamma[mu];
      
      //shift=q_in_up_mu
      shift_spin_dw(shift_q_in,q_in,qu.bc,mu);

      NISSA_LOC_VOL_LOOP(ivol)
        {
	  spin temp;
	  
	  //temp=(1-gamma_mu)*q_in_up_mu
	  unsafe_dirac_prod_spin(temp,&gamma_mu,shift_q_in[ivol]);
	  spin_subt(temp,shift_q_in[ivol],temp);
	  
	  //q_out+=temp*g_in_mu
	  spin_summ_the_complex_prod(q_out[ivol],temp,g_in[ivol][mu]);
	}

      //shift=q_in_dw_mu
      shift_spin_up(shift_q_in,q_in,qu.bc,mu);
      //shift=g_in_dw_mu_mu
      shift_spin1field_up(shift_g_in,g_in,gl.bc,mu);

      NISSA_LOC_VOL_LOOP(ivol)
        {
	  spin temp;
	  
	  //temp=(1+gamma_mu)*q_in_dw_mu
	  unsafe_dirac_prod_spin(temp,&gamma_mu,shift_q_in[ivol]);
	  spin_summ(temp,shift_q_in[ivol],temp);
	  
	  //q_out-=temp*g_in_dw_mu_mu
	  spin_subt_the_complex_prod(q_out[ivol],temp,shift_g_in[ivol][mu]);
	}
    }

  //put -i/2
  NISSA_LOC_VOL_LOOP(ivol)
    for(int id=0;id<4;id++)
      complex_prodassign_idouble(q_out[ivol][id],-0.5);
  
  set_borders_invalid(q_out);
  
  nissa_free(shift_q_in);
  nissa_free(shift_g_in);
}

void stochastic_x_space_qqg_vertex(spin *q_out,spin *q_in,quark_info qu,spin1field *g_in,gluon_info gl)
{
  stochastic_x_space_qqg_vertex_source(q_out,q_in,qu,g_in,gl);
  multiply_x_space_twisted_propagator_by_inverting(q_out,q_out,qu);
}

void stochastic_x_space_qqg_vertex(spinspin *q_out,spinspin *q_in,quark_info qu,spin1field *g_in,gluon_info gl)
{
  //source and temp prop
  spin *tsource=nissa_malloc("tsource",loc_vol+bord_vol,spin);
  spin *tprop=nissa_malloc("tprop",loc_vol,spin);
    
  //loop over the source index
  for(int id_so=0;id_so<4;id_so++)
    {
      //prepare the source
      NISSA_LOC_VOL_LOOP(ivol)
        for(int id_si=0;id_si<4;id_si++)
          complex_copy(tsource[ivol][id_si],q_in[ivol][id_si][id_so]);
      
      //operate on single source index
      stochastic_x_space_qqg_vertex(tprop,tsource,qu,g_in,gl);
      NISSA_LOC_VOL_LOOP(ivol)
        for(int id_si=0;id_si<4;id_si++)
          complex_copy(q_out[ivol][id_si][id_so],tprop[ivol][id_si]);
    }
  
  set_borders_invalid(q_out);
  
  nissa_free(tsource);
  nissa_free(tprop);
}
