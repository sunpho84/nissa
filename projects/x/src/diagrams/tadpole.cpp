#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../propagators/tlSym_gluon_propagator.hpp"
#include "../propagators/twisted_propagator.hpp"
#include "../routines/shift.hpp"
#include "../types/types.hpp"

void compute_tadpole_diagram_in_mom_space(spinspin *q_tad,quark_info qu,gluon_info gl)
{
  //compute tadpole integral
  spinspin loc_tad;
  spinspin_put_to_zero(loc_tad);
  NISSA_LOC_VOL_LOOP(imom)
    {
      spin1prop g_prop;
      mom_space_tlSym_gluon_propagator_of_imom(g_prop,gl,imom);
      spinspin_summassign(loc_tad,g_prop);
    }
  spin1prop glb_tad;
  MPI_Allreduce(loc_tad,glb_tad,32,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  spinspin_prodassign_double(glb_tad,1.0/glb_vol);
  
  NISSA_LOC_VOL_LOOP(imom)
    {
      spinspin_put_to_zero(q_tad[imom]);
      for(int mu=0;mu<4;mu++)
	{
	  //compute momentum and its trigo
	  double ap_mu=M_PI*(2*glb_coord_of_loclx[imom][mu]+qu.bc[mu])/glb_size[mu];
	  double c=cos(ap_mu),s=sin(ap_mu);
	  
	  //compute weights for id and gamma
	  complex wid,wga;
	  complex_prod_double(wid,glb_tad[mu][mu],-c/2);
	  safe_complex_prod_idouble(wga,glb_tad[mu][mu],s/2);
	  
	  //add the id and gamma
	  spinspin_dirac_summ_the_prod_complex(q_tad[imom],&(base_gamma[0]),wid);
	  spinspin_dirac_summ_the_prod_complex(q_tad[imom],&(base_gamma[map_mu[mu]]),wga);
	}
    }
}

void compute_tadpole_diagram_in_x_space(spinspin *q_tad,quark_info qu,gluon_info gl)
{
  vector_reset(q_tad);
  
  //compute the propagator
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol+bord_vol,spin1prop);
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);
  NISSA_LOC_VOL_LOOP(ivol)
    if(glblx_of_loclx[ivol]!=0)
      spinspin_put_to_zero(g_prop[ivol]);
  
  spin1prop *g_prop_sh=nissa_malloc("g_prop_sh",loc_vol,spin1prop);
  for(int mu=0;mu<4;mu++)
    {
      shift_spinspin_source_up(g_prop_sh,g_prop,qu.bc,mu);
      
      NISSA_LOC_VOL_LOOP(ivol)
      {
	spinspin temp;
	unsafe_spinspin_complex_prod(temp,opg[mu],g_prop_sh[ivol][mu][mu]);
	spinspin_prodassign_double(temp,-0.25);
	spinspin_summassign(q_tad[ivol],temp);
      }

      shift_spinspin_source_dw(g_prop_sh,g_prop,qu.bc,mu);
      
      NISSA_LOC_VOL_LOOP(ivol)
      {
	spinspin temp;
	unsafe_spinspin_complex_prod(temp,omg[mu],g_prop_sh[ivol][mu][mu]);
	spinspin_prodassign_double(temp,-0.25);
	spinspin_summassign(q_tad[ivol],temp);
      }
    }
  
  nissa_free(g_prop_sh);
  nissa_free(g_prop);
}

void finish_tadpole_computation(spinspin *q_out,quark_info qu)
{
  NISSA_LOC_VOL_LOOP(imom)
  {
    spinspin q_prop,t;
    mom_space_twisted_propagator_of_imom(q_prop,qu,imom);
    
    unsafe_spinspin_prod_spinspin(t, q_prop,q_out[imom]);
    unsafe_spinspin_prod_spinspin(q_out[imom], t,q_prop);
    spinspin_prodassign_double(q_out[imom],glb_vol2);
  }

  pass_spinspin_from_mom_to_x_space(q_out,q_out,qu.bc);
}

void compute_tadpole_twisted_propagator_in_x_space(spinspin *q_out,quark_info qu,gluon_info gl)
{
  compute_tadpole_diagram_in_x_space(q_out,qu,gl);
  pass_spinspin_from_x_to_mom_space(q_out,q_out,qu.bc);
  
  finish_tadpole_computation(q_out,qu);
}

void compute_tadpole_twisted_propagator_in_mom_space(spinspin *q_out,quark_info qu,gluon_info gl)
{
  compute_tadpole_diagram_in_mom_space(q_out,qu,gl);

  finish_tadpole_computation(q_out,qu);
}
