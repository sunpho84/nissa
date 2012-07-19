#include <math.h>

#include "nissa.h"

#include "../types/types.h"
#include "../propagators/tlSym_gluon_propagator.h"
#include "../routines/shift.h"
#include "../routines/fourier.h"
#include "../propagators/twisted_propagator.h"

void compute_tadpole_diagram_in_mom_space(spinspin *q_tad,quark_info qu,gluon_info gl)
{
  //compute tadpole integral
  spinspin loc_tad;
  spinspin_put_to_zero(loc_tad);
  nissa_loc_vol_loop(imom)
    {
      spin1prop g_prop;
      mom_space_tlSym_gluon_propagator_of_imom(g_prop,gl,imom);
      spinspin_summassign(loc_tad,g_prop);
    }
  spin1prop glb_tad;
  MPI_Allreduce(loc_tad,glb_tad,32,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  spinspin_prodassign_double(glb_tad,1.0/glb_vol);
  
  nissa_loc_vol_loop(imom)
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
	  spinspin_dirac_summ_the_prod_complex(q_tad[imom],&(base_gamma[nissa_map_mu[mu]]),wga);
	}
    }
}

void compute_tadpole_diagram_in_x_space(spinspin *q_tad,spin1prop g_prop)
{
  memset(q_tad,0,sizeof(spinspin)*loc_vol);
  
  for(int mu=0;mu<4;mu++)
    {
      coords x_up,x_dw;
      memset(x_dw,0,sizeof(coords));
      memset(x_up,0,sizeof(coords));
      x_up[mu]=1;
      x_dw[mu]=glb_size[mu]-1;
      
      int dl,dr;
      int ul,ur;
      get_loclx_and_rank_of_coord(&dl,&dr,x_dw);
      get_loclx_and_rank_of_coord(&ul,&ur,x_up);
      
      if(rank==dr)
	{
	  unsafe_spinspin_complex_prod(q_tad[dl],nissa_omg[mu],g_prop[mu][mu]);
	  spinspin_prodassign_double(q_tad[dl],-0.25);
	}
      
      if(rank==ur)
	{
	  unsafe_spinspin_complex_prod(q_tad[ul],nissa_opg[mu],g_prop[mu][mu]);
	  spinspin_prodassign_double(q_tad[ul],-0.25);
	}
    }
}

void compute_tadpole_diagram_in_x_space(spinspin *q_tad,gluon_info gl)
{
  //compute the propagator
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);
  
  //compute propagator in the origin
  spin1prop g_prop_or;
  coords origin={0,0,0,0};
  compute_x_space_propagator_to_sink_from_source(g_prop_or,g_prop,gl.bc,origin,origin);
  
  compute_tadpole_diagram_in_x_space(q_tad,g_prop_or);
  
  nissa_free(g_prop);
}

void compute_tadpole_twisted_propagator_in_x_space(spinspin *q_out,quark_info qu,gluon_info gl)
{
  compute_tadpole_diagram_in_x_space(q_out,gl);
  pass_spinspin_from_x_to_mom_space(q_out,q_out,qu.bc);

  nissa_loc_vol_loop(imom)
  {
    spinspin s,t;
    mom_space_twisted_propagator_of_imom(s,qu,imom);
    unsafe_spinspin_prod_spinspin(t, s,q_out[imom]);
    unsafe_spinspin_prod_spinspin(q_out[imom], t,s);
    spinspin_prodassign_double(q_out[imom],glb_vol2);
  }

  pass_spinspin_from_mom_to_x_space(q_out,q_out,qu.bc);
}
