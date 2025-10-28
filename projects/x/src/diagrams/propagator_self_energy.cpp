#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../propagators/twisted_propagator.hpp"
#include "../propagators/tlSym_gluon_propagator.hpp"
#include "../types/types_routines.hpp"
#include "../routines/shift.hpp"

void summ_the_contribution_of_self_energy_twisted_diagram_in_x_space(spinspin *q_out,spinspin *osi,spinspin *q,spin1prop *g,int nu,int mu,spinspin *oso,double weight)
{
  NISSA_LOC_VOL_LOOP(ivol)
    {
      spinspin t;
      unsafe_spinspin_prod_spinspin(t,q[ivol],oso[mu]);
      safe_spinspin_prod_spinspin(t,osi[nu],t);
      spinspin_prod_double(t,t,weight);
      spinspin_summ_the_complex_prod(q_out[ivol],t,g[ivol][nu][mu]);
    }
}

void compute_self_energy_twisted_diagram_in_x_space(spinspin *q_out,spinspin *q_prop,quark_info qu,spin1prop *g_prop,gluon_info gl)
{
  spin1prop *g_prop_sh=nissa_malloc("g_prop_sh",loc_vol,spin1prop);
  spinspin *q_prop_sh=nissa_malloc("q_prop_sh",loc_vol,spinspin);
  
  //reset the corr
  memset(q_out,0,sizeof(spinspin)*loc_vol);
  
  //loop on sink and source lorentz index
  for(int nu=0;nu<4;nu++)
    for(int mu=0;mu<4;mu++)
      {
	//  term 1: +(1-gnu) G(x,-mu) S(x+nu,-mu) (1-gmu)
	//memcpy(g_prop_sh,g_prop,sizeof(spin1prop)*loc_vol);
	shift_spinspin_source_dw(g_prop_sh,g_prop,gl.bc,mu);
	shift_spinspin_sink_up(q_prop_sh,q_prop,qu.bc,nu);
	shift_spinspin_source_dw(q_prop_sh,q_prop_sh,qu.bc,mu);
        summ_the_contribution_of_self_energy_twisted_diagram_in_x_space(q_out,omg,q_prop_sh,g_prop_sh,nu,mu,omg,-0.25);
	
	//  term 2: -(1+gnu) G(x-nu,-mu) S(x-nu,-mu) (1-gmu)
	shift_spinspin_sink_dw(g_prop_sh,g_prop,gl.bc,nu);
	shift_spinspin_source_dw(g_prop_sh,g_prop_sh,gl.bc,mu);
	shift_spinspin_sink_dw(q_prop_sh,q_prop,qu.bc,nu);
	shift_spinspin_source_dw(q_prop_sh,q_prop_sh,qu.bc,mu);
        summ_the_contribution_of_self_energy_twisted_diagram_in_x_space(q_out,opg,q_prop_sh,g_prop_sh,nu,mu,omg,+0.25);
	
	//  term 3: -(1-gnu) G(x,0)      S(x+nu,mu) (1+gmu)
	//memcpy(g_prop_sh,g_prop,sizeof(spin1prop)*loc_vol);
	//memcpy(g_prop_sh,g_prop_sh,sizeof(spin1prop)*loc_vol);
	shift_spinspin_sink_up(q_prop_sh,q_prop,qu.bc,nu);
	shift_spinspin_source_up(q_prop_sh,q_prop_sh,qu.bc,mu);
        summ_the_contribution_of_self_energy_twisted_diagram_in_x_space(q_out,omg,q_prop_sh,g_prop,nu,mu,opg,+0.25);
	
	//  term 4: +(1+gnu) G(x-nu,0)   S(x-nu,mu) (1+gmu)
	shift_spinspin_sink_dw(g_prop_sh,g_prop,gl.bc,nu);
	//memcpy(g_prop_sh,g_prop_sh,sizeof(spin1prop)*loc_vol);
	shift_spinspin_sink_dw(q_prop_sh,q_prop,qu.bc,nu);
	shift_spinspin_source_up(q_prop_sh,q_prop_sh,qu.bc,mu);
        summ_the_contribution_of_self_energy_twisted_diagram_in_x_space(q_out,opg,q_prop_sh,g_prop_sh,nu,mu,opg,-0.25);
      }
  
  nissa_free(q_prop_sh);
  nissa_free(g_prop_sh);
}
  
void compute_self_energy_twisted_diagram_in_x_space(spinspin *q_out,quark_info qu,gluon_info gl)
{
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);
  
  compute_self_energy_twisted_diagram_in_x_space(q_out,q_prop,qu,g_prop,gl);
  
  nissa_free(q_prop);
  nissa_free(g_prop);
}

void compute_self_energy_twisted_diagram_in_mom_space(spinspin *q_out,spinspin *q_prop,quark_info qu,spin1prop *g_prop,gluon_info gl)
{
  if(nranks>1) CRASH("implemented only in scalar");
  
  memset(q_out,0,sizeof(spinspin)*loc_vol);
  
  bool zmp=1;
  for(int mu=0;mu<4;mu++) zmp&=(gl.bc[mu]==0);

  if(zmp)
    NISSA_LOC_VOL_LOOP(ip)
    {
      NISSA_LOC_VOL_LOOP(iq)
      {
	spinspin w[4];

	coords ir_mu;
	//compute the weights and find r                                                                                                                                                  
	for(int mu=0;mu<4;mu++)
	  {
	    //r=p-q                                                                                                                                                                       
	    ir_mu[mu]=(glb_coord_of_loclx[ip][mu]-glb_coord_of_loclx[iq][mu]+glb_size[mu])%glb_size[mu];
	    
	    double p_mu=M_PI*(2*glb_coord_of_loclx[ip][mu]+qu.bc[mu])/glb_size[mu];
	    double q_mu=M_PI*(2*glb_coord_of_loclx[iq][mu]+gl.bc[mu])/glb_size[mu];
	    
	    double s_mu=(2*p_mu-q_mu)/2;
	    double co=cos(s_mu),si=sin(s_mu);
	    
	    spinspin_dirac_prod_double(w[mu],&(base_gamma[map_mu[mu]]),co);
	    spinspin_dirac_summ_the_prod_idouble(w[mu],&(base_gamma[0]),-si);
	  }
	
	int ir=glblx_of_coord(ir_mu);
	for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++)
	    {
	      spinspin temp;
	      unsafe_spinspin_prod_spinspin(temp,w[mu],q_prop[ir]);
	      safe_spinspin_prod_spinspin(temp,temp,w[nu]);
	      spinspin_summ_the_complex_prod(q_out[ip],temp,g_prop[iq][mu][nu]);
	    }
      }
      
      spinspin_prodassign_double(q_out[ip],-1);
    }
  else
    CRASH("Non periodic boundary conditions not implemented yet!");  
}

void compute_self_energy_twisted_diagram_in_mom_space(spinspin *q_out,quark_info qu,gluon_info gl)
{
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol+bord_vol,spin1prop);
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol+bord_vol,spinspin);
  
  compute_mom_space_twisted_propagator(q_prop,qu);
  compute_mom_space_tlSym_gluon_propagator(g_prop,gl);
  
  compute_self_energy_twisted_diagram_in_mom_space(q_out,q_prop,qu,g_prop,gl);
  
  nissa_free(q_prop);
  nissa_free(g_prop);
}

void compute_self_energy_twisted_propagator_in_x_space(spinspin *q_out,quark_info qu,gluon_info gl)
{
  compute_self_energy_twisted_diagram_in_x_space(q_out,qu,gl);
  pass_spinspin_from_x_to_mom_space(q_out,q_out,qu.bc);

  NISSA_LOC_VOL_LOOP(imom)
  {
    spinspin s,t;
    mom_space_twisted_propagator_of_imom(s,qu,imom);
    unsafe_spinspin_prod_spinspin(t,q_out[imom],s);
    unsafe_spinspin_prod_spinspin(q_out[imom],s,t);
    spinspin_prodassign_double(q_out[imom],glb_vol2);
  }
  
  pass_spinspin_from_mom_to_x_space(q_out,q_out,qu.bc);
}


void compute_self_energy_twisted_propagator_in_x_space_tough_way(spinspin *q_out,quark_info qu,gluon_info gl)
{
  spinspin *sigma2=nissa_malloc("sigma2",loc_vol+bord_vol,spinspin);
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_self_energy_twisted_diagram_in_x_space(sigma2,qu,gl);
  
  memset(q_out,0,sizeof(spinspin)*loc_vol);
  
  //iout -- A -- B -- 0
  //     qa   si   qb
  
  NISSA_LOC_VOL_LOOP(iout)
  {
    printf("%d\n",iout);
    NISSA_LOC_VOL_LOOP(A)
      NISSA_LOC_VOL_LOOP(B)
        {
	  spinspin qa,si,qb;
	  compute_x_space_propagator_to_sink_from_source(qa,q_prop,qu.bc,glb_coord_of_loclx[iout],glb_coord_of_loclx[A]);
	  compute_x_space_propagator_to_sink_from_source(si,sigma2,qu.bc,glb_coord_of_loclx[A],glb_coord_of_loclx[B]);
	  compute_x_space_propagator_to_sink_from_source(qb,q_prop,qu.bc,glb_coord_of_loclx[B],glb_coord_of_loclx[0]);
	  
	  spinspin siqb;
	  unsafe_spinspin_prod_spinspin(siqb,si,qb);
	  spinspin_summ_the_spinspin_prod(q_out[iout],qa,siqb);
	}
  }
  nissa_free(q_prop);
  nissa_free(sigma2);
}

