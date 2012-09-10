#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/propagators/tlSym_gluon_propagator.h"
#include "../src/diagrams/meson_exchange.h"
#include "../src/diagrams/propagator_self_energy.h"
#include "../src/types/types_routines.h"
#include "../src/routines/read_and_write.h"
#include "../src/vertex/vertex.h"
#include "../src/routines/correlations.h"
#include "../src/stochastic/stochastic_twisted_propagator.h"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.h"

//perfor 1d fft multiplying by V**ord
void transf(double *corrx,complex *corrp,int ord)
{
  //fft acts on complex, so buffering
  complex buf[glb_size[0]];
  
  //parameters for fft
  int mu=0;
  int ncompl_per_site=1;
  int sign=-1;
  int normalize=0;
  fft1d(buf,corrp,ncompl_per_site,mu,sign,normalize);
  
  //normalize
  for(int t=0;t<glb_size[0];t++)
    corrx[t]=buf[t][0]*pow(glb_vol,ord)/glb_size[0]*3;
}

void compute(double *lead,double *self,double *exch,quark_info qu,gluon_info gl)
{
  if(nissa_nranks>1) crash("works only on scalar");

  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);

  complex corrp[glb_size[0]];
  
  compute_mom_space_twisted_propagator(q_prop,qu);
  compute_mom_space_tlSym_gluon_propagator(g_prop,gl);
  
  ////////////////////////////////////////////////////////////// leading order //////////////////////////////////////////////////////
  
  //reset the corr in p space
  memset(corrp,0,glb_size[0]*sizeof(complex));
  
  coords pc={0,0,0,0};
  for(pc[0]=0;pc[0]<glb_size[0];pc[0]++)
    {
      int p=glblx_of_coord(pc);
      
      nissa_loc_vol_loop(q)
      {
	int ppq=glblx_of_comb(p,+1,q,+1);
	
	spinspin d1,d2;
	unsafe_spinspin_prod_dirac(d1,q_prop[ppq],base_gamma+5);
	unsafe_spinspin_prod_dirac(d2,q_prop[q],base_gamma+5);
	summ_the_trace_prod_spinspins(corrp[glb_coord_of_loclx[p][0]],d1,d2);
      }
    }
  
  transf(lead,corrp,1);

  
  //////////////////////////////////////////////////////////////// self //////////////////////////////////////////////////////////////
  
  //temporary we do it with external routine
  spinspin *self_prop=nissa_malloc("self_prop",loc_vol,spinspin);
  compute_self_energy_twisted_diagram_in_x_space(self_prop,qu,gl);
  pass_spinspin_from_x_to_mom_space(self_prop,self_prop,qu.bc);
  
  //put the legs
  nissa_loc_vol_loop(imom)
  {
    safe_spinspin_prod_spinspin(self_prop[imom],self_prop[imom],q_prop[imom]);
    safe_spinspin_prod_spinspin(self_prop[imom],q_prop[imom],self_prop[imom]);
    spinspin_prodassign_double(self_prop[imom],glb_vol*glb_vol);
  }
    
  //reset the corr in p space
  memset(corrp,0,glb_size[0]*sizeof(complex));
  
  for(pc[0]=0;pc[0]<glb_size[0];pc[0]++)
    {
      int p=glblx_of_coord(pc);
      
      nissa_loc_vol_loop(q)
      {
	int ppq=glblx_of_comb(p,+1,q,+1);
	
	spinspin d1,d2;
	unsafe_spinspin_prod_dirac(d1,self_prop[ppq],base_gamma+5);
	unsafe_spinspin_prod_dirac(d2,q_prop[q],base_gamma+5);
	summ_the_trace_prod_spinspins(corrp[glb_coord_of_loclx[p][0]],d1,d2);
      }
    }
  
  transf(self,corrp,1);
  
  nissa_free(self_prop);
  
  ////////////////////////////////////////////////////// self_bis ////////////////////////////////////////
  
  //check this on a small volume with previous implementation
  //if this passes we are 50% from exchange
  
  spinspin *self_bis_prop=nissa_malloc("self_bis_prop",loc_vol,spinspin);
  
  /*
  nissa_loc_vol_loop(ip)
  {
    nissa_loc_vol_loop(iq)
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
          
	  spinspin_dirac_prod_double(w[mu],&(base_gamma[nissa_map_mu[mu]]),co);
	  spinspin_dirac_summ_the_prod_idouble(w[mu],&(base_gamma[0]),-si);
	}
      
      int ir=glblx_of_coord(ir_mu);
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  {
	    spinspin temp;
	    unsafe_spinspin_prod_spinspin(temp,w[mu],q_prop[ir]);
	    safe_spinspin_prod_spinspin(temp,temp,w[nu]);
	    spinspin_summ_the_complex_prod(self_bis_prop[ip],temp,g_prop[iq][mu][nu]);
	  }
    }
    
    spinspin_prodassign_double(self_bis_prop[ip],-1);
  }
  */
  compute_self_energy_twisted_diagram_in_mom_space(self_bis_prop,qu,gl);
  
  //recomputed for check
  self_prop=nissa_malloc("self_prop",loc_vol,spinspin);
  compute_self_energy_twisted_diagram_in_x_space(self_prop,qu,gl);
  pass_spinspin_from_x_to_mom_space(self_prop,self_prop,qu.bc);

  nissa_loc_vol_loop(imom)
  {
    master_printf("-----%d\n",imom);
    print_spinspin(self_prop[imom]);
    print_spinspin(self_bis_prop[imom]);
    master_printf("%lg\n",self_prop[imom][0][0][0]-self_bis_prop[imom][0][0][0]);
  }
  nissa_free(self_prop);

  //put the legs
  nissa_loc_vol_loop(imom)
  {
    safe_spinspin_prod_spinspin(self_bis_prop[imom],self_bis_prop[imom],q_prop[imom]);
    safe_spinspin_prod_spinspin(self_bis_prop[imom],q_prop[imom],self_bis_prop[imom]);
    spinspin_prodassign_double(self_bis_prop[imom],glb_vol*glb_vol);
  }
    
  //reset the corr in p space
  memset(corrp,0,glb_size[0]*sizeof(complex));
  
  for(pc[0]=0;pc[0]<glb_size[0];pc[0]++)
    {
      int p=glblx_of_coord(pc);
      
      nissa_loc_vol_loop(q)
      {
	int ppq=glblx_of_comb(p,+1,q,+1);
	
	spinspin d1,d2;
	unsafe_spinspin_prod_dirac(d1,self_bis_prop[ppq],base_gamma+5);
	unsafe_spinspin_prod_dirac(d2,q_prop[q],base_gamma+5);
	summ_the_trace_prod_spinspins(corrp[glb_coord_of_loclx[p][0]],d1,d2);
      }
    }

  transf(exch,corrp,1);
  
  nissa_free(self_bis_prop);  
  
  nissa_free(q_prop);
  nissa_free(g_prop);
}

int main(int narg,char **arg)
{
  //basic initialization
  int T=4,L=4;
  init_nissa();
  init_grid(T,L);
  
  //quark
  double kappa=0.115;
  double mass=0.00;
  double quark_theta[4]={0,0,0,0};
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double alpha=1;
  double gluon_theta[4]={0,0,0,0};
  gluon_info gl=create_Wilson_gluon_info(alpha,gluon_theta);
  
  double lead[glb_size[0]];
  double self[glb_size[0]];
  double exch[glb_size[0]];
  
  compute(lead,self,exch,qu,gl);
  
  /////////////////////////////// print correlation /////////////////////////////

  master_printf("Lead\n");
  for(int t=0;t<glb_size[0];t++)
    master_printf("%d %16.16lg\n",t,lead[t]);
  
  master_printf("Self\n");
  for(int t=0;t<glb_size[0];t++)
    master_printf("%d %16.16lg\n",t,self[t]);
  
  master_printf("Exch\n");
  for(int t=0;t<glb_size[0];t++)
    master_printf("%d %16.16lg\n",t,exch[t]);
  
  /////////////////////////////////// close ////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
