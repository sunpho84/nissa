#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/propagators/tlSym_gluon_propagator.h"
#include "../src/diagrams/meson_exchange.h"
#include "../src/types/types_routines.h"
#include "../src/routines/read_and_write.h"
#include "../src/routines/correlations.h"
#include "../src/stochastic/stochastic_twisted_propagator.h"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.h"


void compute(double *corr,quark_info qu,gluon_info gl)
{
  if(nissa_nranks>1) crash("works only on scalar");

  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);

  compute_mom_space_twisted_propagator(q_prop,qu);
  compute_mom_space_tlSym_gluon_propagator(g_prop,gl);

  //reset the corr in p space
  complex corrp[glb_size[0]];
  memset(corrp,0,glb_size[0]*sizeof(complex));
  
  /*
  coords pc={0,0,0,0};
  for(pc[0]=0;pc[0]<glb_size[0];pc[0]++)
    {
      int p=glblx_of_coord(pc);

      nissa_loc_vol_loop(q)
        {
	  int pmq=site_comb(p,+1,q,+1);
	  
	  spinspin d1,d2;
	  unsafe_spinspin_prod_dirac(d1,q_prop[pmq],base_gamma+5);
	  unsafe_spinspin_prod_dirac(d2,q_prop[q],base_gamma+5);
	  summ_the_trace_prod_spinspins(corrp[glb_coord_of_loclx[p][0]],d1,d2);
	}
    }
  */
  
  coords pc={0,0,0,0};
  for(pc[0]=0;pc[0]<glb_size[0];pc[0]++)
    {
      int p=glblx_of_coord(pc);
      
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  {
	    master_printf("part %d/%d\n",(pc[0]*4+mu)*4+nu,4*4*glb_size[0]);
	    
	    nissa_loc_vol_loop(q)
	      nissa_loc_vol_loop(r)
	      {
		int ppq=site_comb(p,+1,q,+1);
		int ppr=site_comb(p,+1,r,+1);
		int rmq=site_comb(r,+1,q,-1);
		
		spinspin vmu,vnu;
		mom_space_qq_vertex_function(vmu,site_comb(ppq,+1,ppr,+1),qu,mu);
		mom_space_qq_vertex_function(vnu,site_comb(q,+1,r,+1),qu,nu);
		
		spinspin up;
		unsafe_spinspin_prod_spinspin(up,vmu,q_prop[ppq]);
		safe_spinspin_prod_spinspin(up,q_prop[ppr],up);
		
		spinspin dw;
		unsafe_spinspin_prod_spinspin(dw,vnu,q_prop[r]);
		safe_spinspin_prod_spinspin(dw,q_prop[q],dw);
		
		spinspin up_w;
		unsafe_spinspin_complex_prod(up_w,up,g_prop[rmq][mu][nu]);
		
		int ig=0;
		
		spinspin gup_w,gupg_w;
		unsafe_dirac_prod_spinspin(gup_w,base_gamma+ig,up_w);
		safe_spinspin_prod_dirac(gupg_w,gup_w,base_gamma+ig);
		
		complex temp;
		summ_the_trace_prod_spinspins(corrp[pc[0]],gupg_w,dw);
	      }
	  }
    }
  
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
    corr[t]=buf[t][0]*pow(glb_vol,1)/glb_size[0]*3;
  
  nissa_free(q_prop);
  nissa_free(g_prop);
}

int main(int narg,char **arg)
{
  //basic initialization
  int T=8,L=8;
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
  
  double *corr=nissa_malloc("corr",glb_size[0],double);  
  
  compute(corr,qu,gl);
  
  /////////////////////////////// print correlation /////////////////////////////

  FILE *fout=stdout;//open_file("/tmp/self.dat","w");
  for(int t=0;t<glb_size[0];t++)
    master_fprintf(fout,"%d %16.16lg\n",t,glb_reduce_double(corr[t]));  
  if(rank==0) fclose(fout);
  
  /////////////////////////////////// close ////////////////////////////////////
  
  nissa_free(corr);

  close_nissa();
  
  return 0;
}
