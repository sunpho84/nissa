#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/propagators/tlSym_gluon_propagator.hpp"
#include "../src/diagrams/propagator_self_energy.hpp"
#include "../src/diagrams/meson_exchange.hpp"
#include "../src/diagrams/tadpole.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/routines/read_and_write.hpp"
#include "../src/vertex/vertex.hpp"
#include "../src/routines/correlations.hpp"
#include "../src/stochastic/stochastic_twisted_propagator.hpp"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.hpp"

int igamma=0;
int REIM=0;

int comp_tree=0;
int comp_exch=0;
int comp_tad=0;
int comp_self=1;

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
    corrx[t]=buf[t][REIM]*pow(glb_vol,ord)/glb_size[0]*3;
}

void compute(double *lead,double *self,double *exch,double *tad,quark_info qu,gluon_info gl)
{
  if(nranks>1) CRASH("works only on scalar");

  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);

  complex corrp[glb_size[0]];
  
  compute_mom_space_twisted_propagator(q_prop,qu);
  compute_mom_space_tlSym_gluon_propagator(g_prop,gl);

  coords pc={0,0,0,0};
  
  ////////////////////////////////////////////////////////////// leading order //////////////////////////////////////////////////////
  
  if(comp_tree)
    {
      //reset the corr in p space
      memset(corrp,0,glb_size[0]*sizeof(complex));
      
      for(pc[0]=0;pc[0]<glb_size[0];pc[0]++)
	{
	  int p=glblx_of_coord(pc);
	  
	  NISSA_LOC_VOL_LOOP(q)
	  {
	    int ppq=glblx_of_comb(p,+1,q,+1);
	    
	    spinspin d1,d2;
	    unsafe_spinspin_prod_dirac(d1,q_prop[ppq],base_gamma+igamma);
	    unsafe_spinspin_prod_dirac(d2,q_prop[q],base_gamma+igamma);
	    summ_the_trace_prod_spinspins(corrp[glb_coord_of_loclx[p][0]],d1,d2);
	  }
	}
      
      transf(lead,corrp,1);
    }
  
  //////////////////////////////////////////////////////////////// self //////////////////////////////////////////////////////////////
  
  if(comp_self)
    {
      spinspin *self_prop=nissa_malloc("self_prop",loc_vol,spinspin);
      compute_self_energy_twisted_propagator_in_x_space(self_prop,qu,gl);
      
      corr16 *corr_self=nissa_malloc("corr",loc_vol,corr16);
      pass_spinspin_from_mom_to_x_space(q_prop,q_prop,qu.bc);
      compute_all_2pts_qdagq_correlations(corr_self,q_prop,self_prop);
      pass_spinspin_from_x_to_mom_space(q_prop,q_prop,qu.bc);

      memset(self,0,sizeof(double)*glb_size[0]);
      NISSA_LOC_VOL_LOOP(ivol)
	self[glb_coord_of_loclx[ivol][0]]+=corr_self[ivol][igamma][REIM]*3;
      
      nissa_free(self_prop);
      nissa_free(corr_self);
    }
  
  ////////////////////////////////////////////////////// self_bis ////////////////////////////////////////
  
  //#define bis
#ifdef bis
    //check this on a small volume with previous implementation
    //if this passes we are 50% from exchange
    if(comp_self)
      {
	spinspin *self_prop_bis=nissa_malloc("self_prop_bis",loc_vol,spinspin);
	
	compute_self_energy_twisted_diagram_in_mom_space(self_prop_bis,qu,gl);
	
	vector_reset(self_prop_bis);
	NISSA_LOC_VOL_LOOP(ip)
	{
	  NISSA_LOC_VOL_LOOP(iq)
	  {
	    spinspin w[4];
	    
	    //compute the weights and find r
	    for(int mu=0;mu<4;mu++)
	      {
		double p_mu=M_PI*(2*glb_coord_of_loclx[ip][mu]+qu.bc[mu])/glb_size[mu];
		double q_mu=M_PI*(2*glb_coord_of_loclx[iq][mu]+gl.bc[mu])/glb_size[mu];
		double r_mu=p_mu-q_mu;
		
		mom_space_qq_vertex_function(w[mu],p_mu,r_mu,qu,mu);
	      }
	    
	    int ir=glblx_of_diff(ip,iq);
	    for(int mu=0;mu<4;mu++)
	      for(int nu=0;nu<4;nu++)
		{
		  //multiply by the vertices
		  spinspin temp;
		  unsafe_spinspin_prod_spinspin(temp,w[mu],q_prop[ir]);
		  safe_spinspin_prod_spinspin(temp,temp,w[nu]);
		  
		  //add the gluon line
		  spinspin_summ_the_complex_prod(self_prop_bis[ip],temp,g_prop[iq][mu][nu]);
		}
	  }
	}
	
	//put the legs
	NISSA_LOC_VOL_LOOP(imom)
	{
	  safe_spinspin_prod_spinspin(self_prop_bis[imom],self_prop_bis[imom],q_prop[imom]);
	  safe_spinspin_prod_spinspin(self_prop_bis[imom],q_prop[imom],self_prop_bis[imom]);
	}
	
	//reset the corr in p space
	memset(corrp,0,glb_size[0]*sizeof(complex));
	
	for(pc[0]=0;pc[0]<glb_size[0];pc[0]++)
	  {
	    int p=glblx_of_coord(pc);
	    
	    NISSA_LOC_VOL_LOOP(q)
	    {
	      int ppq=glblx_of_comb(p,+1,q,+1);
	      
	      spinspin d1,d2;
	      unsafe_spinspin_prod_dirac(d1,self_prop_bis[ppq],base_gamma+igamma);
	      unsafe_spinspin_prod_dirac(d2,q_prop[q],base_gamma+igamma);
	      summ_the_trace_prod_spinspins(corrp[glb_coord_of_loclx[p][0]],d1,d2);
	    }
	  }
	
	transf(exch,corrp,1);
	
	nissa_free(self_prop_bis);  
      }

#endif

  //////////////////////////////////////////////////// true exchange ////////////////////////////////////////////
  
    if(comp_exch)
      {
#define true_exch
#ifdef true_exch
	
	for(pc[0]=0;pc[0]<=glb_size[0]/2;pc[0]++)
	  compute_meson_exchange_correction_of_mom_gamma(corrp[pc[0]],qu,gl,glblx_of_coord(pc),igamma);
	
	//symmetrize
	for(int t=1;t<glb_size[0]/2;t++)
	  {
	    int dest=(glb_size[0]-t)%glb_size[0];
	    complex_copy(corrp[dest],corrp[t]);
	  }
	
	transf(exch,corrp,1);
	
#else
	
	corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
	compute_meson_exchange_correction_analyticallyA(corr,qu,gl);
	
	memset(exch,0,sizeof(double)*glb_size[0]);
	NISSA_LOC_VOL_LOOP(ivol)
	  exch[glb_coord_of_loclx[ivol][0]]+=corr[ivol][igamma][0];
	
	nissa_free(corr);
	
#endif
      }
    
  //////////////////////////////////////////////////// tadpole ////////////////////////////////////////////
  
    if(comp_tad)
      {
	MASTER_PRINTF("Computing tadpole corrections\n");
	
	//compute tadpole propagator
	spinspin *tad_prop=nissa_malloc("tad_prop",loc_vol,spinspin);
	compute_tadpole_twisted_propagator_in_x_space(tad_prop,qu,gl);
	
	pass_spinspin_from_mom_to_x_space(q_prop,q_prop,qu.bc);
	
	//compute correlators and write
	corr16 *corr_tad=nissa_malloc("corr",loc_vol,corr16);
	compute_all_2pts_qdagq_correlations(corr_tad,q_prop,tad_prop);
	
	pass_spinspin_from_x_to_mom_space(q_prop,q_prop,qu.bc);
	
	memset(tad,0,sizeof(double)*glb_size[0]);
	NISSA_LOC_VOL_LOOP(ivol)
	  tad[glb_coord_of_loclx[ivol][0]]+=corr_tad[ivol][igamma][REIM]*3;
	
	//free
	nissa_free(corr_tad);
	nissa_free(tad_prop);
      }

    nissa_free(q_prop);
    nissa_free(g_prop);
}

void exch_stoch_comp(double *exch,quark_info qu,gluon_info gl)
{
  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  compute_meson_exchange_correction_stochastically(corr,qu,gl);

  memset(exch,0,sizeof(double)*glb_size[0]);
  
  NISSA_LOC_VOL_LOOP(ivol)
    exch[glb_coord_of_loclx[ivol][0]]+=corr[ivol][igamma][0];
  
  nissa_free(corr);
}

int main(int narg,char **arg)
{
  //basic initialization
  int T=64,L=32;
  init_nissa(narg,arg);
  init_grid(T,L);
  start_loc_rnd_gen(100);  

  //quark
  double kappa=0.125;
  double mass=0.00;
  double quark_theta[4]={1,0,0,0};
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double alpha=1;
  double gluon_theta[4]={0,0,0,0};
  gluon_info gl=create_Wilson_gluon_info(alpha,gluon_theta);
  
  double lead[glb_size[0]];
  double self[glb_size[0]];
  double exch[glb_size[0]];
  double tad[glb_size[0]];
  
  compute(lead,self,exch,tad,qu,gl);
  
  /////////////////////////////// print correlation /////////////////////////////
  
  if(comp_tree)
    {
      MASTER_PRINTF("Lead\n");
      for(int t=0;t<glb_size[0];t++)
	MASTER_PRINTF("%d %16.16lg\n",t,lead[t]);
    }
  
  if(comp_self)
    {
      MASTER_PRINTF("Self\n");
      for(int t=0;t<glb_size[0];t++)
	MASTER_PRINTF("%d %16.16lg\n",t,self[t]);
    }

  if(comp_exch)
    {
      MASTER_PRINTF("Exch\n");
      for(int t=0;t<glb_size[0];t++)
	MASTER_PRINTF("%d %16.16lg\n",t,exch[t]);
    }

  if(comp_tad)
    {
      MASTER_PRINTF("Tad\n");
      for(int t=0;t<glb_size[0];t++)
	MASTER_PRINTF("%d %16.16lg\n",t,tad[t]);
    }
  
  double stexch_sum[glb_size[0]];
  double stexch_sum2[glb_size[0]];
  memset(stexch_sum,0,sizeof(double)*glb_size[0]);
  memset(stexch_sum2,0,sizeof(double)*glb_size[0]);
  int n=100;
  for(int i=1;i<=n;i++)
    {
      //compute stocashtically
      exch_stoch_comp(exch,qu,gl);
      
      //add to the sum and to sum2, print ave
      if(i==n) MASTER_PRINTF("Stoch_exch after %d\n",i);
      for(int t=0;t<glb_size[0];t++)
	{
	  stexch_sum[t]+=exch[t];
	  stexch_sum2[t]+=exch[t]*exch[t];
	  
	  double ave=stexch_sum[t]/i;
	  double err=stexch_sum2[t]/i;
	  err=sqrt((err-ave*ave)/(i-1));
	  
	  if(i==n) MASTER_PRINTF("%d %lg %lg\n",t,ave*3,err*3);
	}
    }
  
  /////////////////////////////////// close ////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
