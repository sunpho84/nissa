#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../propagators/twisted_propagator.hpp"
#include "../propagators/tlSym_gluon_propagator.hpp"
#include "../routines/correlations.hpp"
#include "../vertex/vertex.hpp"
#include "../stochastic/stochastic_twisted_propagator.hpp"
#include "../stochastic/stochastic_tlSym_gluon_propagator.hpp"

/*
     ___/____A____/___
    /   \   {_    \   \
   /         _}        \
  X         {_          O
   \          }        /
    \___\____{_____\__/
        /    B     /
 */

void compute_meson_exchange_correction_of_mom_gamma(complex out,quark_info qu,gluon_info gl,int ip,int ig)
{
  //compute props
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  compute_mom_space_twisted_propagator(q_prop,qu);
  compute_mom_space_tlSym_gluon_propagator(g_prop,gl);
  
  //fill the table of momentum
  momentum_t *mom=nissa_malloc("mom",loc_vol,momentum_t);
  NISSA_LOC_VOL_LOOP(imom)
    for(int mu=0;mu<4;mu++)
      mom[imom][mu]=2*M_PI*glb_coord_of_loclx[imom][mu]/glb_size[mu];

  complex_put_to_zero(out);

  spinspin *O_line=nissa_malloc("O_line",loc_vol,spinspin);
  spinspin *X_line=nissa_malloc("X_line",loc_vol,spinspin);

  NISSA_LOC_VOL_LOOP(iq)
    {
      int ippq=glblx_of_summ(ip,iq);
      
      //compute the quark prop entering the origin and going out
      unsafe_dirac_prod_spinspin(O_line[iq],base_gamma+ig,q_prop[iq]);
      safe_spinspin_prod_spinspin(O_line[iq],q_prop[ippq],O_line[iq]);
      
      //and the quark entering X and going out
      unsafe_dirac_prod_spinspin(X_line[iq],base_gamma+ig,q_prop[ippq]);
      safe_spinspin_prod_spinspin(X_line[iq],q_prop[iq],X_line[iq]);
    }
  
  NISSA_LOC_VOL_LOOP(iq)
  {
    //compute p+q
    momentum_t ppq;
    for(int mu=0;mu<4;mu++)
      ppq[mu]=mom[ip][mu]+mom[iq][mu];
    
    spinspin q_line;
    spinspin_put_to_zero(q_line);
    
    NISSA_LOC_VOL_LOOP(ir)
    {
      //compute p+r and r-q
      //int ippr=glblx_of_summ(ip,ir);
      int irmq=glblx_of_diff(ir,iq);
      momentum_t ppr;
      //momentum_t rmq;
      for(int mu=0;mu<4;mu++)
	{
	  ppr[mu]=mom[ip][mu]+mom[ir][mu];
	  //rmq[mu]=mom[ir][mu]-mom[iq][mu];
	}
      
      //compute the vertex functions
      spinspin w_up[4],w_dw[4];
      for(int mu=0;mu<4;mu++)
	{
	  mom_space_qq_vertex_function(w_up[mu],ppq[mu],ppr[mu],qu,mu);
	  mom_space_qq_vertex_function(w_dw[mu],mom[ir][mu],mom[iq][mu],qu,mu);
	}
      
      //integrate the gluon line
      for(int mu=0;mu<4;mu++)
	{
	  spinspin temp_mu;
	  
	  unsafe_spinspin_prod_spinspin(temp_mu,X_line[ir],w_up[mu]);
	  for(int nu=0;nu<4;nu++)
	    if(gl.alpha!=1 || nu==mu)
	      {
		spinspin temp;
		unsafe_spinspin_prod_spinspin(temp,w_dw[nu],temp_mu);
		spinspin_summ_the_complex_prod(q_line,temp,g_prop[irmq][mu][nu]);
	      }
	}
    }
    
    summ_the_trace_prod_spinspins(out,q_line,O_line[iq]);
  }

  nissa_free(mom);
  
  nissa_free(O_line);
  nissa_free(X_line);
  
  nissa_free(q_prop);
  nissa_free(g_prop);
}

void compute_meson_exchange_correction_stochastically(corr16 *corr,spinspin *q_prop,quark_info qu,gluon_info gl)
{
  spin1field *phi=nissa_malloc("phi",loc_vol+bord_vol,spin1field);
  spin1field *eta=nissa_malloc("eta",loc_vol+bord_vol,spin1field);
  generate_stochastic_source_and_tlSym_gluon_propagator(phi,eta,gl);
  spinspin *prop_phi=nissa_malloc("prop_phi",loc_vol,spinspin);
  spinspin *prop_eta=nissa_malloc("prop_eta",loc_vol,spinspin);

  generate_stochastic_A_twisted_propagator(prop_phi,q_prop,qu,phi,gl);
  generate_stochastic_A_twisted_propagator(prop_eta,q_prop,qu,eta,gl);
  
  compute_all_2pts_qdagq_correlations(corr,prop_phi,prop_eta);
  
  nissa_free(phi);
  nissa_free(eta);
  nissa_free(prop_phi);
  nissa_free(prop_eta);
}

//one source
void compute_meson_exchange_correction_stochastically(corr16 *corr,quark_info qu,gluon_info gl)
{
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  
  compute_meson_exchange_correction_stochastically(corr,q_prop,qu,gl);
  nissa_free(q_prop);
}

//multiple sources
void compute_meson_exchange_correction_stochastically(corr16 *zm_ave,corr16 *zm_err,corr16 *ave,corr16 *err,quark_info qu,gluon_info gl,int n)
{
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  corr16 *zm_corr=nissa_malloc("zm_corr",glb_size[0],corr16);
  
  //reset the summs
  if(ave!=NULL) memset(ave,0,sizeof(corr16)*loc_vol);
  if(err!=NULL) memset(err,0,sizeof(corr16)*loc_vol);
  if(zm_ave!=NULL) memset(zm_ave,0,sizeof(corr16)*glb_size[0]);
  if(zm_err!=NULL) memset(zm_err,0,sizeof(corr16)*glb_size[0]);

  //summs n estimates
  for(int i=0;i<n;i++)
    {
      MASTER_PRINTF("Elaborating estimate %d of %d\n",i,n);
      
      compute_meson_exchange_correction_stochastically(corr,q_prop,qu,gl);
      memset(zm_corr,0,sizeof(corr16)*glb_size[0]);
      NISSA_LOC_VOL_LOOP(ivol)
        {
	  int t=glb_coord_of_loclx[ivol][0];
	  for(int ig=0;ig<16;ig++)
	    for(int ri=0;ri<2;ri++)
	      {
		zm_corr[t][ig][ri]+=corr[ivol][ig][ri];
		if(ave!=NULL)
		  {
		    ave[ivol][ig][ri]+=corr[ivol][ig][ri];
		    if(err!=NULL) err[ivol][ig][ri]+=corr[ivol][ig][ri]*corr[ivol][ig][ri];
		  }
	      }
	}
      
      //compute zm ave and err
      if(zm_ave!=NULL)
	for(int t=0;t<glb_size[0];t++)
	  for(int ig=0;ig<16;ig++)
	    for(int ri=0;ri<2;ri++)
	      {
		zm_corr[t][ig][ri]=glb_reduce_double(zm_corr[t][ig][ri]);
		zm_ave[t][ig][ri]+=zm_corr[t][ig][ri];
		if(zm_err!=NULL) zm_err[t][ig][ri]+=zm_corr[t][ig][ri]*zm_corr[t][ig][ri];
	      }
    }
  
  //average the full x dist
  if(ave!=NULL)
    NISSA_LOC_VOL_LOOP(ivol)
      for(int ig=0;ig<16;ig++)
	for(int ri=0;ri<2;ri++)
	  {
	    ave[ivol][ig][ri]/=n;
	    if(err!=NULL)
	      {
		err[ivol][ig][ri]/=n;
		err[ivol][ig][ri]=sqrt((err[ivol][ig][ri]-ave[ivol][ig][ri]*ave[ivol][ig][ri])/(n-1));
	      }
	  }
  
  //average the zero momentum
  if(zm_ave!=NULL)
    for(int t=0;t<glb_size[0];t++)
      for(int ig=0;ig<16;ig++)
	for(int ri=0;ri<2;ri++)
	  {
	    zm_ave[t][ig][ri]/=n;
	    if(zm_err!=NULL)
	      {
		zm_err[t][ig][ri]/=n;
		zm_err[t][ig][ri]=sqrt((zm_err[t][ig][ri]-zm_ave[t][ig][ri]*zm_ave[t][ig][ri])/(n-1));
	      }
	  }
  
  nissa_free(corr);
  nissa_free(zm_corr);
  nissa_free(q_prop);
}

//multiple sources
void compute_meson_exchange_correction_stochastically(corr16 *ave,quark_info qu,gluon_info gl,int n)
{
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol+bord_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  
  //reset the summs
  vector_reset(ave);

  //summs n estimates
  for(int i=0;i<n;i++)
    {
      MASTER_PRINTF("Elaborating estimate %d of %d\n",i,n);
      
      compute_meson_exchange_correction_stochastically(corr,q_prop,qu,gl);

      NISSA_LOC_VOL_LOOP(ivol)
	for(int ig=0;ig<16;ig++)
	  for(int ri=0;ri<2;ri++)
	    ave[ivol][ig][ri]+=corr[ivol][ig][ri]/n;
    }
  
  nissa_free(corr);
  nissa_free(q_prop);
}
