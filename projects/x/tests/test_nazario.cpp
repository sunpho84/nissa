#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/propagators/tlSym_gluon_propagator.hpp"
#include "../src/types/types_routines.hpp"

void compute_amp_self_x(spinspin *self_prop,spinspin *q_prop,spin1prop *g_prop)
{
  vector_reset(self_prop);
  for(int mu=0;mu<4;mu++)
    {
      dirac_matr g=base_gamma[map_mu[mu]];

      NISSA_LOC_VOL_LOOP(ivol)
        {
	  spinspin t1,t2;
	  unsafe_spinspin_prod_dirac(t1,q_prop[ivol],&g);
	  unsafe_dirac_prod_spinspin(t2,&g,t1);
	  spinspin_summ_the_complex_prod(self_prop[ivol],t2,g_prop[ivol][mu][mu]);
	}
    }
}

void compute_amp_self_mom(spinspin *self_prop,spinspin *q_prop,spin1prop *g_prop,quark_info qu,gluon_info gl)
{
  pass_spinspin_from_x_to_mom_space(q_prop,q_prop,qu.bc);
  master_printf("Passed quark prop to mom\n");
  pass_spinspin_from_x_to_mom_space(g_prop,g_prop,gl.bc);
  master_printf("Passed gluon prop to mom\n");
  
  vector_reset(self_prop);
  for(int mu=0;mu<4;mu++)
    {
      dirac_matr g=base_gamma[map_mu[mu]];
      
      NISSA_LOC_VOL_LOOP(imom)
        {
	  NISSA_LOC_VOL_LOOP(jmom)
            {
	      coords k;
	      for(int nu=0;nu<4;nu++) k[nu]=(glb_size[nu]+glb_coord_of_loclx[imom][nu]-glb_coord_of_loclx[jmom][nu])%glb_size[nu];
	      int kmom=glblx_of_coord(k);
	      
	      spinspin t1,t2;
	      safe_spinspin_prod_dirac(t1,q_prop[kmom],&g);
	      unsafe_dirac_prod_spinspin(t2,&g,t1);
	      spinspin_summ_the_complex_prod(self_prop[imom],t2,g_prop[jmom][mu][mu]);
	    }
	  master_printf("mu %d/4, imom=%d/%d\n",mu,imom+1,glb_vol);
	}
    }
  
  pass_spinspin_from_mom_to_x_space(self_prop,self_prop,qu.bc);
  pass_spinspin_from_mom_to_x_space(q_prop,q_prop,qu.bc);
  pass_spinspin_from_mom_to_x_space(g_prop,g_prop,gl.bc);
}

void compute_without_trick(double *corr,spinspin *q_prop,spin1prop *g_prop,quark_info qu,gluon_info gl)
{
  //////////////////// self energy propagator in x space////////////////////////
  
  spinspin *self_prop=nissa_malloc("self_prop",loc_vol,spinspin);
  
  //alternatively, compute in mom or x space
  compute_amp_self_x(self_prop,q_prop,g_prop);
  //compute_amp_self_mom(self_prop,q_prop,g_prop,qu,gl);
  
  ////////////////////// external legs in mom space //////////////////////////////
  
  pass_spinspin_from_x_to_mom_space(self_prop,self_prop,qu.bc);
  
  NISSA_LOC_VOL_LOOP(imom)
    {
      spinspin prop;
      mom_space_twisted_propagator_of_imom(prop,qu,imom);

      safe_spinspin_prod_spinspin(self_prop[imom],prop,self_prop[imom]);
      safe_spinspin_prod_spinspin(self_prop[imom],self_prop[imom],prop);
      spinspin_prodassign_double(self_prop[imom],pow(glb_vol,2));
    }

  pass_spinspin_from_mom_to_x_space(self_prop,self_prop,qu.bc);
  
  ////////////////////////////// compute trace with id //////////////////////////
  
  vector_reset(corr);
  NISSA_LOC_VOL_LOOP(ivol)
    {
      dirac_matr g=base_gamma[0];
      complex c={0,0};
      safe_spinspin_prod_dirac(self_prop[ivol],self_prop[ivol],&g);
      trace_spinspin(c,self_prop[ivol]);
      
      int t=glb_coord_of_loclx[ivol][0];
      corr[t]+=c[0];
    }
  
  nissa_free(self_prop);
}

void compute_with_trick(double *corr,spinspin *q_prop,spin1prop *g_prop,quark_info qu,gluon_info gl)
{
  pass_spinspin_from_x_to_mom_space(q_prop,q_prop,qu.bc);
  master_printf("Passed quark prop to mom\n");
  pass_spin1prop_from_x_to_mom_space(g_prop,g_prop,gl.bc);
  master_printf("Passed gluon prop to mom\n");
  
  complex *corr_mom=nissa_malloc("temp",glb_size[0],complex);
  vector_reset(corr_mom);
  
  for(int mu=0;mu<4;mu++)
    {
      dirac_matr g=base_gamma[map_mu[mu]];

      coords p={0,0,0,0};
      for(p[0]=0;p[0]<glb_size[0];p[0]++)
	{
	  int imom=glblx_of_coord(p);
	  NISSA_LOC_VOL_LOOP(jmom)
	    {
	      coords k;
	      for(int nu=0;nu<4;nu++) k[nu]=(glb_size[nu]+glb_coord_of_loclx[imom][nu]-glb_coord_of_loclx[jmom][nu])%glb_size[nu];
	      int kmom=glblx_of_coord(k);

	      spinspin t1,t2;
	      safe_spinspin_prod_dirac(t1,q_prop[kmom],&g);
	      unsafe_dirac_prod_spinspin(t2,&g,t1);
	      
	      spinspin self_prop;
	      spinspin_summ_the_complex_prod(self_prop,t2,g_prop[jmom][mu][mu]);
	      
	      safe_spinspin_prod_spinspin(self_prop,q_prop[imom],self_prop);
	      safe_spinspin_prod_spinspin(self_prop,self_prop,q_prop[imom]);
	      
	      complex t3;
	      trace_spinspin(t3,self_prop);
	      
	      complex_summassign(corr_mom[p[0]],t3);
	    }
	  master_printf("mu=%d/4, t=%d/%d\n",mu+1,p[0]+1,glb_size[0]);
	}
    }
  
  for(int t=0;t<glb_size[0];t++)
    master_printf("%lg %lg\n",corr_mom[t][0],corr_mom[t][1]);

  
  //temporal fourier transform
  const int use_nissa_fft=0;
  if(use_nissa_fft)
    {
      //fft acts on complex, so buffering
      complex buf[glb_size[0]];
      
      //parameters for fft
      int mu=0;
      int ncompl_per_site=1;
      int sign=-1;
      int normalize=0;
      fft1d(buf,corr_mom,ncompl_per_site,mu,sign,normalize);
      
      //normalize
      for(int t=0;t<glb_size[0];t++)
	corr[t]=buf[t][0]*pow(glb_vol,3)/glb_size[0];
    }
  else
    {
      vector_reset(corr);
      for(int t=0;t<glb_size[0];t++)
	{
	  for(int ip=0;ip<glb_size[0];ip++)
	    {
	      double p=M_PI*(2*ip+qu.bc[0])/glb_size[0];
	      complex ph={cos(p*t),sin(-p*t)};
	      
	      complex t1;
	      unsafe_complex_prod(t1,corr_mom[ip],ph);
	      
	      corr[t]+=t1[0];
	    }
	  
	  corr[t]*=pow(glb_vol,3)/glb_size[0];
	}
    }
  
  nissa_free(corr_mom);
}

int main(int narg,char **arg)
{
  //basic initialization
  int T=48,L=24;
  init_nissa(narg,arg);
  init_grid(T,L);
  
  //quark
  double kappa=0.115;
  double mass=0.00;
  double quark_theta[4]={0,0,0,0};
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  master_printf("Computed quark prop\n");
  
  //gluon
  double alpha=1;
  double gluon_theta[4]={0,0,0,0};
  gluon_info gl=create_Wilson_gluon_info(alpha,gluon_theta);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);
  master_printf("Computed gluon prop\n");
  
  double *corr=nissa_malloc("corr",glb_size[0],double);  
  
  compute_without_trick(corr,q_prop,g_prop,qu,gl);
  //compute_with_trick(corr,q_prop,g_prop,qu,gl);
  
  /////////////////////////////// print correlation /////////////////////////////

  FILE *fout=open_file("/tmp/self.dat","w");
  for(int t=0;t<glb_size[0];t++)
    master_fprintf(fout,"%d %16.16lg\n",t,glb_reduce_double(corr[t]));  
  if(rank==0) fclose(fout);
  
  /////////////////////////////////// close ////////////////////////////////////
  
  nissa_free(corr);
  nissa_free(q_prop);
  nissa_free(g_prop);

  close_nissa();
  
  return 0;
}
