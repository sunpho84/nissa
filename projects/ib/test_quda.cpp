#include <nissa.hpp>

#include "base/quda_bridge.hpp"

using namespace nissa;

double rel_diff_norm(spincolor *test,spincolor *ref)
{
  double_vector_subtassign((double*)test,(double*)ref,locVol*sizeof(spincolor)/sizeof(double));
  const double norm2_diff=double_vector_glb_norm2(test,locVol);
  const double norm2_ref=double_vector_glb_norm2(ref,locVol);
  const double res=sqrt(norm2_diff/norm2_ref);
  
  return res;
}

void in_main(int narg,char **arg)
{
  const int T=16,L=16;
  
  init_grid(T,L);
  
  start_loc_rnd_gen(235235);
  
  spincolor *in=nissa_malloc("in",locVol,spincolor);
  spincolor *out=nissa_malloc("out",locVol,spincolor);
  spincolor *tmp=nissa_malloc("tmp",locVol,spincolor);
  
  /// First test: load a spincolor and unload it
  generate_undiluted_source(in,RND_Z2,-1);
  
  quda_iface::remap_nissa_to_quda(tmp,in);
  quda_iface::remap_quda_to_nissa(out,tmp);
  
  master_printf("testing map and unmap, residue: %lg\n",rel_diff_norm(out,in));
  
  /// Second test: apply the dirac operator
  
  quad_su3 *conf=nissa_malloc("conf",locVol+bord_vol,quad_su3);
  spincolor *out_nissa=nissa_malloc("out_nissa",locVol,spincolor);
  
  generate_hot_lx_conf(conf);
  master_printf("plaq: %lg\n",global_plaquette_lx_conf(conf));
  
  vector_reset(in);
  in[0][0][0][0]=1.0;
  
  const double kappa=0.125,mu=0.0;
  quda_iface::apply_tmD(out,conf,kappa,mu,in);
  apply_tmQ(out_nissa,conf,kappa,mu,in);
  
  safe_dirac_prod_spincolor(out_nissa,base_gamma+5,out_nissa);
  
  master_printf("comparing\n");
  for(int ivol=0;ivol<locVol;ivol++)
    for(int id=0;id<NDIRAC;id++)
      for(int ic=0;ic<NCOL;ic++)
	for(int ri=0;ri<2;ri++)
	  {
	    const double n=out_nissa[ivol][id][ic][ri];
	    const double q=out[ivol][id][ic][ri];
	    if(fabs(n)>1e-10 or fabs(q)>1e-10)
	      master_printf("out,[nissa,quda][ivol=%d,id=%d,ic=%d,ri=%d]: %lg %lg\n",
			    ivol,id,ic,ri,n,q);
	  }
  master_printf("testing tmD, residue: %lg\n",rel_diff_norm(out,out_nissa));
  
  nissa_free(tmp);
  nissa_free(in);
  nissa_free(out);
  nissa_free(out_nissa);
  nissa_free(conf);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
