#include <nissa.hpp>

using namespace nissa;

void in_main(int narg,char **arg)
{
  const int T=16,L=16;
  
  init_grid(T,L);
  
  start_loc_rnd_gen(235235);
  
  spincolor *in=nissa_malloc("in",loc_vol,spincolor);
  spincolor *out=nissa_malloc("out",loc_vol,spincolor);
  double *tmp=nissa_malloc("tmp",loc_vol*sizeof(spincolor)/sizeof(double),double);
  
  /// First test: load a spincolor and unload it
  generate_undiluted_source(in,RND_Z2,-1);
  
  quda_iface::remap_nissa_to_quda(tmp,in);
  quda_iface::remap_quda_to_nissa(out,tmp);
  
  double_vector_subtassign((double*)out,(double*)in,loc_vol*sizeof(spincolor)/sizeof(double));
  const double norm2_diff=double_vector_glb_norm2(out,loc_vol);
  const double norm2_in=double_vector_glb_norm2(in,loc_vol);
  const double res=sqrt(norm2_diff/norm2_in);
  master_printf("testing map and unmap, residue: %lg\n",res);
  
  /// Second test: apply the dirac operator
  
  quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);
  
  generate_cold_lx_conf(conf);
  
  quda_iface::apply_tmD(out,conf,0.125,0.0,in);
  
  nissa_free(tmp);
  nissa_free(in);
  nissa_free(out);
  nissa_free(conf);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
