#include "nissa.hpp"

using namespace nissa;

//compute the matrix element of the conserved current between two propagators
THREADABLE_FUNCTION_4ARG(calc_cur, spin1field*,cur, spincolor*,source, quad_su3*,conf, spincolor*,prop)
{
  GET_THREAD_ID();
  
  vector_reset(cur);
  
  communicate_lx_spincolor_borders(source);
  communicate_lx_spincolor_borders(prop);
  communicate_lx_quad_su3_borders(conf);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<NDIM;mu++)
      {
	// int ifw=loclx_neighup[ivol][mu];
	// spincolor f;
	
	// //piece psi_ivol U_ivol psi_fw
	// unsafe_su3_prod_spincolor(f,conf[ivol][mu],prop[ifw]);
	// spincolor_scalar_prod(cur[ivol][mu],source[ivol],f);
	
	int ifw=loclx_neighup[ivol][mu];
	spincolor f,Gf;
	complex c;
	
	//piece psi_ivol U_ivol psi_fw
	unsafe_su3_prod_spincolor(f,conf[ivol][mu],prop[ifw]);
	spincolor_copy(Gf,f);
	dirac_subt_the_prod_spincolor(Gf,base_gamma+igamma_of_mu[mu],f);
	spincolor_scalar_prod(c,source[ivol],Gf);
	complex_summ_the_prod_idouble(cur[ivol][mu],c,-0.5);
	
	//piece psi_fw U_ivol^dag psi_ivol
	unsafe_su3_dag_prod_spincolor(f,conf[ivol][mu],prop[ivol]);
	spincolor_copy(Gf,f);
	dirac_summ_the_prod_spincolor(Gf,base_gamma+igamma_of_mu[mu],f);
	spincolor_scalar_prod(c,source[ifw],Gf);
	complex_summ_the_prod_idouble(cur[ivol][mu],c,+0.5);
      }
  set_borders_invalid(cur);
}
THREADABLE_FUNCTION_END

void in_main(int narg,char **arg)
{
  init_grid(16,8);
  
  spincolor *source=nissa_malloc("source",loc_vol,spincolor);
  spincolor *prop=nissa_malloc("prop",loc_vol,spincolor);
  quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);
  
  read_ildg_gauge_conf(conf,"conf");
  momentum_t th={1,0,0,0};
  put_boundaries_conditions(conf,th,0,0);
  
  read_real_vector(source,"source.0000.00000","scidac-binary-data");
  read_real_vector(prop,"source.0000.00000.inverted","scidac-binary-data");
  
  spin1field *cur=nissa_malloc("cur",loc_vol,spin1field);
  calc_cur(cur,source,conf,prop);
  
  FILE *fout=open_file("cur.txt","w");
  master_fprintf(fout,"# [loops_em]  norm square of source     = %.16lg\n",double_vector_glb_norm2(source,loc_vol));
  master_fprintf(fout,"# [loops_em]  norm square of propagator = %.16lg\n",double_vector_glb_norm2(prop,loc_vol));
  NISSA_LOC_VOL_LOOP(ivol)
  {
    master_fprintf(fout,"# [loops_em] x\t");
    for(int mu=0;mu<NDIM;mu++) master_fprintf(fout,"%d ",loc_coord_of_loclx[ivol][mu]);
    master_fprintf(fout,"\n");
    
    for(int mu=0;mu<NDIM;mu++) master_fprintf(fout,"%d" "\t" "%.16lg" "\t" "%.16lg\n",mu,cur[ivol][mu][RE],cur[ivol][mu][IM]);
  }
  close_file(fout);
  
  nissa_free(cur);
  nissa_free(conf);
  nissa_free(prop);
  nissa_free(source);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
