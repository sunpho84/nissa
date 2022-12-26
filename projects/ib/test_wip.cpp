#include <nissa.hpp>

// #include "base/quda_bridge.hpp"

using namespace nissa;

// double rel_diff_norm(spincolor *test,spincolor *ref)
// {
//   double_vector_subtassign((double*)test,(double*)ref,locVol*sizeof(spincolor)/sizeof(double));
//   const double norm2_diff=double_vector_glb_norm2(test,locVol);
//   const double norm2_ref=double_vector_glb_norm2(ref,locVol);
//   const double res=sqrt(norm2_diff/norm2_ref);
  
//   return res;
// }

CUDA_MANAGED double* pt;
CUDA_MANAGED int io;

namespace nissa
{
  template <typename IMin,
	    typename IMax,
	    typename F>
  __global__
  void cuda_generic_kernel2(const IMin min,
			    const IMax max,
			    const F& f)
  {
    const auto i=min+blockIdx.x*blockDim.x+threadIdx.x;
    if(i<max)
      f(i);
  }
  
  template <typename IMin,
	    typename IMax,
	    typename F>
  void cuda_parallel_for2(const int line,
			 const char* file,
			 const IMin min,
			 const IMax max,
			 const F& f)
  {
    const auto length=(max-min);
    const dim3 block_dimension(NUM_THREADS);
    const dim3 grid_dimension((length+block_dimension.x-1)/block_dimension.x);
    
    double initTime=0;
    extern int rank,verbosity_lv;
    const bool print=(verbosity_lv>=1// 2
		      and rank==0);
    if(print)
      {
	printf("at line %d of file %s launching kernel on loop [%ld,%ld) using blocks of size %d and grid of size %d\n",
	   line,file,(int64_t)min,(int64_t)max,block_dimension.x,grid_dimension.x);
	initTime=take_time();
      }
    
    if(length>0)
      {
	cuda_generic_kernel2<<<grid_dimension,block_dimension>>>(min,max,f);
	cudaDeviceSynchronize();
      }
    
    if(print)
      printf(" finished in %lg s\n",take_time()-initTime);
  }
}

void in_main(int narg,char **arg)
{
  const int T=16,L=16;
  
  init_grid(T,L);
  
  LxField<quad_su3> conf("conf");
  
  master_printf("allocated in %p\n",conf._data);
  {
    auto c=conf.getWritable();
    master_printf("allocated in %p\n",c._data);
    double& e=c[locVol-1][3][2][2][1];
    master_printf("end: %p, should be %p\n",&e,c._data+locVol*4*3*3*2);
    
  }

  cuda_parallel_for2(
      44, "/home/francesco/QCD/SORGENTI/nissa_origi/projects/ib/test_wip.cpp",
      0, 1, [confa = conf.getReadable()]__device__(const int &ivol) {
        pt = confa._data;
        io = confa.externalSize;
      });
  master_printf("data: %p external_size: %d\n",pt,io);
  // start_loc_rnd_gen(235235);
  
  // spincolor *in=nissa_malloc("in",locVol+bord_vol,spincolor);
  // spincolor *out=nissa_malloc("out",locVol+bord_vol,spincolor);
  // spincolor *tmp=nissa_malloc("tmp",locVol+bord_vol,spincolor);
  
  // /// First test: load a spincolor and unload it
  // generate_undiluted_source(in,RND_Z2,-1);
  
  // quda_iface::remap_nissa_to_quda(tmp,in);
  // quda_iface::remap_quda_to_nissa(out,tmp);
  
  // master_printf("testing map and unmap, residue: %lg\n",rel_diff_norm(out,in));
  
  // /// Second test: apply the dirac operator
  
  // quad_su3 *conf=nissa_malloc("conf",locVol+bord_vol+edge_vol,quad_su3);
  // spincolor *out_nissa=nissa_malloc("out_nissa",locVol+bord_vol,spincolor);
  
  // generate_hot_lx_conf(conf);
  // master_printf("plaq: %lg\n",global_plaquette_lx_conf(conf));
  
  // vector_reset(in);
  // in[0][0][0][0]=1.0;
  
  // const double kappa=0.1325,csw=1.345,mu=0.243;
  // quda_iface::apply_tmD(out,conf,kappa,csw,mu,in);
  
  // clover_term_t *Cl=nissa_malloc("Cl",locVol,clover_term_t);
  // inv_clover_term_t *invCl=nissa_malloc("invCl",locVol,inv_clover_term_t);
  // clover_term(Cl,csw,conf);
  // invert_twisted_clover_term(invCl,mu*tau3[0],kappa,Cl);
  // apply_tmclovQ(out_nissa,conf,kappa,Cl,mu,in);
  // nissa_free(invCl);
  // nissa_free(Cl);
  
  // safe_dirac_prod_spincolor(out_nissa,base_gamma[5],out_nissa);
  
  // master_printf("comparing\n");
  // for(int ivol=0;ivol<locVol;ivol++)
  //   for(int id=0;id<NDIRAC;id++)
  //     for(int ic=0;ic<NCOL;ic++)
  // 	for(int ri=0;ri<2;ri++)
  // 	  {
  // 	    const double n=out_nissa[ivol][id][ic][ri];
  // 	    const double q=out[ivol][id][ic][ri];
  // 	    if(fabs(n)>1e-10 or fabs(q)>1e-10)
  // 	      master_printf("out,[nissa,quda][ivol=%d,id=%d,ic=%d,ri=%d]: %lg %lg\n",
  // 			    ivol,id,ic,ri,n,q);
  // 	  }
  // master_printf("testing tmD, residue: %lg\n",rel_diff_norm(out,out_nissa));
  
  // nissa_free(tmp);
  // nissa_free(in);
  // nissa_free(out);
  // nissa_free(out_nissa);
  // nissa_free(conf);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
