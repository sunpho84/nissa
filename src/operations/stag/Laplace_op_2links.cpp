#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3_op.hpp"

#include "stag.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  using namespace stag;
 
  //covariant shift backward: i=i+mu
  void cshift_bw(color *out,quad_su3 *conf,int mu,color *in,bool reset_first=true)
  {
    GET_THREAD_ID();
    
    communicate_lx_color_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	if(reset_first) color_put_to_zero(out[ivol]);
	su3_summ_the_prod_color(out[ivol],conf[ivol][mu],in[loclx_neighup[ivol][mu]]);
      }
    set_borders_invalid(out);
  }
  
  //covariant shift forward: i=i-mu
  void cshift_fw(color *out,quad_su3 *conf,int mu,color *in,bool reset_first=true)
  {
    GET_THREAD_ID();
    
    communicate_lx_color_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	if(reset_first) color_put_to_zero(out[ivol]);
	su3_dag_summ_the_prod_color(out[ivol],conf[loclx_neighdw[ivol][mu]][mu],in[loclx_neighdw[ivol][mu]]);
      }
    set_borders_invalid(out);
  }
  
  //multiply by the 2-links Laplace operator
  THREADABLE_FUNCTION_3ARG(Laplace_operator_2_links, color*,out, quad_su3*,conf, color*,in)
  {
    color *temp=nissa_malloc("temp",loc_vol+bord_vol,color);
    int nentries=loc_vol*sizeof(color)/sizeof(double);
    
    vector_reset(out);
    
    for(int mu=0;mu<NDIM;mu++)
      {
	cshift_bw(temp,conf,mu,in);
	cshift_bw(out,conf,mu,temp,false);
	cshift_fw(temp,conf,mu,in);
	cshift_fw(out,conf,mu,temp,false);
      }
    double_vector_prodassign_double((double*)out,0.25,nentries);
    double_vector_summassign_double_vector_prod_double((double*)out,(double*)in,-NDIM/2.0,nentries);
    
    nissa_free(temp);
  }
  THREADABLE_FUNCTION_END
}
