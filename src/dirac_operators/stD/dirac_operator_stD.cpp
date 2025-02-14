#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/debug.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  //Adams operator https://arxiv.org/pdf/1103.6191.pdf
  void apply_Adams(eo_ptr<color> out,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double m,double m_Adams,eo_ptr<color> temp,eo_ptr<color> in)
  {
    CRASH("reimplement");
    
    // // out = g5 X id * in
    // apply_stag_op(out,conf,u1b,stag::GAMMA_5,stag::IDENTITY,in);
    
    // // temp = D * in
    // add_backfield_with_stagphases_to_conf(conf,u1b);
    // apply_stD(temp,conf,m,in);
    // rem_backfield_with_stagphases_from_conf(conf,u1b);
    
    // for(int eo=0;eo<2;eo++)
    //   {
	
    // 	// temp = i * D * in
    // 	NISSA_PARALLEL_LOOP(ivol,0,locVolh)
    // 	  for(int ic=0;ic<NCOL;ic++)
    // 	    assign_complex_prod_i(temp[eo][ivol][ic]);
    // 	NISSA_PARALLEL_LOOP_END;
	
    // 	// out = (i * D - m * g5 X id) * in
    // 	double_vector_summ_double_vector_prod_double((double*)out[eo],(double*)temp[eo],(double*)out[eo],-m_Adams,2*NCOL*locVolh);
    //   }
  }
  
  //AdamsII operator https://arxiv.org/pdf/1008.2833.pdf
  void apply_AdamsII(eo_ptr<color> out,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double m,double m_Adams,eo_ptr<color> temp,eo_ptr<color> in)
  {
    CRASH("reimplement");
    
    // // out = D * in
    // add_backfield_with_stagphases_to_conf(conf,u1b);
    // apply_stD(out,conf,m,in);
    // rem_backfield_with_stagphases_from_conf(conf,u1b);
    
    // // temp = (Gamma5 x Gamma5) * D * in
    // apply_stag_op(temp,conf,u1b,stag::GAMMA_5,stag::GAMMA_5,out);
    
    // // out = g5 X id * in
    // apply_stag_op(out,conf,u1b,stag::GAMMA_5,stag::IDENTITY,in);
    
    // // out = (g5 X g5 * D - m * g5 X id) * in
    // for(int eo=0;eo<2;eo++)
    //   double_vector_summ_double_vector_prod_double((double*)out[eo],(double*)temp[eo],(double*)out[eo],-m_Adams,2*NCOL*locVolh);
  }
}
