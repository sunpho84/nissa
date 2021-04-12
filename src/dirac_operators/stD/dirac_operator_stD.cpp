#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <math.h>

#include <base/bench.hpp>
#include <base/debug.hpp>
#include <base/vectors.hpp>
#include <communicate/borders.hpp>
#include <linalgs/linalgs.hpp>
#include <new_types/su3_op.hpp>

#include <dirac_operators/stD/dirac_operator_stD_portable.cpp>
#include <dirac_operators/stD/dirac_operator_stD_32_portable.cpp>
#include <measures/fermions/stag.hpp>

namespace nissa
{
  //return the even part of the application of D to a vector
  void evn_apply_stD(color *out,eo_ptr<quad_su3> conf,double m,eo_ptr<color> in,double sign=1)
  {
    apply_stDeo_half(out,conf,in[ODD]);
    double_vector_linear_comb((double*)out,(double*)in[EVN],m,(double*)out,sign*2,2*NCOL*locVolh.nastyConvert());
  }
  
  //return the odd part of the application of D to a vector
  void odd_apply_stD(color *out,eo_ptr<quad_su3> conf,double m,eo_ptr<color> in,double sign=1)
  {
    apply_st2Doe(out,conf,in[EVN]);
    double_vector_linear_comb((double*)out,(double*)in[ODD],m,(double*)out,sign*0.5,2*NCOL*locVolh());
  }
  
  //return the result of the application of D to a vector
  void apply_stD(eo_ptr<color> out,eo_ptr<quad_su3> conf,double m,eo_ptr<color> in)
  {
    evn_apply_stD(out[EVN],conf,m,in);
    odd_apply_stD(out[ODD],conf,m,in);
  }
  
  void evn_apply_stD_dag(color *out,eo_ptr<quad_su3> conf,double m,eo_ptr<color> in)
  {
    evn_apply_stD(out,conf,m,in,-1);
  }
  
  void odd_apply_stD_dag(color *out,eo_ptr<quad_su3> conf,double m,eo_ptr<color> in)
  {
    odd_apply_stD(out,conf,m,in,-1);
  }
  
  //Adams operator https://arxiv.org/pdf/1103.6191.pdf
  void apply_Adams(eo_ptr<color> out,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double m,double m_Adams,eo_ptr<color> temp,eo_ptr<color> in)
  {
    // out = g5 X id * in
    apply_stag_op(out,conf,u1b,stag::GAMMA_5,stag::IDENTITY,in);
    
    // temp = D * in
    add_backfield_with_stagphases_to_conf(conf,u1b);
    apply_stD(temp,conf,m,in);
    rem_backfield_with_stagphases_from_conf(conf,u1b);
    
    for(int eo=0;eo<2;eo++)
      {
	
	// temp = i * D * in
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  for(int ic=0;ic<NCOL;ic++)
	    assign_complex_prod_i(temp[eo][ieo.nastyConvert()][ic]);
	NISSA_PARALLEL_LOOP_END;
	
	// out = (i * D - m * g5 X id) * in
	double_vector_summ_double_vector_prod_double((double*)out[eo],(double*)temp[eo],(double*)out[eo],-m_Adams,2*NCOL*locVolh.nastyConvert());
      }
  }
  
  //AdamsII operator https://arxiv.org/pdf/1008.2833.pdf
  void apply_AdamsII(eo_ptr<color> out,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double m,double m_Adams,eo_ptr<color> temp,eo_ptr<color> in)
  {
    // out = D * in
    add_backfield_with_stagphases_to_conf(conf,u1b);
    apply_stD(out,conf,m,in);
    rem_backfield_with_stagphases_from_conf(conf,u1b);
    
    // temp = (Gamma5 x Gamma5) * D * in
    apply_stag_op(temp,conf,u1b,stag::GAMMA_5,stag::GAMMA_5,out);
    
    // out = g5 X id * in
    apply_stag_op(out,conf,u1b,stag::GAMMA_5,stag::IDENTITY,in);
    
    // out = (g5 X g5 * D - m * g5 X id) * in
    for(int eo=0;eo<2;eo++)
      double_vector_summ_double_vector_prod_double((double*)out[eo],(double*)temp[eo],(double*)out[eo],-m_Adams,2*NCOL*locVolh.nastyConvert());
  }
}
