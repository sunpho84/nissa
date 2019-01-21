#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3_op.hpp"

#include "dirac_operator_stDLeb_portable.cpp"
#include "dirac_operator_stD_portable.cpp"
#include "dirac_operator_stD_32_portable.cpp"
#include "measures/fermions/stag.hpp"

namespace nissa
{
  //return the even part of the application of D to a vector
  void evn_apply_stD(color *out,quad_su3 **conf,double m,color **in,double sign=1)
  {
    apply_stDeo_half(out,conf,in[ODD]);
    double_vector_linear_comb((double*)out,(double*)in[EVN],m,(double*)out,sign*2,2*NCOL*loc_volh);
  }
  
  //return the odd part of the application of D to a vector
  void odd_apply_stD(color *out,quad_su3 **conf,double m,color **in,double sign=1)
  {
    apply_st2Doe(out,conf,in[EVN]);
    double_vector_linear_comb((double*)out,(double*)in[ODD],m,(double*)out,sign*0.5,2*NCOL*loc_volh);
  }
  
  //return the result of the application of D to a vector
  void apply_stD(color **out,quad_su3 **conf,double m,color **in)
  {
    evn_apply_stD(out[EVN],conf,m,in);
    odd_apply_stD(out[ODD],conf,m,in);
  }
  
  void evn_apply_stD_dag(color *out,quad_su3 **conf,double m,color **in)
  {evn_apply_stD(out,conf,m,in,-1);}
  void odd_apply_stD_dag(color *out,quad_su3 **conf,double m,color **in)
  {odd_apply_stD(out,conf,m,in,-1);}

  void apply_Adams(color **out, quad_su3 **conf, quad_u1 **u1b, double m, double m_twisted, color **temp, color **in){
    
    apply_stD(temp,conf,m,in);
    

    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh){
      assign_complex_prod_i(temp[EVN][ivol][0]);
      assign_complex_prod_i(temp[EVN][ivol][1]);
      assign_complex_prod_i(temp[EVN][ivol][2]);
      assign_complex_prod_i(temp[ODD][ivol][0]);
      assign_complex_prod_i(temp[ODD][ivol][1]);
      assign_complex_prod_i(temp[ODD][ivol][2]);
    }
    set_borders_invalid(temp[EVN]);
    set_borders_invalid(temp[ODD]);

    apply_stag_op(out,conf,u1b,stag::GAMMA_5,stag::IDENTITY,in);

    double_vector_linear_comb((double*)out[EVN],(double*)temp[EVN],1.0,(double*)out[EVN],-m_twisted,2*NCOL*loc_volh);
    double_vector_linear_comb((double*)out[ODD],(double*)temp[ODD],1.0,(double*)out[ODD],-m_twisted,2*NCOL*loc_volh);

  }
}
