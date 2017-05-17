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
}
