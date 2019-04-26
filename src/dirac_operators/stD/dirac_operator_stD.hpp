#ifndef _DIRAC_OPERATOR_STD_HPP
#define _DIRAC_OPERATOR_STD_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_st2Doe(color *out,quad_su3 **conf,color *in);
  void apply_stD2ee_m2(color *out,quad_su3 **conf,color *temp,double mass2,color *in);
  void apply_stD2Leb_ee_m2(color *out,oct_su3 **conf,color *temp,double mass2,color *in);
  void apply_stD2ee_m2_32(single_color *out,single_quad_su3 **conf,single_color *temp,float mass2,single_color *in);
  void apply_stDeo_half(color *out,quad_su3 **conf,color *in);
  void apply_stDoe(color *out,quad_su3 **conf,color *in);
  void evn_apply_stD(color *out,quad_su3 **conf,double m,color **in,double sign=1);
  void odd_apply_stD(color *out,quad_su3 **conf,double m,color **in,double sign=1);
  void evn_apply_stD_dag(color *out,quad_su3 **conf,double m,color **in);
  void odd_apply_stD_dag(color *out,quad_su3 **conf,double m,color **in);
  void apply_stD(color **out,quad_su3 **conf,double m,color **in);
  
  void apply_Adams(color **out,quad_su3 **conf,quad_u1 **u1b,double m,double m_Adams,color **temp,color **in);
  void apply_AdamsII(color **out,quad_su3 **conf,quad_u1 **u1b,double m,double m_Adams,color **temp,color **in);
}

#endif
