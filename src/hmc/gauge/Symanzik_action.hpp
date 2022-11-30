#ifndef _SYMANZIK_ACTION_HPP
#define _SYMANZIK_ACTION_HPP

#include "base/field.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  constexpr double C1_WILSON=0;
  constexpr double C1_TLSYM=-1.0/12;
  constexpr double C1_IWASAKI=-0.331;
  
  //return c0
  inline double get_C0(double C1)
  {return 1-8*C1;}
  
  void Symanzik_action(double *action,eo_ptr<quad_su3> eo_conf,double beta,double C1);
  
  void Symanzik_action(double *action,
		       const LxField<quad_su3>& conf,
		       const double& beta,
		       const double& C1);
  
  template <class T> void tlSym_action(double *action,T lx_conf,double beta)
  {return Symanzik_action(action,lx_conf,beta,C1_TLSYM);}
  template <class T> void Iwasaki_action(double *action,T lx_conf,double beta)
  {return Symanzik_action(action,lx_conf,beta,C1_IWASAKI);}
}

#endif
