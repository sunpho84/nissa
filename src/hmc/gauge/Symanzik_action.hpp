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
  {
    return 1-8*C1;
  }
  
  double Symanzik_action(const EoField<quad_su3>& eo_conf,
			 const double& beta,
			 const double& C1);
  
  double Symanzik_action(const LxField<quad_su3>& conf,
		       const double& beta,
		       const double& C1);
  
  template <typename T>
  double tlSym_action(T&& lx_conf,
		      const double& beta)
  {
    return Symanzik_action(lx_conf,beta,C1_TLSYM);
  }
  
  template <typename T>
  double Iwasaki_action(T&& lx_conf,
			const double& beta)
  {
    return Symanzik_action(lx_conf,beta,C1_IWASAKI);
  }
}

#endif
