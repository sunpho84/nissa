#ifndef _SYMANZIK_ACTION_HPP
#define _SYMANZIK_ACTION_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  const double C1_WILSON=0;
  const double C1_TLSYM=-1.0/12;
  const double C1_IWASAKI=-0.331;
  
  //return c0
  inline double get_C0(double C1)
  {return 1-8*C1;}
  
  void Symanzik_action(double *action,quad_su3 **eo_conf,double beta,double C1);
  void Symanzik_action(double *action,quad_su3 *lx_conf,double beta,double C1);
  template <class T> void tlSym_action(double *action,T lx_conf,double beta)
  {return Symanzik_action(action,lx_conf,beta,C1_TLSYM);}
  template <class T> void Iwasaki_action(double *action,T lx_conf,double beta)
  {return Symanzik_action(action,lx_conf,beta,C1_IWASAKI);}
}

#endif
