#ifndef _SYMANZIK_ACTION_HPP
#define _SYMANZIK_ACTION_HPP

namespace nissa
{
  const double C1_WILSON=0;
  const double C1_TLSYM=-1.0/12;
  const double C1_IWASAKI=-0.331;
  
  //return c0
  inline double get_C0(double C1,bool stagphases_present)
  {return (1-8*C1)*(1-2*stagphases_present);}
  
  void Symanzik_action(double *action,quad_su3 **eo_conf,double beta,double C1,bool stag_phases_present);
  void Symanzik_action(double *action,quad_su3 *lx_conf,double beta,double C1,bool stag_phases_present);
  template <class T> void tlSym_action(double *action,T lx_conf,double beta,bool stag_phases_present)
  {return Symanzik_action(action,lx_conf,beta,C1_TLSYM,stag_phases_present);}
  template <class T> void Iwasaki_action(double *action,T lx_conf,double beta,bool stag_phases_present)
  {return Symanzik_action(action,lx_conf,beta,C1_IWASAKI,stag_phases_present);}
}

#endif
