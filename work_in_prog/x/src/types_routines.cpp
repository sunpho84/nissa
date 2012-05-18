#include <string.h>

#include "../src/types.h"

gluon_info create_tlSym_gluon_info(double alpha,momentum_t bc,double c1=-1.0/12,double zmp=0)
{
  gluon_info out;
  memcpy(out.bc,bc,sizeof(momentum_t));
  out.alpha=alpha;
  out.zmp=zmp;
  
  return out;
}

gluon_info create_Wilson_gluon_info(double alpha,momentum_t bc,double zmp=0)
{return create_tlSym_gluon_info(alpha,bc,0,zmp);}

quark_info create_twisted_quark_info(double kappa,double mass,momentum_t bc)
{
  quark_info out;
  memcpy(out.bc,bc,sizeof(momentum_t));
  out.kappa=kappa;
  out.mass=mass;
  
  return out;
}

quark_info create_Wilson_quark_info(double kappa,momentum_t bc)
{return create_twisted_quark_info(kappa,0,bc);}
