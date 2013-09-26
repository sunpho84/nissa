#include <string.h>

#include "../../../../src/nissa.h"

#include "../types/types.h"

gluon_info create_tlSym_gluon_info(double alpha,momentum_t bc,double c1=-1.0/12,double zmp=0)
{
  gluon_info out;
  memcpy(out.bc,bc,sizeof(momentum_t));
  out.alpha=alpha;
  out.c1=c1;
  out.zmp=zmp;
  
  return out;
}

gluon_info create_Wilson_gluon_info(double alpha,momentum_t bc,double zmp=0)
{return create_tlSym_gluon_info(alpha,bc,0,zmp);}

quark_info create_twisted_quark_info(double kappa,double mass,momentum_t bc,double zmp=0)
{
  quark_info out;
  memcpy(out.bc,bc,sizeof(momentum_t));
  out.kappa=kappa;
  out.zmp=zmp;
  out.mass=mass;
  
  return out;
}

quark_info create_Wilson_quark_info(double kappa,momentum_t bc)
{return create_twisted_quark_info(kappa,0,bc);}

void get_spin_from_spinspin(spin *out,spinspin *in,int id_so)
{
  int all=check_borders_allocated((void*)in);
  
  int ending;
  if(all)
    {
      communicate_lx_spinspin_borders(in);
      ending=loc_vol+bord_vol;
    }
  else
    ending=loc_vol;
  
  for(int ivol=0;ivol<ending;ivol++)
    for(int id_si=0;id_si<4;id_si++)
      memcpy(out[ivol][id_si],in[ivol][id_si][id_so],sizeof(complex));
  
  if(all) set_borders_valid(out);
  else    set_borders_invalid(out);
}

void put_spin_into_spinspin(spinspin *out,spin *in,int id_so)
{
  NISSA_LOC_VOL_LOOP(ivol)
    for(int id_si=0;id_si<4;id_si++)
      memcpy(out[ivol][id_si][id_so],in[ivol][id_si],sizeof(complex));
  set_borders_invalid(out);
}
