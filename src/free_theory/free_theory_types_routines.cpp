#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "free_theory_types.hpp"

namespace nissa
{
  gauge_info create_tlSym_gauge_info(double alpha,momentum_t bc,double c1=-1.0/12)
  {
    gauge_info out;
    memcpy(out.bc,bc,sizeof(momentum_t));
    out.alpha=alpha;
    out.c1=c1;
    
    return out;
  }
  
  gauge_info create_Wilson_gauge_info(double alpha,momentum_t bc)
  {return create_tlSym_gauge_info(alpha,bc,0);}
  
  tm_quark_info create_twisted_quark_info(double kappa,double mass,momentum_t bc,int r,double zmp=0)
  {
    tm_quark_info out;
    memcpy(out.bc,bc,sizeof(momentum_t));
    out.kappa=kappa;
    out.zmp=zmp;
    out.mass=mass;
    out.r=r;
    
    return out;
  }
  
  tm_quark_info create_Wilson_quark_info(double kappa,momentum_t bc)
  {return create_twisted_quark_info(kappa,0,bc,0);}
  
  void get_spin_from_spinspin(spin *out,spinspin *in,int id_so)
  {
    int all=check_borders_allocated((void*)in,0);
    
    LocLxSite ending;
    if(all)
      {
	communicate_lx_spinspin_borders(in);
	ending=locVol+bord_vol;
      }
    else
      ending=locVol;
    
    for(LocLxSite ivol=0;ivol<ending;ivol++)
      for(int id_si=0;id_si<4;id_si++)
	memcpy(out[ivol.nastyConvert()][id_si],in[ivol.nastyConvert()][id_si][id_so],sizeof(complex));
    
    if(all) set_borders_valid(out);
    else    set_borders_invalid(out);
  }
  
  void put_spin_into_spinspin(spinspin *out,spin *in,int id_so)
  {
    NISSA_LOC_VOL_LOOP(ivol)
      for(int id_si=0;id_si<4;id_si++)
	memcpy(out[ivol.nastyConvert()][id_si][id_so],in[ivol.nastyConvert()][id_si],sizeof(complex));
    set_borders_invalid(out);
  }
}
