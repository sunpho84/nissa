#include "new_types.hpp"

namespace cuda
{
  void float_gauge_field::load_from_lx_conf(nissa::quad_su3 *conf)
  {
    float4 *buf=new float4[12*nissa::loc_vol];
    
    //not really working: mu
    /*
    for(int ivol=0;ivol<nissa::loc_vol;ivol++)
      for(int ifloat=0;ifloat<12;ifloat++)
	buf[ifloat*nissa::loc_vol+ivol]=((double*)(conf[ivol]))[ifloat];
    */
    delete [] buf;
  }
  
}
