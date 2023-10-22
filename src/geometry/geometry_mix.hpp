#ifndef _GEOMETRY_MIX_HPP
#define _GEOMETRY_MIX_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <unistd.h>

#include <base/bench.hpp>
#include <base/field.hpp>
#include <geometry/geometry_eo.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  void paste_eo_parts_into_lx_vector_internal(void *out_lx,eo_ptr<void> in_eo,size_t bps);
  
  template <typename LX,
	    typename EO>
  void paste_eo_parts_into_lx_vector(FieldFeat<LX>& _outLx,
				     const EO& inEo)
  {
    START_TIMING(remap_time,nremap);
    
    //paste
    FOR_BOTH_PARITIES(par,
    {
      PAR(0,locVolh,
	  CAPTURE(par,outLx=(~_outLx).getWritable(),
		  TO_READ(inEo)),
	  eo,
	  {
	    for(int internalDeg=0;internalDeg<inEo[par].nInternalDegs;internalDeg++)
	      outLx(loclx_of_loceo[par][eo],internalDeg)=inEo[par](eo,internalDeg);
	  });
    });
    
    STOP_TIMING(remap_time);
  }
  
  /////////////////////////////////////////////////////////////////
  
  template <typename EO,
	    typename LX>
  void split_lx_vector_into_eo_parts(EO&& outEo,
				     const LX& inLx)
  {
    START_TIMING(remap_time,nremap);
    
    PAR(0,locVol,
	CAPTURE(TO_WRITE(outEo),
		TO_READ(inLx)),locLx,
	{
	  for(int internalDeg=0;internalDeg<LX::nInternalDegs;internalDeg++)
	    outEo[loclx_parity[locLx]](loceo_of_loclx[locLx],internalDeg)=inLx(locLx,internalDeg);
	});
    
    STOP_TIMING(remap_time);
  }
  
  /////////////////////////////////////////////////////////////////
  
  template <typename EoO,
	    typename LX>
  void get_evn_or_odd_part_of_lx_vector(EoO&& outEo,
					const LX& inLx,
					const FieldCoverage& par)
  {
    START_TIMING(remap_time,nremap);
    
    //get
    PAR(0,locVolh,
	CAPTURE(par,
		TO_WRITE(outEo),
		TO_READ(inLx)),locEo,
	{
	  for(int internalDeg=0;internalDeg<inLx.nInternalDegs;internalDeg++)
	    outEo(locEo,internalDeg)=inLx(loclx_of_loceo[par][locEo],internalDeg);
	});
    STOP_TIMING(remap_time);
  }
  
  /////////////////////////////////////////////////////////////////
  
  void remap_vector_internal(char *out,char *in,size_t bps,int *dest_of_source,int length);
  
  // template <class T> void remap_Leb_ev_or_od_to_loc_vector(T *out,T *in,int par)
  // {remap_vector_internal((char*)out,(char*)in,sizeof(T),loceo_of_Lebeo[par],loc_volh);}
  // template <class T> void remap_Lebeo_to_loceo_vector(T **out,T **in)
  // {for(int eo=0;eo<2;eo++) remap_Leb_ev_or_od_to_loc_vector(out[eo],in[eo],eo);}
  // template <class T> void remap_loc_ev_or_od_to_Leb_vector(T *out,T *in,int par)
  // {remap_vector_internal((char*)out,(char*)in,sizeof(T),Lebeo_of_loceo[par],loc_volh);}
  // template <class T> void remap_loceo_to_Lebeo_vector(T **out,T **in)
  // {for(int eo=0;eo<2;eo++) remap_loc_ev_or_od_to_Leb_vector(out[eo],in[eo],eo);}
  
  // void remap_loceo_conf_to_Lebeo_oct(oct_su3 *out,quad_su3 **in,int par);
}

#endif
