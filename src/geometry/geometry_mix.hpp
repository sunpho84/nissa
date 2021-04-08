#ifndef _GEOMETRY_MIX_HPP
#define _GEOMETRY_MIX_HPP

#include <unistd.h>

#include "geometry_eo.hpp"
#include "geometry_Leb.hpp"

#include "new_types/su3.hpp"

namespace nissa
{
  void paste_eo_parts_into_lx_vector_internal(void *out_lx,eo_ptr<void> in_eo,int64_t bps);
  
  template <class T> void paste_eo_parts_into_lx_vector(T *out_lx,eo_ptr<T> in_eo)
  {
    paste_eo_parts_into_lx_vector_internal(out_lx,{in_eo[EVN],in_eo[ODD]},sizeof(T));
  }
  
  /////////////////////////////////////////////////////////////////
  
  void split_lx_vector_into_eo_parts_internal(eo_ptr<void> out_eo,void *in_lx,int64_t bps);
  
  template <class T> void split_lx_vector_into_eo_parts(eo_ptr<T> out_eo,T *in_lx)
  {
    split_lx_vector_into_eo_parts_internal({out_eo[EVN],out_eo[ODD]},in_lx,sizeof(T));
  }
  
  /////////////////////////////////////////////////////////////////
  
  void get_evn_or_odd_part_of_lx_vector_internal(void *out_ev_or_od,void *in_lx,int64_t bps,int par);
  
  template <class T> void get_evn_or_odd_part_of_lx_vector(T *out_eo,T *in_lx,int par)
  {
    get_evn_or_odd_part_of_lx_vector_internal(out_eo,in_lx,sizeof(T),par);
  }
  
  /////////////////////////////////////////////////////////////////
  
  void remap_vector_internal(char *out,char *in,int64_t bps,int *dest_of_source,int length);
  
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
