#ifndef _GEOMETRY_MIX_HPP
#define _GEOMETRY_MIX_HPP

#include <unistd.h>

namespace nissa
{
  void paste_eo_parts_into_lx_vector_internal(char *out_lx,char **in_eo,size_t bps);
  template <class T> void paste_eo_parts_into_lx_vector(T *out_lx,T **in_eo)
  {paste_eo_parts_into_lx_vector_internal((char*)out_lx,(char**)in_eo,sizeof(T));}
  
  void split_lx_vector_into_eo_parts_internal(char **out_eo,char *in_lx,size_t bps);
  template <class T> void split_lx_vector_into_eo_parts(T **out_eo,T *in_lx)
  {split_lx_vector_into_eo_parts_internal((char**)out_eo,(char*)in_lx,sizeof(T));}
  
  void get_evn_or_odd_part_of_lx_vector_internal(char *out_ev_or_od,char *in_lx,size_t bps,int par);
  template <class T> void get_evn_or_odd_part_of_lx_vector(T *out_eo,T *in_lx,int par)
  {get_evn_or_odd_part_of_lx_vector_internal((char*)out_eo,(char*)in_lx,sizeof(T),par);}
}

#endif
