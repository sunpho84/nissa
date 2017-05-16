#ifndef _GEOMETRY_MIX_HPP
#define _GEOMETRY_MIX_HPP

#include <unistd.h>

#include "geometry_Leb.hpp"

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
  
  void remap_eo_vector_internal(char **out,char **in,size_t bps,int **remap);
  template <class T> void remap_loceo_to_Lebeo_vector(T **out,T **in)
  {remap_eo_vector_internal((char**)out,(char**)in,sizeof(T),Lebeo_of_loceo);}
  template <class T> void remap_Lebeo_to_loceo_vector(T **out,T **in)
  {remap_eo_vector_internal((char**)out,(char**)in,sizeof(T),loceo_of_Lebeo);}
  
  void remap_lx_vector_internal(char *out,char *in,size_t bps,int *remap);
  template <class T> void remap_loclx_to_Leblx_vector(T *out,T *in)
  {remap_lx_vector_internal((char*)out,(char*)in,sizeof(T),Leblx_of_loclx);}
  template <class T> void remap_Leblx_to_loclx_vector(T *out,T *in)
  {remap_lx_vector_internal((char*)out,(char*)in,sizeof(T),loclx_of_Leblx);}
}

#endif
