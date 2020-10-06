#ifndef _GEOMETRY_LEB_HPP
#define _GEOMETRY_LEB_HPP

#ifndef EXTERN_GEOMETRY_LEB
 #define EXTERN_GEOMETRY_LEB extern
#endif

#define NISSA_DEFAULT_USE_LEB_GEOM 0

#include "geometry/geometry_lx.hpp"

#include <vector>

namespace nissa
{
  EXTERN_GEOMETRY_LEB int *Leblx_of_loclx;
  EXTERN_GEOMETRY_LEB int *loclx_of_Leblx;
  EXTERN_GEOMETRY_LEB int *Leblx_parity;
  EXTERN_GEOMETRY_LEB coords *Leblx_neighup;
  EXTERN_GEOMETRY_LEB coords *Leblx_neighdw;
  EXTERN_GEOMETRY_LEB int *surfLeblx_of_bordLeblx;
  
  EXTERN_GEOMETRY_LEB int *Lebeo_of_loceo[2];
  EXTERN_GEOMETRY_LEB int *loceo_of_Lebeo[2];
  EXTERN_GEOMETRY_LEB coords *Lebeo_neighup[2];
  EXTERN_GEOMETRY_LEB coords *Lebeo_neighdw[2];
  
  EXTERN_GEOMETRY_LEB int Leb_geom_inited;
  EXTERN_GEOMETRY_LEB int use_Leb_geom;
  
  typedef std::vector<std::vector<int> > Leb_factors_t;
  
  void set_Leb_geometry();
  void unset_Leb_geometry();
}

#endif
