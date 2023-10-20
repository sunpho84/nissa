#ifndef _SU3_HPP
#define _SU3_HPP

#include "complex.hpp"
#include "dirac.hpp"
#include "spin.hpp"

#if NCOL == 3
 #define CRASH_IF_NOT_3COL()
#else
 #define CRASH_IF_NOT_3COL() crash("ncol == %d, expected 3",NCOL)
#endif

namespace nissa
{
  typedef complex color0[NCOL];
  
  typedef complex color2[2];
  typedef color2 su2[2];
  
  typedef color0 halfspincolor[NDIRAC/2];
  typedef halfspincolor color_halfspincolor[NCOL];
  typedef color_halfspincolor halfspincolor_halfspincolor[NDIRAC/2];
  
  typedef color0 spincolor[NDIRAC];
  typedef spin0 colorspin[NCOL];
  
  typedef colorspin spincolorspin[NDIRAC];
  typedef spincolorspin colorspincolorspin[NCOL];
  
  typedef spinspin colorspinspin[NCOL];
  
  typedef color0 su3[NCOL];
  typedef su3 quad_su3[NDIM];
  typedef su3 oct_su3[2*NDIM];
  
  typedef colorspinspin su3spinspin[NCOL];
  
  typedef su3 as2t_su3[NDIM*(NDIM+1)/2];
  typedef su3 clover_term_t[4];
  typedef halfspincolor_halfspincolor inv_clover_term_t[2];
  
  typedef single_complex single_color[NCOL];
  typedef single_color single_su3[NCOL];
  typedef single_color single_halfspincolor[2];
  typedef single_color single_spincolor[NDIRAC];
  typedef single_su3 single_quad_su3[NDIRAC];
}

#endif
