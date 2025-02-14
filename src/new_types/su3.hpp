#ifndef _SU3_HPP
#define _SU3_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "complex.hpp"
#include "dirac.hpp"
#include "spin.hpp"

#define NCOL 3

#define UNROLL_FOR_ALL_COLS(IC)			\
  UNROLL_FOR(IC,0,NCOL)

namespace nissa
{
  template <typename F>
  using Color=
    Complex<F>[NCOL];
  
  using color=
    Color<double>;
  
  typedef complex color2[2];
  typedef color2 su2[2];
  
  typedef color halfspincolor[NDIRAC/2];
  typedef halfspincolor color_halfspincolor[NCOL];
  typedef color_halfspincolor halfspincolor_halfspincolor[NDIRAC/2];
  
  typedef color spincolor[NDIRAC];
  typedef spin colorspin[NCOL];
  
  typedef colorspin spincolorspin[NDIRAC];
  typedef spincolorspin colorspincolorspin[NCOL];
  
  typedef spinspin colorspinspin[NCOL];
  
  template <typename F>
  using Su3=
    Color<F>[NCOL];
  
  using su3=
    Su3<double>;
  
  template <typename F>
  using QuadSu3=
    Su3<F>[NDIM];
  
  using quad_su3=
    QuadSu3<double>;
  
  typedef su3 quad_su3[NDIM];
  typedef su3 oct_su3[2*NDIM];
  
  typedef colorspinspin su3spinspin[NCOL];
  
  typedef su3 as2t_su3[NDIM*(NDIM+1)/2];
  typedef su3 clover_term_t[4];
  typedef halfspincolor_halfspincolor inv_clover_term_t[2];
}

#endif
