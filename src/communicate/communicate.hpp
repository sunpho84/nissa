#ifndef COMMUNICATE_H
#define COMMUNICATE_H

#include "borders.hpp"
#include "edges.hpp"

#include "base/global_variables.hpp"
#include "base/macros.hpp"

/*
  Order in memory of borders for a 3^4 lattice.
  Border contains all sites having a single coordinate equal to -1, or equal to loc_size.
  The number of such sites is bord_vol, and they are divided in two groups.
  First group is relative to the "backward" border (all sites having just one local coordinate -1),
  second group is relative to "forward" border (all sites having just one coordinate loc_size).
  
  Each group contains the borders of all the 4 directions, each only if really parallelized, in the order (txyz).
  For a 3x3 system the borders are numbered as following:
   
      6 7 8          
     -------         ^
  5 | X X X | B      |
  4 | X X X | A      0
  3 | X X X | 9      |
     -------         X---1--->
      0 1 2          
  
  This is the shape and ordering of the border in the memory, for a 3^4 lattice
 _____________________________________________________________________________________________________________________ 
|___________________________________________________dir____=_____0____________________________________________________|
|_______________x__=__0_______________|||_______________x__=__1_______________|||_______________x__=__2_______________|
|____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|
| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
 --------------------------------------------------------------------------------------------------------------------- 

 _____________________________________________________________________________________________________________________ 
|___________________________________________________dir____=_____1____________________________________________________|
|_______________t__=__0_______________|||_______________t__=__1_______________|||_______________t__=__2_______________|
|____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|
| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
 --------------------------------------------------------------------------------------------------------------------- 

 _____________________________________________________________________________________________________________________ 
|___________________________________________________dir____=_____2____________________________________________________|
|_______________t__=__0_______________|||_______________t__=__1_______________|||_______________t__=__2_______________|
|____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|
| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
 --------------------------------------------------------------------------------------------------------------------- 

 _____________________________________________________________________________________________________________________ 
|___________________________________________________dir____=_____3____________________________________________________|
|_______________t__=__0_______________|||_______________t__=__1_______________|||_______________t__=__2_______________|
|____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|
| y | y | y || y | y | y || y | y | y ||| y | y | y || y | y | y || y | y | y ||| y | y | y || y | y | y || y | y | y |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
 --------------------------------------------------------------------------------------------------------------------- 

*/

namespace nissa
{  
#define DEFINE_EO_BORDERS_ROUTINES(TYPE)				\
  inline void NAME3(communicate_ev_and_od,TYPE,borders)(TYPE **s)	\
  {communicate_ev_and_od_borders((void**)s,NAME3(lx,TYPE,comm));}	\
  inline void NAME3(communicate_ev_or_od,TYPE,borders)(TYPE *s,int eo)	\
  {communicate_ev_or_od_borders(s,NAME3(eo,TYPE,comm),eo);}		\
  inline void NAME3(start_communicating_ev_or_od,TYPE,borders)(TYPE *s,int eo) \
  {start_communicating_ev_or_od_borders(NAME3(eo,TYPE,comm),s,eo);}	\
  inline void NAME3(finish_communicating_ev_or_od,TYPE,borders)(TYPE *s) \
  {finish_communicating_ev_or_od_borders(s,NAME3(eo,TYPE,comm));}	\
  inline void NAME3(communicate_ev,TYPE,borders)(TYPE *s)		\
  {communicate_ev_or_od_borders(s,NAME3(eo,TYPE,comm),EVN);}		\
  inline void NAME3(communicate_od,TYPE,borders)(TYPE *s)		\
  {communicate_ev_or_od_borders(s,NAME3(eo,TYPE,comm),ODD);}
  
#define DEFINE_LX_BORDERS_ROUTINES(TYPE)			\
  inline void NAME3(communicate_lx,TYPE,borders)(TYPE *s)	\
  {communicate_lx_borders(s,NAME3(lx,TYPE,comm));}			\
  inline void NAME3(start_communicating_lx,TYPE,borders)(TYPE *s)	\
  {start_communicating_lx_borders(NAME3(lx,TYPE,comm),s);}		\
  inline void NAME3(finish_communicating_lx,TYPE,borders)(TYPE *s)	\
  {finish_communicating_lx_borders(s,NAME3(lx,TYPE,comm));}
  
#define DEFINE_BORDERS_ROUTINES(TYPE)		\
  DEFINE_LX_BORDERS_ROUTINES(TYPE)		\
  DEFINE_EO_BORDERS_ROUTINES(TYPE)
  
  DEFINE_BORDERS_ROUTINES(spin)
  DEFINE_BORDERS_ROUTINES(color)
  DEFINE_BORDERS_ROUTINES(spincolor)
  DEFINE_BORDERS_ROUTINES(spincolor_128)
  DEFINE_BORDERS_ROUTINES(halfspincolor)
  DEFINE_BORDERS_ROUTINES(colorspinspin)
  DEFINE_BORDERS_ROUTINES(spinspin)
  DEFINE_BORDERS_ROUTINES(su3spinspin)
  DEFINE_BORDERS_ROUTINES(su3)
  DEFINE_BORDERS_ROUTINES(quad_su3)
  DEFINE_BORDERS_ROUTINES(single_color)
  DEFINE_BORDERS_ROUTINES(single_quad_su3)
  DEFINE_BORDERS_ROUTINES(single_halfspincolor)
}

#endif
