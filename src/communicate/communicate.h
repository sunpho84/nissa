#ifndef COMMUNICATE_H
#define COMMUNICATE_H

#include "buffered_borders.h"
#include "edges.h"
#include "unbuffered_borders.h"

#include "../base/global_variables.h"

inline void communicate_ev_or_od_borders(color *s,int eo)
{buffered_communicate_ev_or_od_borders(s,buffered_eo_color_comm,eo);}
inline void buffered_start_communicating_ev_or_od_color_borders(color *s,int eo)
{buffered_start_communicating_ev_or_od_borders(buffered_eo_color_comm,s,eo);}
inline void buffered_finish_communicating_ev_or_od_color_borders(color *s)
{buffered_finish_communicating_ev_or_od_borders(s,buffered_eo_color_comm);}

inline void communicate_lx_spincolor_borders(spincolor *s)
{buffered_communicate_lx_borders(s,buffered_lx_spincolor_comm);}
inline void buffered_start_communicating_lx_spincolor_borders(spincolor *s)
{buffered_start_communicating_lx_borders(buffered_lx_spincolor_comm,s);}
inline void buffered_finish_communicating_lx_spincolor_borders(spincolor *s)
{buffered_finish_communicating_lx_borders(s,buffered_lx_spincolor_comm);}

inline void communicate_lx_colorspinspin_borders(colorspinspin *s)
{buffered_communicate_lx_borders(s,buffered_lx_colorspinspin_comm);}
inline void buffered_start_communicating_lx_colorspinspin_borders(colorspinspin *s)
{buffered_start_communicating_lx_borders(buffered_lx_colorspinspin_comm,s);}
inline void buffered_finish_communicating_lx_colorspinspin_borders(colorspinspin *s)
{buffered_finish_communicating_lx_borders(s,buffered_lx_colorspinspin_comm);}

inline void communicate_lx_su3spinspin_borders(su3spinspin *s)
{buffered_communicate_lx_borders(s,buffered_lx_su3spinspin_comm);}

inline void communicate_lx_quad_su3_borders(quad_su3 *data)
{buffered_communicate_lx_borders(data,buffered_lx_quad_su3_comm);}

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

#endif
