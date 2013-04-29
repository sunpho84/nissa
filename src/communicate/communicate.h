#ifndef COMMUNICATE_H
#define COMMUNICATE_H

#include "buffered_borders.h"
#include "edges.h"
#include "unbuffered_borders.h"

#include "../base/global_variables.h"

inline void communicate_lx_colorspinspin_borders(colorspinspin *s)
{buffered_communicate_lx_borders(s,buffered_lx_colorspinspin_comm);}
inline void buffered_start_communicating_lx_colorspinspin_borders(colorspinspin *s)
{buffered_start_communicating_lx_borders(buffered_lx_colorspinspin_comm,s);}
inline void buffered_finish_communicating_lx_colorspinspin_borders(colorspinspin *s)
{buffered_finish_communicating_lx_borders(s,buffered_lx_colorspinspin_comm);}

inline void communicate_lx_spincolor_borders(spincolor *s)
{buffered_communicate_lx_borders(s,buffered_lx_spincolor_comm);}
inline void buffered_start_communicating_lx_spincolor_borders(spincolor *s)
{buffered_start_communicating_lx_borders(buffered_lx_spincolor_comm,s);}
inline void buffered_finish_communicating_lx_spincolor_borders(spincolor *s)
{buffered_finish_communicating_lx_borders(s,buffered_lx_spincolor_comm);}

inline void communicate_lx_su3spinspin_borders(su3spinspin *s)
{buffered_communicate_lx_borders(s,buffered_lx_su3spinspin_comm);}

inline void communicate_lx_quad_su3_borders(quad_su3 *data)
{buffered_communicate_lx_borders(data,buffered_lx_quad_su3_comm);}

#endif
