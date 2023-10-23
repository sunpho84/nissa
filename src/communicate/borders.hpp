#ifndef _BORDERS_HPP
#define _BORDERS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "communicate/communicate.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/float_128.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
//   void comm_start(comm_t &comm,int *dir_comm=NULL,int tot_size=-1);
//   void communicate_ev_and_od_borders(eo_ptr<void> vec,comm_t &comm);
//   void communicate_Leb_ev_and_od_borders(eo_ptr<void> vec,comm_t &comm);
//   void communicate_ev_or_od_borders(void *vec,comm_t &comm,int eo);
//   void communicate_Leb_ev_or_od_borders(void *vec,comm_t &comm,int eo);
//   void communicate_lx_borders(void *vec,comm_t &comm);
//   void communicate_Leblx_borders(void *vec,comm_t &comm);
//   void comm_wait(comm_t &comm);
//   void finish_communicating_ev_and_od_borders(eo_ptr<void> vec,comm_t &comm);
//   void finish_communicating_Leb_ev_and_od_borders(eo_ptr<void> vec,comm_t &comm);
//   void finish_communicating_ev_or_od_borders(void *vec,comm_t &comm);
//   void finish_communicating_Leb_ev_or_od_borders(void *vec,comm_t &comm);
//   void finish_communicating_lx_borders(void *vec,comm_t &comm);
//   void finish_communicating_Leblx_borders(void *vec,comm_t &comm);
//   void start_communicating_ev_and_od_borders(comm_t &comm,eo_ptr<void> vec);
//   void start_communicating_Leb_ev_and_od_borders(comm_t &comm,eo_ptr<void> vec);
//   void start_communicating_ev_or_od_borders(comm_t &comm,void *vec,int eo);
//   void start_communicating_Leb_ev_or_od_borders(comm_t &comm,void *vec,int eo);
//   void start_communicating_lx_borders(comm_t &comm,void *vec);
//   void start_communicating_Leblx_borders(comm_t &comm,void *vec);
//   void fill_buffered_sending_buf_with_ev_and_od_vec(comm_t &comm,eo_ptr<void> vec);
//   void fill_buffered_sending_buf_with_ev_or_od_vec(comm_t &comm,void *vec,int eo);
//   void fill_buffered_sending_buf_with_lx_vec(comm_t &comm,void *vec);
//   void fill_ev_and_od_bord_with_buffered_receiving_buf(eo_ptr<void> vec,comm_t &comm);
//   void fill_ev_or_od_bord_with_buffered_receiving_buf(void *vec,comm_t &comm);
//   void fill_lx_bord_with_buffered_receiving_buf(void *vec,comm_t &comm);
//   void comm_set(comm_t &comm);
//   void set_eo_comm(comm_t &comm,int nbytes_per_site);
//   void set_lx_comm(comm_t &comm,int nbytes_per_site);
//   void set_lx_or_eo_comm(comm_t &comm,int lx_eo,int nbytes_per_site);
//   void comm_unset(comm_t &comm);
  
//   DEFINE_BORDERS_ROUTINES(spin)
//   DEFINE_BORDERS_ROUTINES(spin1field)
//   DEFINE_BORDERS_ROUTINES(color)
//   DEFINE_BORDERS_ROUTINES(spincolor)
//   DEFINE_BORDERS_ROUTINES(spincolor_128)
//   DEFINE_BORDERS_ROUTINES(halfspincolor)
//   DEFINE_BORDERS_ROUTINES(colorspinspin)
//   DEFINE_BORDERS_ROUTINES(spinspin)
//   DEFINE_BORDERS_ROUTINES(su3spinspin)
//   DEFINE_BORDERS_ROUTINES(su3)
//   DEFINE_BORDERS_ROUTINES(quad_su3)
//   DEFINE_BORDERS_ROUTINES(oct_su3)
//   DEFINE_BORDERS_ROUTINES(single_color)
//   DEFINE_BORDERS_ROUTINES(single_quad_su3)
//   DEFINE_BORDERS_ROUTINES(single_halfspincolor)
}

#endif
