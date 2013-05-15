#ifndef _BUFFERED_BORDERS_H
#define _BUFFERED_BORDERS_H

#include "../new_types/new_types_definitions.h"

void comm_start(comm_t &comm,int *dir_comm=NULL,int tot_size=-1);
void communicate_ev_and_od_borders(void **vec,comm_t &comm);
void communicate_ev_or_od_borders(void *vec,comm_t &comm,int eo);
void communicate_lx_borders(void *vec,comm_t &comm);
void comm_wait(comm_t &comm);
void finish_communicating_ev_and_od_borders(void **vec,comm_t &comm);
void finish_communicating_ev_or_od_borders(void *vec,comm_t &comm);
void finish_communicating_lx_borders(void *vec,comm_t &comm);
void start_communicating_ev_and_od_borders(comm_t &comm,void **vec);
void start_communicating_ev_or_od_borders(comm_t &comm,void *vec,int eo);
void start_communicating_lx_borders(comm_t &comm,void *vec);
void fill_buffered_sending_buf_with_ev_and_od_vec(comm_t &comm,void **vec);
void fill_buffered_sending_buf_with_ev_or_od_vec(comm_t &comm,void *vec,int eo);
void fill_buffered_sending_buf_with_lx_vec(comm_t &comm,void *vec);
void fill_ev_and_od_bord_with_buffered_receiving_buf(void **vec,comm_t &comm);
void fill_ev_or_od_bord_with_buffered_receiving_buf(void *vec,comm_t &comm);
void fill_lx_bord_with_buffered_receiving_buf(void *vec,comm_t &comm);
void comm_set(comm_t &comm);
void set_eo_comm(comm_t &comm,int nbytes_per_site);
void set_lx_comm(comm_t &comm,int nbytes_per_site);
void set_lx_or_eo_comm(comm_t &comm,int lx_eo,int nbytes_per_site);
void comm_unset(comm_t &comm);

#endif
