#ifndef _BUFFERED_BORDERS_H
#define _BUFFERED_BORDERS_H

#include "../new_types/new_types_definitions.h"

void buffered_comm_start(buffered_comm_t *in,int *dir_comm=NULL,int tot_size=-1);
void buffered_comm_wait(buffered_comm_t *in);
void buffered_communicate_ev_and_od_borders(void **vec,buffered_comm_t *a,int nbytes_per_site);
void buffered_communicate_ev_or_od_borders(void *vec,buffered_comm_t *a,int nbytes_per_site,int eo);
void buffered_communicate_lx_borders(void *vec,buffered_comm_t *a,int nbytes_per_site);
void buffered_finish_communicating_ev_and_od_borders(void **vec,buffered_comm_t *a,int nbytes_per_site);
void buffered_start_communicating_ev_and_od_borders(buffered_comm_t *a,void **vec,size_t nbytes_per_site);
void fill_buffered_sending_buf_with_ev_and_od_vec(buffered_comm_t *a,void **vec,size_t nbytes_per_site);
void fill_buffered_sending_buf_with_ev_or_od_vec(buffered_comm_t *a,void *vec,size_t nbytes_per_site,int eo);
void fill_buffered_sending_buf_with_lx_vec(buffered_comm_t *a,void *vec,size_t nbytes_per_site);
void fill_ev_and_od_bord_with_buffered_receiving_buf(void **vec,buffered_comm_t *a,size_t nbytes_per_site);
void fill_ev_or_od_bord_with_buffered_receiving_buf(void *vec,buffered_comm_t *a,size_t nbytes_per_site);
void fill_lx_bord_with_buffered_receiving_buf(void *vec,buffered_comm_t *a,size_t nbytes_per_site);
void finish_buffered_communicating_ev_or_od_borders(void *vec,buffered_comm_t *a,int nbytes_per_site);
void finish_buffered_communicating_lx_borders(void *vec,buffered_comm_t *a,int nbytes_per_site);
void set_buffered_comm(buffered_comm_t &in);
void set_eo_buffered_comm(buffered_comm_t &in,int nbytes_per_site);
void set_lx_buffered_comm(buffered_comm_t &in,int nbytes_per_site);
void set_lx_or_eo_buffered_comm(buffered_comm_t &in,int nbytes_per_site,int lx_eo);
void start_buffered_communicating_ev_or_od_borders(buffered_comm_t *a,void *vec,int nbytes_per_site,int eo);
void start_buffered_communicating_lx_borders(buffered_comm_t *a,void *vec,int nbytes_per_site);
void unset_buffered_comm(buffered_comm_t *in);
#endif
