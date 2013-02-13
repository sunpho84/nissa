#ifndef _SPI_H
#define _SPI_H
int spi_start_communicating_lx_borders(spi_comm_t &a,void *vec,int nbytes_per_site);
void fill_ev_or_od_bord_with_spi_receiving_buf(void *vec,spi_comm_t &a,int nbytes_per_site);
void fill_lx_bord_with_spi_receiving_buf(void *vec,spi_comm_t &a,int nbytes_per_site);
void fill_spi_sending_buf_with_ev_or_od_vec(spi_comm_t &a,void *vec,int nbytes_per_site,int eo);
void fill_spi_sending_buf_with_lx_vec(spi_comm_t &a,void *vec,int nbytes_per_site);
void get_spi_coord();
void init_spi();
void set_eo_spi_comm(spi_comm_t &in,int nbytes_per_site) ;
void set_lx_or_eo_spi_comm(spi_comm_t &in,int nbytes_per_site,int lx_eo);
void set_lx_spi_comm(spi_comm_t &in,int nbytes_per_site) ;
void set_spi_comm(spi_comm_t &in,int buf_size,int *payload_address_offset,int *message_length,MUHWI_Destination *dest,int *rec_payload_offset);
void set_spi_geometry();
void set_spi_neighbours();
void spi_communicate_ev_or_od_borders(void *vec,spi_comm_t &a,int nbytes_per_site,int eo);
void spi_communicate_lx_borders(void *vec,spi_comm_t &a,int nbytes_per_site);
void spi_comm_wait(spi_comm_t &in);
void spi_finish_communicating_ev_or_od_borders(void *vec,spi_comm_t &a,int nbytes_per_site);
void spi_finish_communicating_lx_borders(void *vec,spi_comm_t &a,int nbytes_per_site);
void spi_global_barrier();
void spi_start_comm(spi_comm_t &in);
int spi_start_communicating_ev_or_od_borders(spi_comm_t &a,void *vec,int nbytes_per_site,int eo);
void unset_spi_comm(spi_comm_t &in);
#endif
