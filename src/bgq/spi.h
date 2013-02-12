#ifndef _SPI_H
#define _SPI_H
void allocate_spi_comm(spi_comm_t &in,int buf_size);
void get_spi_coord();
void init_spi();
void set_lx_spi_comm(spi_comm_t &in,int nbytes_per_site);
void set_spi_geometry();
void set_spi_neighbours();
void spi_comm_wait(spi_comm_t &in);
void spi_global_barrier();
void spi_start_comm(spi_comm_t &in);
void unset_spi_comm(spi_comm_t &in);
#endif
