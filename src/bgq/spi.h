#ifndef _SPI_H
#define _SPI_H

#include "../new_types/new_types_definitions.h"

void get_spi_coord();
void init_spi();
void set_spi_geometry();
void set_spi_neighbours();
void spi_global_barrier();
void spi_descriptor_setup(buffered_comm_t &in);
void spi_comm_start(buffered_comm_t &in,int *dir_comm,int tot_size);
void spi_comm_wait(buffered_comm_t &in);
void spi_descriptor_unset(buffered_comm_t &in);

#endif
