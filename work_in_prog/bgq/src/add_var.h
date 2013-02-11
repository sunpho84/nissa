#ifndef _ADD_VAR_H
#define _ADD_VAR_H

#include <spi/include/kernel/MU.h>
#include <spi/include/mu/InjFifo.h>

//type to hold the 5D coordinates
typedef uint8_t coords_5D[5];

//spi rank coordinates
coords_5D spi_rank_coord;

//neighbours int the 4 dirs
MUHWI_Destination_t spi_neigh[2][4];

//spi fifo
uint64_t *spi_fifo[8];
MUSPI_InjFifoSubGroup_t spi_fifo_sg_ptr;

//spi descriptors
char *spi_descriptors;

#endif
