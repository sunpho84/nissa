#ifndef _ADD_VAR_H
#define _ADD_VAR_H

#include <spi/include/kernel/MU.h>
#include <spi/include/mu/InjFifo.h>
#include <spi/include/mu/GIBarrier.h>

//////////////// new types ///////////////

//type to hold the 5D coordinates
typedef uint8_t coords_5D[5];

//structure used to hold spi buffers 
struct spi_comm_t
{
  //size of the buffers, buffers
  uint64_t buf_size;
  char *send_buf,*recv_buf;
  //counter for received bytes
  volatile uint64_t recv_counter;
  //descriptors
  MUHWI_Descriptor_t *descriptors;
};

/////////////// new vars /////////////////

//spi rank coordinates
coords_5D spi_rank_coord;

//neighbours int the 4 dirs
MUHWI_Destination_t spi_neigh[2][4];

//spi fifo
uint64_t *spi_fifo[8];
MUSPI_InjFifoSubGroup_t spi_fifo_sg_ptr;

//spi barrier
MUSPI_GIBarrier_t spi_barrier;

#endif
