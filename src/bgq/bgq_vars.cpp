#include <spi/include/kernel/MU.h>
#include <spi/include/mu/InjFifo.h>
#include <spi/include/mu/GIBarrier.h>

#include "bgq_types.h"

//to be used to declare extern
#ifndef EXTERN
#define EXTERN
#endif

//spi rank coordinates
EXTERN coords_5D spi_rank_coord;

//neighbours int the 4 dirs
EXTERN MUHWI_Destination_t spi_neigh[2][4];

//spi fifo and counters for bytes
EXTERN uint64_t *spi_fifo[8],spi_desc_count[8];
EXTERN MUSPI_InjFifoSubGroup_t spi_fifo_sg_ptr;

//spi barrier
EXTERN MUSPI_GIBarrier_t spi_barrier;

#undef EXTERN
