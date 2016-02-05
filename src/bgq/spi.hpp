#ifndef _SPI_HPP
#define _SPI_HPP

#ifdef SPI
 #include <firmware/include/personality.h>
 #include <spi/include/mu/Addressing_inlines.h>
 #include <spi/include/mu/Descriptor_inlines.h>
 #include <spi/include/mu/GIBarrier.h>
 #include <spi/include/kernel/MU.h>
 #include <spi/include/kernel/location.h>
#endif

#include <stdint.h>

#include "communicate/communicate.hpp"

#ifndef EXTERN_SPI
 #define EXTERN_SPI extern
#endif

//spi fifo and counters for bytes
#define NSPI_FIFO 8

namespace nissa
{
  typedef uint8_t coords_5D[5];
  
  //flag to remember if spi has been initialized
  EXTERN_SPI int spi_inited;
  
  //spi rank coordinates
  EXTERN_SPI coords_5D spi_rank_coord;
#ifdef SPI
  EXTERN_SPI coords_5D spi_dir_is_torus,spi_dir_size;
#endif
  
  //destination coords
#ifdef SPI
  EXTERN_SPI MUHWI_Destination spi_dest[8];
#endif
  
  EXTERN_SPI coords_5D spi_dest_coord[8];
  
  //neighbours in the 4 dirs
#ifdef SPI
  EXTERN_SPI MUHWI_Destination_t spi_neigh[2][4];
  #endif
  
  EXTERN_SPI uint64_t *spi_fifo[NSPI_FIFO],spi_desc_count[NSPI_FIFO];
#ifdef SPI
  EXTERN_SPI MUSPI_InjFifoSubGroup_t spi_fifo_sg_ptr;
#endif
  EXTERN_SPI uint64_t spi_fifo_map[8];
  EXTERN_SPI uint8_t spi_hint_ABCD[8],spi_hint_E[8];
  
  //spi barrier
#ifdef SPI
  EXTERN_SPI MUSPI_GIBarrier_t spi_barrier;
  
  //bats
  EXTERN_SPI MUSPI_BaseAddressTableSubGroup_t spi_bat_gr;
#endif
  EXTERN_SPI uint32_t spi_bat_id[2];
  
  //physical address
  EXTERN_SPI uint64_t spi_send_buf_phys_addr;
  
  void get_spi_coord();
  void init_spi();
  void set_spi_geometry();
  void set_spi_neighbours();
  void spi_global_barrier();
  void spi_descriptor_setup(comm_t &in);
  void spi_comm_start(comm_t &in,int *dir_comm,int tot_size);
  void spi_comm_wait(comm_t &in);
  void spi_descriptor_unset(comm_t &in);
}

#undef EXTERN_SPI

#endif
