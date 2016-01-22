#ifndef _SPI_HPP
#define _SPI_HPP

#ifdef SPI
 #include <spi/include/kernel/MU.h>
#endif

#include <stdint.h>

#include "communicate/communicate.hpp"

namespace nissa
{
  typedef uint8_t coords_5D[5];
  
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

#endif
