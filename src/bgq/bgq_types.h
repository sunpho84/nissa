#ifndef _BGQ_TYPES_H
#define _BGQ_TYPES_H

#include <spi/include/kernel/MU.h>

//////////////// new types /////////////////

//type to hold the 5D coordinates
typedef uint8_t coords_5D[5];

//structure used to hold spi buffers 
struct spi_comm_t
{
  //communication in progress
  int comm_in_prog;
  //size of the buffers, buffers
  uint64_t buf_size;
  char *send_buf,*recv_buf;
  //counter for received bytes
  volatile uint64_t recv_counter;
  //descriptors
  MUHWI_Descriptor_t *descriptors;
  //bat
  MUSPI_BaseAddressTableSubGroup_t spi_bat_gr;
};

#endif
