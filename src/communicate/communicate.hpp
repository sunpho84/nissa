#ifndef COMMUNICATE_HPP
#define COMMUNICATE_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>
#include <stdio.h>
#include <stdint.h>

#define NISSA_DEFAULT_WARN_IF_NOT_COMMUNICATED 0
#define NISSA_DEFAULT_USE_ASYNC_COMMUNICATIONS 1

/*
  Order in memory of borders for a 3^4 lattice.
  Border contains all sites having a single coordinate equal to -1, or equal to loc_size.
  The number of such sites is bord_vol, and they are divided in two groups.
  First group is relative to the "backward" border (all sites having just one local coordinate -1),
  second group is relative to "forward" border (all sites having just one coordinate loc_size).
  
  Each group contains the borders of all the 4 directions, each only if really parallelized, in the order (txyz).
  For a 3x3 system the borders are numbered as following:
  
      6 7 8
     -------         ^
  5 | X X X | B      |
  4 | X X X | A      0
  3 | X X X | 9      |
     -------         X---1--->
      0 1 2
  
  This is the shape and ordering of the border in the memory, for a 3^4 lattice
 _____________________________________________________________________________________________________________________
|___________________________________________________dir____=_____0____________________________________________________|
|_______________x__=__0_______________|||_______________x__=__1_______________|||_______________x__=__2_______________|
|____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|
| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
 ---------------------------------------------------------------------------------------------------------------------
 
 _____________________________________________________________________________________________________________________
|___________________________________________________dir____=_____1____________________________________________________|
|_______________t__=__0_______________|||_______________t__=__1_______________|||_______________t__=__2_______________|
|____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|
| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
 ---------------------------------------------------------------------------------------------------------------------
 
 _____________________________________________________________________________________________________________________
|___________________________________________________dir____=_____2____________________________________________________|
|_______________t__=__0_______________|||_______________t__=__1_______________|||_______________t__=__2_______________|
|____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|
| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
 ---------------------------------------------------------------------------------------------------------------------
 
 _____________________________________________________________________________________________________________________
|___________________________________________________dir____=_____3____________________________________________________|
|_______________t__=__0_______________|||_______________t__=__1_______________|||_______________t__=__2_______________|
|____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|
| y | y | y || y | y | y || y | y | y ||| y | y | y || y | y | y || y | y | y ||| y | y | y || y | y | y || y | y | y |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
 ---------------------------------------------------------------------------------------------------------------------
 
*/

#ifndef EXTERN_COMMUNICATE
 #define EXTERN_COMMUNICATE extern
#endif

#ifdef SPI
 #include <spi/include/kernel/MU.h>
#endif

namespace nissa
{
  
#ifdef USE_MPI
  //out and in buffer
  struct comm_t
  {
    //bgq specific structures, in alternative to ordinary MPI
#ifdef SPI
    //descriptors
    MUHWI_Descriptor_t *descriptors;
    MUHWI_Destination spi_dest[8];
#else
    //destinations and source ranks
    int send_rank[8],recv_rank[8];
    //requests and message
    MPI_Request requests[16];
    int nrequest,imessage;
#endif
    
    //communication in progress
    int comm_in_prog;
    //local size
    uint64_t nbytes_per_site;
    //size of the message
    uint64_t tot_mess_size;
    //offsets
    int send_offset[8],message_length[8],recv_offset[8];
    
    //constructor
    bool initialized;
    comm_t(){initialized=false;}
  };
#endif
  
  EXTERN_COMMUNICATE int ncomm_allocated;
  EXTERN_COMMUNICATE int comm_in_prog;
  EXTERN_COMMUNICATE int warn_if_not_communicated;
  EXTERN_COMMUNICATE int use_async_communications;
  
  //buffers
  EXTERN_COMMUNICATE uint64_t recv_buf_size,send_buf_size;
  EXTERN_COMMUNICATE char *recv_buf,*send_buf;
  
#define DEFINE_COMM(T) EXTERN_COMMUNICATE comm_t NAME3(lx,T,comm),NAME3(eo,T,comm)
  
  DEFINE_COMM(spin);
  DEFINE_COMM(spin1field);
  DEFINE_COMM(color);
  DEFINE_COMM(spincolor);
  DEFINE_COMM(spincolor_128);
  DEFINE_COMM(halfspincolor);
  DEFINE_COMM(colorspinspin);
  DEFINE_COMM(spinspin);
  DEFINE_COMM(su3spinspin);
  DEFINE_COMM(su3);
  DEFINE_COMM(as2t_su3);
  DEFINE_COMM(quad_su3);
  DEFINE_COMM(oct_su3);
  DEFINE_COMM(single_color);
  DEFINE_COMM(single_halfspincolor);
  DEFINE_COMM(single_quad_su3);
}

#undef EXTERN_COMMUNICATE

#endif
