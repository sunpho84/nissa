#ifndef _SYSTEMRANDOM_HPP
#define _SYSTEMRANDOM_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <fcntl.h>
#include <mpi.h>
#include <unistd.h>

#include <base/debug.hpp>
#include <routines/rank.hpp>

namespace nissa
{
  /// Read from /dev/urandom
  template <typename T>
  void get_system_random(T &t)
  {
    const int size=sizeof(T);
    
    if(is_master_rank())
      {
	const char path[]="/dev/urandom";
	int fd=open(path,O_RDONLY);
	if(fd==-1) crash("Opening %s",path);
	
	int rc=read(fd,&t,size);
        if(rc!=size) crash("reading %zu bytes from %s, obtained: %d",size,path,rc);
	if(close(fd)==-1) crash("Closing %s",path);
      }
    MPI_Bcast(&t,size,MPI_CHAR,master_rank,MPI_COMM_WORLD);
  }
}

#endif
