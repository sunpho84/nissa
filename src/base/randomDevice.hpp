#ifndef _RANDOMDEVICE_HPP
#define _RANDOMDEVICE_HPP

#include <cstdio>
#include <fcntl.h>
#include <unistd.h>

#include <base/debug.hpp>
#include <routines/mpi_routines.hpp>

namespace nissa
{
  //read from /dev/urandom
  template <typename T>
  void get_system_random(T& t)
  {
    const int size=sizeof(T);
    
    if(is_master_rank())
      {
	const char path[]="/dev/urandom";
	const int fd=open(path,O_RDONLY);
	if(fd==-1) crash("Opening %s",path);
	
	const int rc=read(fd,&t,size);
        if(rc!=size) crash("reading %zu bytes from %s, obtained: %d",size,path,rc);
	if(close(fd)==-1) crash("Closing %s",path);
    }
    
    MPI_Bcast(&t,size,MPI_CHAR,master_rank,MPI_COMM_WORLD);
  }
}

#endif
