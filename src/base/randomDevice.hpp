#ifndef _RANDOMDEVICE_HPP
#define _RANDOMDEVICE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file base/randomDevice.hpp

#include <cstdio>
#include <fcntl.h>
#include <unistd.h>

#include <base/debug.hpp>
#include <metaprogramming/concepts.hpp>
#include <routines/mpiRoutines.hpp>

namespace nissa
{
  /// Read from /dev/urandom
  template <TriviallyCopyable T>
  void get_system_random(T& t)
  {
    const int size=sizeof(T);
    
    if(isMasterRank())
      {
	constexpr char path[]=
	  "/dev/urandom";
	
	const int fd=
	  open(path,O_RDONLY);
	if(fd==-1) crash("Opening %s",path);
	
	
        if(const int rc=read(fd,&t,size);rc!=size)
	  crash("reading %zu bytes from %s, obtained: %d",size,path,rc);
	
	if(close(fd)==-1) crash("Closing %s",path);
    }
    
    mpiBcast(t);
  }
}

#endif
