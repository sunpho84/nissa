#ifndef COMMUNICATE_HPP
#define COMMUNICATE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <stdint.h>
#include <mpi.h>

#include <base/memoryManager.hpp>
#include <routines/ios.hpp>
#include <routines/mpiRoutines.hpp>

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

namespace nissa
{
  /// Handle to manage async communications
  struct CommHandle
  {
    /// Requests
    MpiRequests requests;
    
    /// Buffer
    char* buf;
    
    /// Buffer size
    size_t n;
    
    MemoryType CurMT;
    
    /// Allocates the buffer
    void allocate(const size_t& _n,
		  const MemoryType& MT)
    {
      if(buf)
	CRASH("Already allocated");
      
      n=_n;
      CurMT=MT;
      
      buf=getMemoryManager(CurMT)->provide<char>(n);
    }
    
    /// Free the buffers
    void free()
    {
      if(not buf)
	CRASH("Not allocated");
      
      getMemoryManager(CurMT)->release(buf);
      n=0;
    }
    
    /// Move the buffers to the given memory space
    void makeAvailableOn(const MemoryType& OMT)
    {
      if(OMT!=CurMT)
	{
	  CommHandle tmp;
	  tmp.allocate(n,OMT);
	  nissa::memcpy(OMT,CurMT,tmp.buf,buf,n);
	  std::swap(*this,tmp);
	}
    }
    
    /// Default constructor
    CommHandle() :
      buf(nullptr),
      n(0)
    {
    }
    
    /// Construct with the knowledge of n and MT
    CommHandle(const size_t& n,
	       const MemoryType& MT)
    {
      allocate(n,MT);
    }
    
    /// Move constructor
    CommHandle(CommHandle&& oth) :
      requests(std::move(oth.requests)),
      buf(oth.buf),
      n(oth.n),
      CurMT(oth.CurMT)
    {
      oth.buf=nullptr;
    }
    
    /// Move assign
    CommHandle& operator=(CommHandle&& oth)
    {
      std::swap(buf,oth.buf);
      std::swap(CurMT,oth.CurMT);
      std::swap(n,oth.n);
      
      return *this;
    }
    
    /// Destructor
    ~CommHandle()
    {
      if(requests.size())
	mpiWaitAll(requests);
      
      if(buf)
	free();
    }
  };
  
  /// Handle to communicate in and out communications
  struct CommHandles
  {
    /// Sending part
    CommHandle send;
    
    /// Receiving part
    CommHandle recv;
  };
}

#endif
