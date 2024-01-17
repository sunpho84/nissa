#ifndef COMMUNICATE_HPP
#define COMMUNICATE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <stdint.h>
#include <mpi.h>
#include <vector>

#include <base/memoryManager.hpp>
#include <routines/ios.hpp>

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
  namespace resources
  {
    /// Receive buffer size
    inline uint64_t curRecvBufSize{0};
    
    /// Send buffer size
    inline uint64_t curSendBufSize{0};
    
    /// Receive buffer
    inline char* _recvBuf{nullptr};
    
    /// Send buffer
    inline char* _sendBuf{nullptr};
    
    /// Gets the buffer after checking if large enough or allocating if not allocated at all
    template <typename T=char>
    inline T* _getCommBuf(const char* name,
			  uint64_t& curBufSize,
			  char* & buf,
			  const uint64_t& reqBufSizeInUnitsOfT)
    {
      if(const uint64_t reqBufSize=reqBufSizeInUnitsOfT*sizeof(T);
	 curBufSize<reqBufSize)
	{
	  masterPrintf("%s buffer requested size %lu but allocated %lu reallocating\n",name,reqBufSize,curBufSize);
	  if(buf)
	    {
	      memoryManager<MemoryType::CPU>()->release(buf,true);
	      masterPrintf("Freed the previously allocated buffer\n");
	    }
	  curBufSize=reqBufSize;
	  buf=memoryManager<MemoryType::CPU>()->template provide<char>(reqBufSize);
	  masterPrintf("Allocated the new buffer\n");
	}
      
      return (T*)buf;
    }
  }
  
  /// Gets the receiving buffer
  template <typename T=char>
  inline T* getRecvBuf(const uint64_t& recvBufSizeInUnitsOfT)
  {
    using namespace resources;
    
    return _getCommBuf<T>("Receiving",curRecvBufSize,_recvBuf,recvBufSizeInUnitsOfT);
  }
  
  /// Gets the sending buffer
  template <typename T=char>
  inline T* getSendBuf(const uint64_t& sendBufSizeInUnitsOfT)
  {
    using namespace resources;
    
    return _getCommBuf<T>("Sending",curSendBufSize,_sendBuf,sendBufSizeInUnitsOfT);
  }
  
  /// Frees the communication buffers
  inline void freeCommunicationBuffers()
  {
    using namespace resources;
    
    for(auto* buf : {_recvBuf,_sendBuf})
      memoryManager<MemoryType::CPU>()->release(buf,true);
    
    for(auto* size : {&curRecvBufSize,&curSendBufSize})
      (*size)=0;
  }
}

#endif
