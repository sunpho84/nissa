#ifndef _MPITYPES_HPP
#define _MPITYPES_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#ifdef USE_MPI
# include <mpi.h>
#endif
#include <expr/comp.hpp>

namespace nissa
{
  DECLARE_DYNAMIC_COMP(MpiRankCoord);
  
  DECLARE_DYNAMIC_COMP(MpiRank);
  
  /// Master rank
  constexpr MpiRank masterRank=0;
  
  PROVIDE_RESOURCE(thisRank,MpiRank);
  
  /// Check if this is the master rank
  INLINE_FUNCTION
  bool isMasterRank()
  {
    return thisRank==masterRank;
  }
  
  /// Barrier across all ranks
  inline void mpiRanksBarrier()
  {
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  
  /// Abort execution
  __attribute__((noreturn))
  inline void mpiAbort(const int& err) 
  {
    printf("on rank %ld aborting\n",thisRank());
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(0);
  }
  
  /// Broadcast a passed variable
  template <TriviallyCopyable T>
  void mpiBcast(T& t,
		const MpiRank& rank=masterRank)
    requires(not std::is_pointer_v<T>)
  {
#ifdef USE_MPI
    MPI_Bcast(&t,sizeof(T),MPI_CHAR,rank(),MPI_COMM_WORLD);
#endif
  }
  
  /// Returns the passed variable after broadcast
  template <TriviallyCopyable T>
  T getMpiBcast(const T& t,
	     const MpiRank& rank=masterRank)
    requires(not std::is_pointer_v<T>)
  {
    T temp;
    mpiBcast(temp,rank);
    
    return temp;
  }
  
  /// Return the name of the processor
  inline std::string mpiGetProcessorName()
  {
#ifdef USE_MPI
    int resultlen=MPI_MAX_PROCESSOR_NAME;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(name,&resultlen);
    
    return name;
#else
    return "proc";
#endif
  }
  
}

#endif
