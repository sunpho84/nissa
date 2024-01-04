#ifndef _MPIROUTINES_HPP
#define _MPIROUTINES_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <mpi.h>
#include <algorithm>
#include <concepts>
#include <vector>

#include <base/universe.hpp>
#include <expr/comp.hpp>
#include <routines/mpiRank.hpp>

namespace nissa
{
  /// Mpi Rank coordinates
  using MpiRankCoords=
    Coords<MpiRankCoord>;
  
  PROVIDE_RESOURCE(nRanks,MpiRank);
  PROVIDE_RESOURCE(nRanksPerDir,MpiRankCoords);
  PROVIDE_RESOURCE(isDirParallel,Coords<bool>);
  PROVIDE_RESOURCE(thisRankCoords,MpiRankCoords);
  PROVIDE_RESOURCE(neighRanks,StackTens<CompsList<Ori,Dir>,MpiRank>);
  
  /// Gets the number of ranks
  INLINE_FUNCTION
  MpiRank getMpiNRanks()
  {
    int res;
    
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&res);
#else
    res=1;
#endif
    
    return res;
  }
  
  /// Gets present rank
  INLINE_FUNCTION
  MpiRank getMpiRank()
  {
    int res;
    
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&res);
#else
    res=masterRank();
#endif
    
    return res;
  }
  
  /// Send and receive a trivial type
  template <TriviallyCopyable T>
  INLINE_FUNCTION
  T mpiSendrecv(const MpiRank& rankTo,
		const T& toSend,
		const MpiRank& rankFr)
  {
    T toRecv;
    
#ifdef USE_MPI
    MPI_Sendrecv(&toSend,sizeof(T),MPI_CHAR,rankTo(),0,
		 &toRecv,sizeof(T),MPI_CHAR,rankFr(),0,
		 MPI_COMM_WORLD,MPI_STATUS_IGNORE);
#else
    toRecv=toSend;
#endif
    
    return toRecv;
  }
  
  /// Specialization of sendrecv for std::vector
  template <TriviallyCopyable T>
  INLINE_FUNCTION
  std::vector<T> mpiSendrecv(const MpiRank& rankTo,
			     const std::vector<T>& toSend,
			     const MpiRank& rankFr)
  {
    const size_t nToSend=
      toSend.size();
    
    const size_t nToRecv=
      mpiSendrecv(rankTo,nToSend,rankFr);
    
    std::vector<T> toRecv(nToRecv);
    
#ifdef USE_MPI
    MPI_Sendrecv(&toSend[0],nToSend*sizeof(T),MPI_CHAR,rankTo(),0,
		 &toRecv[0],nToRecv*sizeof(T),MPI_CHAR,rankFr(),0,
		 MPI_COMM_WORLD,MPI_STATUS_IGNORE);
#else
    toRecv=toSend;
#endif
    
    return toRecv;
  }
  
  /// Exec only on master rank
  template <typename F,
	    typename...Args>
  INLINE_FUNCTION
  auto onMasterRank(const F& f,
		    Args&&...args)
    -> decltype(f(std::forward<Args>(args)...))
  {
    if(isMasterRank())
      return f(std::forward<Args>(args)...);
    else
      return {};
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
  
  /// Exec only on master rank and gets the result
  template <typename F,
	    typename...Args>
  INLINE_FUNCTION
  auto getExecOnMasterRank(const F& f,
			   Args&&...args)
    -> decltype(f(std::forward<Args>(args)...))
  {
    return getMpiBcast(onMasterRank(f,std::forward<Args>(args)...));
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
  
  /// Barrier across all ranks
  inline void mpiRanksBarrier()
  {
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  
  /// Initialize mpi
  inline void mpiInitThread(int narg,
			    char **arg)
  {
#ifdef USE_MPI
# ifdef USE_OPENMP
    int provided;
    MPI_Init_thread(&narg,&arg,MPI_THREAD_SERIALIZED,&provided);
# else
    MPI_Init(&narg,&arg);
 #endif
#endif
  }
  
  /// Abort execution
  inline void mpiAbort(const int& err)
  {
#ifdef USE_MPI
    printf("on rank %ld aborting\n",thisRank());
    MPI_Abort(MPI_COMM_WORLD,0);
#else
    exit(0);
#endif
  }
  
#define PROVIDE_MPI_DATATYPE_OF(T,MPI_T)	\
  /*! MPI Datatype corresponding to T */	\
  inline MPI_Datatype _mpiDatatypeOf(T*)	\
  {						\
    return MPI_T;				\
  }
  
  PROVIDE_MPI_DATATYPE_OF(int64_t,MPI_LONG)
  PROVIDE_MPI_DATATYPE_OF(double,MPI_DOUBLE)
  
#undef PROVIDE_MPI_DATATYPE_OF
  
  /// Instantiates the correct datatype, given the type
  template <typename T>
  INLINE_FUNCTION
  MPI_Datatype mpiDatatypeOf()
  {
    return _mpiDatatypeOf((T*)nullptr);
  }
  
}

#endif
