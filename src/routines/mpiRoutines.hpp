#ifndef _MPIROUTINES_HPP
#define _MPIROUTINES_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#ifdef USE_MPI
# include <mpi.h>
#endif

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
  
  /// Exec only on master rank and gets the result
  template <typename F,
	    typename...Args>
  INLINE_FUNCTION
  auto getExecOnMasterRank(const F& f,
			   Args&&...args)
    -> decltype(f(std::forward<Args>(args)...))
  {
    return mpiGetBcast(onMasterRank(f,std::forward<Args>(args)...));
  }
  
  /// Initialize mpi
  inline void mpiInit(int narg,
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
  
  /// Finalize mpi
  inline void mpiFinalize()
  {
#ifdef USE_MPI
    MPI_Finalize();
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
  
  /// Decrypt the MPI error
  inline void internalDecryptMpiError(const int& line,
				      const char *file,
				      const int& rc,
				      const char *templ,...)
  {
#ifdef USE_MPI
    if(rc!=MPI_SUCCESS and isMasterRank())
      {
	char err[1024];
	int len=1024;
	MPI_Error_string(rc,err,&len);
	
	va_list ap;
	va_start(ap,templ);
	char mess[1024];
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	internalCrash(line,file,"%s, MPI raised error: %s",mess,err);
      }
#endif
  }
  
#define decryptMpiError(...)					\
  internalDecryptMpiError(__LINE__,__FILE__,__VA_ARGS__)
}

#endif
