#ifndef _MPIROUTINES_HPP
#define _MPIROUTINES_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#ifdef USE_MPI
# include <mpi.h>
#endif

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
  
  /////////////////////////////////////////////////////////////////
  
#ifdef USE_MPI
  
  //// MPI Datatype corresponding to n contiguous bytes
  inline MPI_Datatype mpiDatatypeOfNContiguousBytes(const size_t& n)
  {
    /// Lookup table for contiguous data of custom length
    static std::map<size_t,MPI_Datatype> table;
    
    if(auto probe=table.find(n);probe!=table.end())
      return probe->second;
    else
      {
	auto& res=table[n];
	MPI_Type_contiguous(n,MPI_CHAR,&res);
	MPI_Type_commit(&res);
	
	return res;
      }
  }
  
  namespace internal
  {
    //// MPI Datatype corresponding to T, if trivially copiable
    template <TriviallyCopyable T>
    inline MPI_Datatype _mpiDatatypeOf(T*)
    {
      return mpiDatatypeOfNContiguousBytes(sizeof(T));
    }
    
#define PROVIDE_MPI_DATATYPE_OF(T,MPI_T)	\
    /*! MPI Datatype corresponding to T */	\
    inline MPI_Datatype _mpiDatatypeOf(T*)	\
    {						\
      return MPI_T;				\
    }
    
    PROVIDE_MPI_DATATYPE_OF(char,MPI_CHAR)
    PROVIDE_MPI_DATATYPE_OF(int32_t,MPI_INT)
    PROVIDE_MPI_DATATYPE_OF(uint32_t,MPI_UINT32_T)
    PROVIDE_MPI_DATATYPE_OF(int64_t,MPI_LONG)
    PROVIDE_MPI_DATATYPE_OF(uint64_t,MPI_UINT64_T)
    PROVIDE_MPI_DATATYPE_OF(float,MPI_FLOAT)
    PROVIDE_MPI_DATATYPE_OF(double,MPI_DOUBLE)
    
#undef PROVIDE_MPI_DATATYPE_OF
  }
  
    /// Instantiates the correct datatype, given the type
    template <typename T>
    INLINE_FUNCTION
    MPI_Datatype mpiDatatypeOf()
    {
      return internal::_mpiDatatypeOf((T*)nullptr);
    }
  
#endif
  
  /////////////////////////////////////////////////////////////////
  
  /// Send and receive a trivial type
  template <TriviallyCopyable T>
  INLINE_FUNCTION
  T mpiSendrecv(const MpiRank& rankTo,
		const T& toSend,
		const MpiRank& rankFr)
  {
    T toRecv;
    
#ifdef USE_MPI
    MPI_Sendrecv(&toSend,1,mpiDatatypeOf<T>(),rankTo(),0,
		 &toRecv,1,mpiDatatypeOf<T>(),rankFr(),0,
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
    MPI_Sendrecv(&toSend[0],nToSend,mpiDatatypeOfNContiguousBytes(sizeof(T)),rankTo(),0,
		 &toRecv[0],nToRecv,mpiDatatypeOfNContiguousBytes(sizeof(T)),rankFr(),0,
		 MPI_COMM_WORLD,MPI_STATUS_IGNORE);
#else
    toRecv=toSend;
#endif
    
    return toRecv;
  }
  
  /// Decide the direction of the MPI communication
  enum class MpiSendOrRecvFlag{Send,Recv};
  
  using MpiRequest=
#ifdef USE_MPI
  MPI_Request
#else
  std::tuple<void*,int,MpiSendOrRecvFlag,size_t>
#endif
  ;
  
  /// List of requests to be accomplished
  using MpiRequests=
    std::vector<MpiRequest>;
  
  template <TriviallyCopyable T>
  INLINE_FUNCTION
  MpiRequest mpiISendOrRecv(MpiSendOrRecvFlag SR,
			    const MpiRank& othRank,
			    T* p,
			    const size_t& lengthInUnitsOfT,
			    const int& tag)
  {
#ifdef USE_MPI
    
# define ARGS p,lengthInUnitsOfT,mpiDatatypeOfNContiguousBytes(sizeof(T)),othRank(),tag,MPI_COMM_WORLD,&request
    
    MpiRequest request;
    
    if(SR==MpiSendOrRecvFlag::Recv)
      MPI_Irecv(ARGS);
    else
      MPI_Isend(ARGS);
    
    return request;
    
# undef ARGS
    
#else
    return {(void*)p,tag,SR,length};
#endif
  }
  
#define PROVIDE_MPI_ISEND_OR_IRECV(OPT)					\
  template <TriviallyCopyable T>					\
  INLINE_FUNCTION							\
  MpiRequest mpiI ## OPT (const MpiRank& othRank,			\
			  T* p,						\
			  const size_t& length,				\
			  const int& tag)				\
  {									\
    return mpiISendOrRecv(MpiSendOrRecvFlag::OPT,othRank,p,length,tag);	\
  }
  
  PROVIDE_MPI_ISEND_OR_IRECV(Send);
  
  PROVIDE_MPI_ISEND_OR_IRECV(Recv);
  
#undef PROVIDE_MPI_ISEND_OR_IRECV
  
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
  
  /// Wait for communications to finish
  inline void mpiWaitAll(std::vector<MpiRequest>& requests)
  {
#ifdef USE_MPI
    VERBOSITY_LV3_MASTER_PRINTF("Entering MPI comm wait\n");
    
    VERBOSITY_LV3_MASTER_PRINTF("Waiting for %zu requests\n",requests.size());
    
    MPI_Waitall(requests.size(),&requests[0],MPI_STATUS_IGNORE);
#else
    CRASH("Some work is due to wire the bits");
#endif
  }
}

#endif
