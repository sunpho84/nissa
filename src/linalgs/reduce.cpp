#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_REDUCE
 #include <linalgs/reduce.hpp>

namespace nissa
{
  
#define DEFINE_LOC_REDUCE_OF(TYPE)					\
  /*! Reduce a vector of */						\
  THREADABLE_FUNCTION_4ARG(loc_reduce, TYPE*,loc_res, TYPE*,buf, int64_t,n, int,nslices) \
  {									\
    GET_THREAD_ID();							\
									\
    if(n%nslices)							\
      crash("number of elements %ld not divisible by number of slices %d",n,nslices); \
									\
    const int64_t nori_per_slice=n/nslices;				\
    int64_t nper_slice=n/nslices;					\
									\
    while(nper_slice>1)							\
      {									\
	const int64_t stride=(nper_slice+1)/2;				\
	const int64_t nreductions_per_slice=nper_slice/2;		\
	const int64_t nreductions=nreductions_per_slice*nslices;	\
									\
	NISSA_PARALLEL_LOOP(ireduction,0,nreductions)			\
	  {								\
	    const int64_t islice=ireduction%nslices;			\
	    const int64_t ireduction_in_slice=ireduction/nslices;	\
	    const int64_t first=ireduction_in_slice+nori_per_slice*islice; \
	    const int64_t second=first+stride;				\
									\
	    TYPE ## _summassign(buf[first],buf[second]);		\
	  }								\
	NISSA_PARALLEL_LOOP_END;					\
	THREAD_BARRIER();						\
									\
	nper_slice=stride;						\
      }									\
									\
    for(int islice=0;islice<nslices;islice++)				\
      NAME2(TYPE,copy)(loc_res[islice],buf[islice*nori_per_slice]);	\
  }									\
  THREADABLE_FUNCTION_END
  
  namespace
  {
    /// Implement oth=first
    inline void int64_t_copy(int64_t& oth,const int64_t& first)
    {
      oth=first;
    }
    
    /// Implement oth+=first
    inline void int64_t_summassign(int64_t& oth,const int64_t& first)
    {
      oth+=first;
    }
    
    /// Implement oth=first
    inline void double_copy(double& oth,const double& first)
    {
      oth=first;
    }
    
    /// Implements oth+=first
    inline void double_summassign(double& oth,const double& first)
    {
      oth+=first;
    }
  }
  
  DEFINE_LOC_REDUCE_OF(int64_t)
  DEFINE_LOC_REDUCE_OF(double)
  DEFINE_LOC_REDUCE_OF(complex)
  DEFINE_LOC_REDUCE_OF(float_128)
  DEFINE_LOC_REDUCE_OF(complex_128)
  
  /// Releases the reduction buffer if allocated
  void deallocate_reduction_buffer()
  {
    if(reducing_buffer!=nullptr)
      {
	master_printf("Freeing reduction buffer, used size: %ld\n",reducing_buffer_size);
	reducing_buffer_size=0;
	
	nissa_free(reducing_buffer);
      }
  }
}
