#ifndef _CHECKSUM_HPP
#define _CHECKSUM_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdint.h>
#include <string.h>

#include "endianness.hpp"
#include "routines/ios.hpp"
#include "geometry/geometry_eo.hpp"



#include "linalgs/reduce.hpp"

namespace nissa
{
  typedef uint32_t checksum[2];
  
  CUDA_HOST_AND_DEVICE
  uint32_t ildg_crc32(uint32_t crc,const unsigned char *buf,size_t len);
  
  void checksum_compute_ildg_data(uint32_t *check,void *data,size_t bps);
  
  namespace checksum_getter
  {
    template <typename T>
    CUDA_HOST_AND_DEVICE
    auto get_data(const T* data,const int& ivol,const size_t bps=sizeof(T))
    {
      return (const char*)data+bps*ivol;
    }
    
    template <typename T>
    CUDA_HOST_AND_DEVICE
    auto get_data(const eo_ptr<T>& data,const int& ivol,const size_t bps=sizeof(T))
    {
      const int eo=loclx_parity[ivol];
      const int loceo=loceo_of_loclx[ivol];
      
      return (const char*)data[eo]+bps*loceo;
    }
  }
  
  template <typename T>
  CUDA_HOST_AND_DEVICE
  uint32_t ildg_crc32_fix_endianness(uint32_t crc,const T *buf,int prec,const size_t bps=sizeof(T))
  {
    if(little_endian)
      {
	uint32_t res=0;
	
	// Swap endianness, compute and put back endianness
	for(int iter=0;iter<2;iter++)
	  {
	    switch(prec)
	      {
	      case 64:
		change_endianness((double*)buf,(double*)buf,bps/sizeof(double),0);
		break;
	      case 32:
		change_endianness((float*)buf,(float*)buf,bps/sizeof(float),0);
		break;
	      default:
#ifndef COMPILING_FOR_DEVICE
		crash("unknown precision %d",prec);
#endif
		break;
	      }
	    
	    if(iter==0)
	      res=ildg_crc32(crc,(unsigned char*)buf,bps);
	  }
	
	return res;
      }
    else
      return ildg_crc32(crc,(unsigned char*)buf,bps);
  }
  
  template <typename T,
	    typename F>
  void locReduce(T *loc_res,T *buf,int64_t n,const int nslices,F f)
  {
    if(n%nslices)
      crash("number of elements %ld not divisible by number of slices %d",n,nslices);
    
    const int64_t nori_per_slice=n/nslices;
    int64_t nper_slice=n/nslices;
    verbosity_lv2_master_printf("n: %lld, nslices: %d, nori_per_slice: %lld nper_slice: %ld\n",nslices,nori_per_slice,nper_slice);
    
    while(nper_slice>1)
      {
	const int64_t stride=(nper_slice+1)/2;
	const int64_t nreductions_per_slice=nper_slice/2;
	const int64_t nreductions=nreductions_per_slice*nslices;
	verbosity_lv3_master_printf("nper_slice: %lld, stride: %lld, nreductions_per_slice: %lld, nreductions: %lld\n",nper_slice,stride,nreductions_per_slice,nreductions);
	
	NISSA_PARALLEL_LOOP(ireduction,0,nreductions)
	  {
	    const int64_t islice=ireduction%nslices;
	    const int64_t ireduction_in_slice=ireduction/nslices;
	    const int64_t first=ireduction_in_slice+nori_per_slice*islice;
	    const int64_t second=first+stride;
	    
	    f(buf[first],buf[second]);
	  }
	NISSA_PARALLEL_LOOP_END;
	THREAD_BARRIER();
	
	nper_slice=stride;
      }
    
    for(int islice=0;islice<nslices;islice++)
      memcpy(loc_res[islice],buf[islice*nori_per_slice],sizeof(T));
  }
  
  CUDA_HOST_AND_DEVICE
  inline void checksumReducer(checksum& res,const checksum& acc)
    {
      for(int i=0;i<2;i++)
	res[i]^=acc[i];
    }
  
  template <typename T>
  void checksum_compute_nissa_data(checksum& check,const T& data,int prec,const size_t bps=sizeof(T))
  {
    const double init_time=take_time();
    
    master_printf("   allocating buffer\n");
    checksum* buff=get_reducing_buffer<checksum>(locVol);
    
    // uint32_t loc_check[2]={0,0};
    
    master_printf("   entering loop\n");
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	const coords_t& X=glbCoordOfLoclx[ivol];
	uint32_t ildg_ivol=X[0];
	for(int mu=NDIM-1;mu>0;mu--) ildg_ivol=ildg_ivol*glbSize[mu]+X[mu];
	uint32_t crc_rank[2]={ildg_ivol%29,ildg_ivol%31};
	
	uint32_t temp=ildg_crc32_fix_endianness(0,checksum_getter::get_data(data,ivol,bps),prec,bps);
	
	for(int i=0;i<2;i++) buff[ivol][i]=temp<<crc_rank[i]|temp>>(32-crc_rank[i]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    master_printf("   starting local reducion\n");
    
    checksum loc_check;
    locReduce(&loc_check,buff,locVol,1,checksumReducer);
    
    master_printf("   starting global reductiond\n");
    MPI_Allreduce(loc_check,check,2,MPI_UNSIGNED,MPI_BXOR,MPI_COMM_WORLD);
    
    master_printf("time to compute checksum: %lg (%s) %zu bps\n",take_time()-init_time,__PRETTY_FUNCTION__,bps);
  }
}

#endif
