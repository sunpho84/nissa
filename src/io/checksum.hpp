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
  struct checksum
  {
    uint32_t data[2];
    
    CUDA_HOST_AND_DEVICE
    inline uint32_t& operator[](const int i)
    {
      return data[i];
    }
    
    CUDA_HOST_AND_DEVICE
    inline const uint32_t& operator[](const int i) const
    {
      return data[i];
    }
  };
  
  void checksum_compute_ildg_data(uint32_t *check,void *data,size_t bps);
  
#define DO1(buf) crc=crc_table[((int)crc^(*buf++))&0xff]^(crc>>8);
#define DO2(buf) DO1(buf);DO1(buf);
#define DO4(buf) DO2(buf);DO2(buf);
#define DO8(buf) DO4(buf);DO4(buf);
  
  CUDA_HOST_AND_DEVICE
  inline uint32_t ildg_crc32(uint32_t crc,const unsigned char *buf,size_t len)
  {
    static uint32_t crc_table[256]={
      0x00000000L,0x77073096L,0xee0e612cL,0x990951baL,0x076dc419L,
      0x706af48fL,0xe963a535L,0x9e6495a3L,0x0edb8832L,0x79dcb8a4L,
      0xe0d5e91eL,0x97d2d988L,0x09b64c2bL,0x7eb17cbdL,0xe7b82d07L,
      0x90bf1d91L,0x1db71064L,0x6ab020f2L,0xf3b97148L,0x84be41deL,
      0x1adad47dL,0x6ddde4ebL,0xf4d4b551L,0x83d385c7L,0x136c9856L,
      0x646ba8c0L,0xfd62f97aL,0x8a65c9ecL,0x14015c4fL,0x63066cd9L,
      0xfa0f3d63L,0x8d080df5L,0x3b6e20c8L,0x4c69105eL,0xd56041e4L,
      0xa2677172L,0x3c03e4d1L,0x4b04d447L,0xd20d85fdL,0xa50ab56bL,
      0x35b5a8faL,0x42b2986cL,0xdbbbc9d6L,0xacbcf940L,0x32d86ce3L,
      0x45df5c75L,0xdcd60dcfL,0xabd13d59L,0x26d930acL,0x51de003aL,
      0xc8d75180L,0xbfd06116L,0x21b4f4b5L,0x56b3c423L,0xcfba9599L,
      0xb8bda50fL,0x2802b89eL,0x5f058808L,0xc60cd9b2L,0xb10be924L,
      0x2f6f7c87L,0x58684c11L,0xc1611dabL,0xb6662d3dL,0x76dc4190L,
      0x01db7106L,0x98d220bcL,0xefd5102aL,0x71b18589L,0x06b6b51fL,
      0x9fbfe4a5L,0xe8b8d433L,0x7807c9a2L,0x0f00f934L,0x9609a88eL,
      0xe10e9818L,0x7f6a0dbbL,0x086d3d2dL,0x91646c97L,0xe6635c01L,
      0x6b6b51f4L,0x1c6c6162L,0x856530d8L,0xf262004eL,0x6c0695edL,
      0x1b01a57bL,0x8208f4c1L,0xf50fc457L,0x65b0d9c6L,0x12b7e950L,
      0x8bbeb8eaL,0xfcb9887cL,0x62dd1ddfL,0x15da2d49L,0x8cd37cf3L,
      0xfbd44c65L,0x4db26158L,0x3ab551ceL,0xa3bc0074L,0xd4bb30e2L,
      0x4adfa541L,0x3dd895d7L,0xa4d1c46dL,0xd3d6f4fbL,0x4369e96aL,
      0x346ed9fcL,0xad678846L,0xda60b8d0L,0x44042d73L,0x33031de5L,
      0xaa0a4c5fL,0xdd0d7cc9L,0x5005713cL,0x270241aaL,0xbe0b1010L,
      0xc90c2086L,0x5768b525L,0x206f85b3L,0xb966d409L,0xce61e49fL,
      0x5edef90eL,0x29d9c998L,0xb0d09822L,0xc7d7a8b4L,0x59b33d17L,
      0x2eb40d81L,0xb7bd5c3bL,0xc0ba6cadL,0xedb88320L,0x9abfb3b6L,
      0x03b6e20cL,0x74b1d29aL,0xead54739L,0x9dd277afL,0x04db2615L,
      0x73dc1683L,0xe3630b12L,0x94643b84L,0x0d6d6a3eL,0x7a6a5aa8L,
      0xe40ecf0bL,0x9309ff9dL,0x0a00ae27L,0x7d079eb1L,0xf00f9344L,
      0x8708a3d2L,0x1e01f268L,0x6906c2feL,0xf762575dL,0x806567cbL,
      0x196c3671L,0x6e6b06e7L,0xfed41b76L,0x89d32be0L,0x10da7a5aL,
      0x67dd4accL,0xf9b9df6fL,0x8ebeeff9L,0x17b7be43L,0x60b08ed5L,
      0xd6d6a3e8L,0xa1d1937eL,0x38d8c2c4L,0x4fdff252L,0xd1bb67f1L,
      0xa6bc5767L,0x3fb506ddL,0x48b2364bL,0xd80d2bdaL,0xaf0a1b4cL,
      0x36034af6L,0x41047a60L,0xdf60efc3L,0xa867df55L,0x316e8eefL,
      0x4669be79L,0xcb61b38cL,0xbc66831aL,0x256fd2a0L,0x5268e236L,
      0xcc0c7795L,0xbb0b4703L,0x220216b9L,0x5505262fL,0xc5ba3bbeL,
      0xb2bd0b28L,0x2bb45a92L,0x5cb36a04L,0xc2d7ffa7L,0xb5d0cf31L,
      0x2cd99e8bL,0x5bdeae1dL,0x9b64c2b0L,0xec63f226L,0x756aa39cL,
      0x026d930aL,0x9c0906a9L,0xeb0e363fL,0x72076785L,0x05005713L,
      0x95bf4a82L,0xe2b87a14L,0x7bb12baeL,0x0cb61b38L,0x92d28e9bL,
      0xe5d5be0dL,0x7cdcefb7L,0x0bdbdf21L,0x86d3d2d4L,0xf1d4e242L,
      0x68ddb3f8l,0x1fda836eL,0x81be16cdL,0xf6b9265bL,0x6fb077e1L,
      0x18b74777L,0x88085ae6L,0xff0f6a70L,0x66063bcaL,0x11010b5cL,
      0x8f659effL,0xf862ae69L,0x616bffd3L,0x166ccf45L,0xa00ae278L,
      0xd70dd2eeL,0x4e048354L,0x3903b3c2L,0xa7672661L,0xd06016f7L,
      0x4969474dL,0x3e6e77dbL,0xaed16a4aL,0xd9d65adcL,0x40df0b66L,
      0x37d83bf0L,0xa9bcae53L,0xdebb9ec5L,0x47b2cf7fL,0x30b5ffe9L,
      0xbdbdf21cL,0xcabac28aL,0x53b39330L,0x24b4a3a6L,0xbad03605L,
      0xcdd70693L,0x54de5729L,0x23d967bfL,0xb3667a2eL,0xc4614ab8L,
      0x5d681b02L,0x2a6f2b94L,0xb40bbe37L,0xc30c8ea1L,0x5a05df1bL,
      0x2d02ef8dL
    };
    
    if(buf==0) return 0L;
    
    crc^=0xffffffffL;
    
    while(len>=8)
      {
	DO8(buf);
	len-=8;
      }
    
    if(len)
      do {DO1(buf);}
      while(--len);
    
    return crc^0xffffffffL;
  }
  
  namespace checksum_getter
  {
    constexpr size_t max_buf_size=sizeof(quad_su3);
    
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
	unsigned char tmp[checksum_getter::max_buf_size];
	
	switch(prec)
	  {
	  case 64:
	    change_endianness((double*)tmp,(double*)buf,bps/sizeof(double),0);
	    break;
	  case 32:
	    change_endianness((float*)tmp,(float*)buf,bps/sizeof(float),0);
		break;
	  default:
#ifndef COMPILING_FOR_DEVICE
	    crash("unknown precision %d",prec);
#endif
	    break;
	  }
	
	return ildg_crc32(crc,tmp,bps);
      }
    else
      return ildg_crc32(crc,(unsigned char*)buf,bps);
  }
  
  template <typename T,
	    typename F>
  void locReduce(T *loc_res,T *buf,int64_t n,const int nslices,F&& f)
  {
    master_printf("%s function\n",__PRETTY_FUNCTION__);
    
    if(n%nslices)
      crash("number of elements %ld not divisible by number of slices %d",n,nslices);
    
    const int64_t nori_per_slice=n/nslices;
    int64_t nper_slice=n/nslices;
    verbosity_lv2_master_printf("n: %lld, nslices: %d, nori_per_slice: %lld nper_slice: %ld\n",nslices,nori_per_slice,nper_slice);
    
    const double init_time=take_time();
    while(nper_slice>1)
      {
	const int64_t stride=(nper_slice+1)/2;
	const int64_t nreductions_per_slice=nper_slice/2;
	const int64_t nreductions=nreductions_per_slice*nslices;
	verbosity_lv3_master_printf("nper_slice: %lld, stride: %lld, nreductions_per_slice: %lld, nreductions: %lld \n",nper_slice,stride,nreductions_per_slice,nreductions);
	
	
	NISSA_PARALLEL_LOOP(ireduction,0,nreductions)
	  {
	    const int64_t islice=ireduction%nslices;
	    const int64_t ireduction_in_slice=ireduction/nslices;
	    const int64_t first=ireduction_in_slice+nori_per_slice*islice;
	    const int64_t second=first+stride;
	    
	    //buf[first]+=buf[second];
	    f(buf[first],buf[second]);
	  }
	  NISSA_PARALLEL_LOOP_END;
	  THREAD_BARRIER();
	  
	  
	nper_slice=stride;
      }
    
    master_printf("reduction ended, took %lg s\n",take_time()-init_time);
    
    for(int islice=0;islice<nslices;islice++)
      loc_res[islice]=buf[islice*nori_per_slice];
  }
  
  template <typename T>
  void checksum_compute_nissa_data(checksum& check,const T& data,int prec,const size_t bps)
  {
    const double init_time=take_time();
    
    master_printf("   allocating buffer\n");
    checksum* buff=get_reducing_buffer<checksum>(locVol);
    master_printf("   finished allocating the buffer, took %lg s\n",take_time()-init_time);
    
    master_printf("   entering loop\n");
    
    if(bps>checksum_getter::max_buf_size)
      crash("please increase the buf size to hold %zu bps",bps);
    
    const double init_fill_time=take_time();
    
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
    
    master_printf("   finished filling the buffer, took %lg s\n",take_time()-init_fill_time);
    
    master_printf("   starting local reduction\n");
    
    checksum loc_check;
    locReduce(&loc_check,buff,locVol,1,[] CUDA_DEVICE (checksum& res,const checksum& acc)  __attribute__((always_inline))
    {
      for(int i=0;i<2;i++)
	res[i]^=acc[i];
    });
    
    {
    master_printf("   testing global reduction\n");
    const double init_glbred_time=take_time();
    int a[2]{},b[2]{};
    MPI_Allreduce(a,b,2,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    master_printf("   finished test, took %lg s\n",take_time()-init_glbred_time);
    }
    
    master_printf("   starting global reduction\n");
    const double init_glbred_time=take_time();
    MPI_Allreduce(loc_check.data,check.data,2,MPI_UNSIGNED,MPI_BXOR,MPI_COMM_WORLD);
    master_printf("   finished glb reducing buffer, took %lg s\n",take_time()-init_glbred_time);
    
    master_printf("time to compute checksum: %lg (%s) %zu bps\n",take_time()-init_time,__PRETTY_FUNCTION__,bps);
  }
}

#endif
