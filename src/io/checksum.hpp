#ifndef _CHECKSUM_HPP
#define _CHECKSUM_HPP

#include <stdint.h>
#include <string.h>

#include "endianness.hpp"
#include "routines/ios.hpp"
#include "geometry/geometry_eo.hpp"

namespace nissa
{
  typedef uint32_t checksum[2];
  
  uint32_t ildg_crc32(uint32_t crc,const unsigned char *buf,size_t len);
  void checksum_compute_ildg_data(uint32_t *check,void *data,size_t bps);
  
  namespace checksum_getter
  {
    template <typename T>
    auto get_data(const T* data,const int& ivol)
    {
      return data+ivol;
    }
    
    template <typename T>
    auto get_data(const eo_ptr<T>& data,const int& ivol)
    {
      const int eo=loclx_parity[ivol];
      const int loceo=loceo_of_loclx[ivol];
      
      return data[eo]+loceo;
    }
  }

  template <typename T>
  uint32_t ildg_crc32_fix_endianness(uint32_t crc,const T *buf,int prec)
  {
    const size_t len=sizeof(T);
    
    if(little_endian)
      {
	unsigned char temp_buf[len];
	switch(prec)
	  {
	  case 64:
	    change_endianness((double*)temp_buf,(double*)buf,len/sizeof(double),0);
	    break;
	  case 32:
	    change_endianness((float*)temp_buf,(float*)buf,len/sizeof(float),0);
	    break;
	  default:
	    crash("unknown precision %d",prec);
	    break;
	  }
	
	return ildg_crc32(crc,temp_buf,len);
      }
    else
      return ildg_crc32(crc,(unsigned char*)buf,len);
  }
  
  template <typename T>
  void checksum_compute_nissa_data(uint32_t *check,const T& data,int prec)
  {
    const double init_time=take_time();
    
    uint32_t loc_check[2]={0,0};
    
    NISSA_LOC_VOL_LOOP(ivol)
      {
	const coords_t& X=glbCoordOfLoclx[ivol];
	uint32_t ildg_ivol=X[0];
	for(int mu=NDIM-1;mu>0;mu--) ildg_ivol=ildg_ivol*glbSize[mu]+X[mu];
	uint32_t crc_rank[2]={ildg_ivol%29,ildg_ivol%31};
	
	uint32_t temp=ildg_crc32_fix_endianness(0,checksum_getter::get_data(data,ivol),prec);
	
	for(int i=0;i<2;i++) loc_check[i]^=temp<<crc_rank[i]|temp>>(32-crc_rank[i]);
      }
    
    MPI_Allreduce(loc_check,check,2,MPI_UNSIGNED,MPI_BXOR,MPI_COMM_WORLD);
    
    master_printf("time to compute checksum: %lg (%s) %zu bps\n",take_time()-init_time,__PRETTY_FUNCTION__,sizeof(T));
  }
}

#endif
