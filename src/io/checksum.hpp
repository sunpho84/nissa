#ifndef _CHECKSUM_HPP
#define _CHECKSUM_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <stdint.h>

#include <base/field.hpp>
#include <io/endianness.hpp>
#include <linalgs/reduce.hpp>
#include <routines/ios.hpp>
#include <geometry/geometry_eo.hpp>

namespace nissa
{
  using Checksum=std::array<uint32_t,2>;
  
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  uint32_t crcValue(uint32_t c)
  {
    constexpr uint32_t poly=0xedb88320L;
    
    for(uint32_t j=0;j<8;j++)
      c=(c&1)?poly^(c>>1):(c>>1);
    
    return c;
  }
  
  /// Pairwise reduce the checksum
  struct ChecksumReduceFunctor
  {
    /// Called on value by value basis
    CUDA_DEVICE INLINE_FUNCTION
    void operator()(uint32_t& res,
		    const uint32_t& acc) const
    {
	res^=acc;
    }
  };
  
  /// Computes the ildg checksum
  template <typename T,
	    SpaceTimeLayout STL>
  Checksum ildgChecksum(const LxField<T,STL>& field)
  {
    using ChecksumCArr=
      uint32_t[2];
    
    LxField<ChecksumCArr> buff("buff");
    
    PAR(0,locVol,
	CAPTURE(TO_WRITE(buff),
		TO_READ(field)),
	ivol,
	{
	  const Coords X=
	    glbCoordOfLoclx[ivol];
	  
	  uint32_t ildg_ivol=X[0];
	  for(int mu=NDIM-1;mu>0;mu--)
	    ildg_ivol=ildg_ivol*glbSize[mu]+X[mu];
	  
	  const uint32_t crc_rank[2]=
	    {ildg_ivol%29,ildg_ivol%31};
	  
	  uint32_t crc=0xffffffffL;
	  for(int iDeg=0;iDeg<(LxField<T>::nInternalDegs);iDeg++)
	    {
	      auto temp=
		field(ivol,iDeg);
	      
	      using Fund=
		decltype(temp);
	      
	      EndiannessMask<BigEndian,nativeEndianness,Fund> mask(temp);
	      
	      for(size_t i=0;i<sizeof(temp);i++)
		crc=crcValue(((int)crc^mask[i])&0xff)^(crc>>8);
	    }
	  crc^=0xffffffffL;
	  
	  for(int i=0;i<2;i++)
	    buff[ivol][i]=
	      crc<<crc_rank[i]|crc>>(32-crc_rank[i]);
	});
    
    ChecksumCArr check;
    buff.reduce(check,ChecksumReduceFunctor(),MPI_BXOR);
    
    return {check[0],check[1]};
  }
}

#endif
