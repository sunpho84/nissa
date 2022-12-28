#ifndef _CHECKSUM_HPP
#define _CHECKSUM_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <stdint.h>
#include <string.h>

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
}

#endif
