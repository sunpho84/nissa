#ifndef _FNV1A64_HPP
#define _FNV1A64_HPP

#include <cstdint>
#include <string_view>

namespace nissa
{
  /// FNV-1a 64-bit offset basis; arithmetic wraps modulo 2^64
  inline constexpr uint64_t fnv1a64Offset=14695981039346656037ULL;
  
  /// Extend an FNV-1a 64-bit hash with the exact bytes of a string view
  constexpr uint64_t fnv1a64Update(uint64_t hash,
                                   const std::string_view bytes)
  {
    for(const unsigned char byte : bytes)
      {
        hash^=byte;
        hash*=1099511628211ULL;
      }
    
    return hash;
  }
  
  /// Hash a complete byte string with FNV-1a 64-bit
  constexpr uint64_t fnv1a64(const std::string_view bytes)
  {
    return fnv1a64Update(fnv1a64Offset,bytes);
  }
  
  static_assert(fnv1a64("")==14695981039346656037ULL);
  static_assert(fnv1a64("a")==12638187200555641996ULL);
}

#endif
