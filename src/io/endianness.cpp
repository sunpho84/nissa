#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_ENDIANNESS
# include "endianness.hpp"

namespace nissa
{
  //check the endianness of the machine
  void check_endianness()
  {
    little_endian=1;
    little_endian=(int)(*(char*)(&little_endian));
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
}
