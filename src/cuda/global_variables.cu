#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "macros.hpp"

#ifdef ONLY_INSTANTIATION
 #define EXTERN extern
#else
 #define EXTERN
#endif

namespace cuda
{
  EXTERN DEVICE_CONSTANT int loc_vol,loc_volh;
}
