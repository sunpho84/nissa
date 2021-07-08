#ifndef _EXPORT_CONF_TO_EXTERNAL_LIB_HPP
#define _EXPORT_CONF_TO_EXTERNAL_LIB_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_DDALPHAAMG
# include "base/DDalphaAMG_bridge.hpp"
#endif

#ifdef USE_QUDA
# include "base/quda_bridge.hpp"
#endif

#include "io/checksum.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  /// Keep track of the exported conf
  bool export_gauge_conf_to_external_lib(quad_su3 *conf);
}

#endif
