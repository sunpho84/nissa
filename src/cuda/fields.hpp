#ifndef _FIELDS_H
#define _FIELDS_H

#include "macros.hpp"

namespace cuda
{
  GLOBAL void field_reset(float *d,const int nper_site);
}

#endif
