#ifndef _TYPES_H
#define _TYPES_H

#include "../../../../src/new_types/new_types_definitions.h"

typedef struct
{
  double kappa;
  double mass;
  momentum_t bc;
  double zmp;
} quark_info;

typedef struct
{
  double alpha;
  double c1;
  momentum_t bc;
  double zmp;
} gluon_info;

typedef complex corr16[16];

#endif
