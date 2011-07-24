#define APPRETTO 1

#include "appretto.h"

#include "base/debug.c"
#include "base/routines.c"
#include "base/alloc.c"
#include "base/new_types_operations.c"
#include "base/init.c"
#include "base/close.c"
#include "base/communicate.c"

#include "IO/writer.c"
#include "IO/reader.c"
#include "IO/input.c"

#include "operations/su3_paths.c"
#include "operations/gaugeconf.c"
#include "operations/contract.c"
#include "operations/smear.c"
#include "operations/fourier_transform.c"
#include "operations/gauge_fixing.c"

#include "dirac_operators/dirac_operator.c"
#include "inverters/invert.c"
