#include <mpi.h>
#include <math.h>
#include <string.h>
#include <lemon.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#include "base/global.c"
#include "base/init.c"
#include "base/close.c"
#include "base/communicate.c"

#include "IO/writer.c"
#include "IO/reader.c"
#include "IO/input.c"

#include "operations/contract.c"
#include "operations/smear.c"
#include "operations/fourier_transform.c"
#include "operations/gauge_fixing.c"

#include "types/dirac.c"
#include "types/su3.c"
#include "types/su3_paths.c"
#include "types/gaugeconf.c"

#include "dirac_operators/dirac_operator.c"
#include "inverters/invert.c"
