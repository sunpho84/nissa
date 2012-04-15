#include "nissa.h"

#include "base/macro.h"
#include "base/new_types.c"
#include "base/global.c"

#include "headers.h"

#include "base/routines.c"
#include "base/vectors.c"
#include "base/new_types_operations.c"
#include "base/init.c"
#include "base/close.c"
#include "base/debug.c"
#include "base/communicate.c"
#include "base/random.c"

#include "linalgs/linalgs.c"

#include "IO/checksum.c"
#include "IO/writer.c"
#include "IO/reader.c"
#include "IO/input.c"

#include "operations/su3_paths.c"
#include "operations/fft.c"
#include "operations/gaugeconf.c"
#include "operations/contract.c"
#include "operations/smear.c"
#include "operations/fourier_transform.c"
#include "operations/gauge_fixing.c"
#include "operations/remap_vector.c"

#include "dirac_operators/dirac_operators.c"
#include "inverters/inverters.c"

#include "hmc/hmc.c"
