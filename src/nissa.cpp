#include "nissa.h"

#include "base/macro.h"
#include "base/new_types.cpp"
#include "base/global.cpp"

#include "headers.h"

#include "base/routines.cpp"
#include "base/vectors.cpp"
#include "base/new_types_operations.cpp"
#include "base/init.cpp"
#include "base/close.cpp"
#include "base/debug.cpp"
#include "base/communicate.cpp"
#include "base/random.cpp"

#include "geometry/geometry_lx.cpp"
#include "geometry/geometry_eo.cpp"
#include "geometry/geometry.h"

#include "linalgs/linalgs.cpp"

#include "IO/checksum.cpp"
#include "IO/writer.cpp"
#include "IO/reader.cpp"
#include "IO/input.cpp"

#include "operations/su3_paths.cpp"
#include "operations/fft.cpp"
#include "operations/gaugeconf.cpp"
#include "operations/contract.cpp"
#include "operations/smear.cpp"
#include "operations/fourier_transform.cpp"
#include "operations/gauge_fixing.cpp"
#include "operations/remap_vector.cpp"

#include "dirac_operators/dirac_operators.cpp"
#include "inverters/inverters.cpp"

#include "hmc/hmc.cpp"
