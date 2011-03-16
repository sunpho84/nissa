#include <mpi.h>
#include <math.h>
#include <lemon.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#include "global.c"
#include "init.c"
#include "close.c"
#include "writer.c"
#include "reader.c"
#include "input.c"
#include "dirac.c"
#include "contract.c"
#include "communicate.c"
#include "su3.c"
#include "smear.c"
#include "su3_paths.c"
#include "gaugeconf.c"
#include "dirac_operator.c"
#include "invert.c"
#include "fourier_transform.c"
