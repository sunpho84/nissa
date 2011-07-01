#include <mpi.h>

extern int rank;
extern int loc_vol,loc_bord;

#include "base/types.cu"
#include "base/global.cu"
#include "base/compression.cu"

#include "init/find_gpu.cu"
#include "init/choose_gpu.cu"

#include "dirac_operators/dirac_operator.cu"
