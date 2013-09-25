#ifndef _GLOBAL_VARIABLES_H
#define _GLOBAL_VARIABLES_H

#ifdef USE_MPI
 #include <mpi.h>
#endif

#define ONLY_INSTANTIATION
#include "global_variables.cpp"

#endif
