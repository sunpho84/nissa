#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "new_types/new_types_definitions.h"
#include "new_types/su3.h"
#include "communicate/communicate.h"
#include "base/global_variables.h"
#include "base/vectors.h"
#include "base/thread_macros.h"
#include "linalgs/linalgs.h"
#include "operations/su3_paths/topological_charge.h"
#ifdef USE_THREADS
 #include "routines/thread.h"
#endif

#include "dirac_operator_WclovQ_portable.cpp"

