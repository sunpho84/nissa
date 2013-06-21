#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/global_variables.h"
#include "../../base/debug.h"
#include "../../base/thread_macros.h"
#include "../../base/vectors.h"
#include "../../communicate/communicate.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/complex.h"
#include "../../new_types/su3.h"
#ifdef USE_THREADS
 #include "../../routines/thread.h"
#endif

#include "dirac_operator_tmDeoimpr_portable.cpp"
