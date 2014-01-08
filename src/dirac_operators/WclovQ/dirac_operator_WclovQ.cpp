#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "base/thread_macros.hpp"
#include "linalgs/linalgs.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "dirac_operator_WclovQ_portable.cpp"

