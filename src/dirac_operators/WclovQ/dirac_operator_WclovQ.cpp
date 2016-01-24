#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/su3_op.hpp"
#include "communicate/borders.hpp"
#include "base/vectors.hpp"
#include "base/thread_macros.hpp"
#include "linalgs/linalgs.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "dirac_operator_WclovQ_portable.cpp"

