#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "linalgs/linalgs.hpp"
#include "measures/gauge/topological_charge.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "dirac_operator_WclovQ_portable.cpp"

