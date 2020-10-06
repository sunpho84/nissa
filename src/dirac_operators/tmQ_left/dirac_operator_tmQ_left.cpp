#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/su3_op.hpp"
#include "communicate/borders.hpp"
#include "base/vectors.hpp"
#include "base/thread_macros.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "dirac_operator_tmQ_left_portable.cpp"

