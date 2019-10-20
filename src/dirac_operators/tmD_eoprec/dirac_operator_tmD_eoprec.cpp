#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "threads/threads.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "dirac_operator_tmD_eoprec_portable.cpp"
