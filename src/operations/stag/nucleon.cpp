#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"

#include "stag.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "nucleon.hpp"

namespace nissa
{
  //compute the local nucleon correlator
  void measure_time_nucleon_corr(quad_su3 **conf,theory_pars_t &theory_pars,nucleon_meas_pars_t &nucleon_meas_pars,int iconf,int conf_created)
  {
  }
}
