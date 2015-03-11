#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/all_to_all.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"

#include "remap_vector.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //constructor
  vector_remap_t::vector_remap_t(int nel_out,void (*index)(int &irank_to,int &iel_to,int iel_fr,void *pars),void *pars)
  {
    all_to_all_scattering_list_t sl;
    for(int iel_out=0;iel_out<nel_out;iel_out++)
      {
	int rank_to,iel_to;
	index(rank_to,iel_to,iel_out,pars);
	if(rank_to>=nranks||rank_to<0) crash("destination rank %d does not exist!",rank_to);
	sl.push_back(std::make_pair(iel_out,iel_to*nranks+rank_to));
      }
    setup_knowing_where_to_send(sl);
  }
}
