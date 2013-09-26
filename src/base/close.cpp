#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>
#include <stdio.h>

#include "debug.hpp"
#include "global_variables.hpp"
#include "random.hpp"
#include "vectors.hpp"
#include "new_types/new_types_definitions.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_Wsklx.hpp"
#ifdef USE_VNODES
 #include "geometry/geometry_vir.hpp"
#endif
#include "hmc/gauge/tree_level_Symanzik_force.hpp"
#include "hmc/gauge/tree_level_Symanzik_action.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  void close_nissa()
  {
    master_printf("Closing nissa\n");
    
    //unset lx geometry
    if(lx_geom_inited) unset_lx_geometry();
    if(Wsklx_order_inited) unset_Wsklx_order();
    
    //unset eo geometry
    if(eo_geom_inited) unset_eo_geometry();
    
#ifdef USE_VNODES
    unset_vir_geometry();
#endif
    
    //stop the random generator
    if(loc_rnd_gen_inited) stop_loc_rnd_gen();
    
    //print information over the maximum amount of memory used
    master_printf("Maximal memory used during the run: %d bytes (",max_required_memory);
    if(rank==0) fprintf_friendly_filesize(stdout,max_required_memory);
    master_printf(") per rank\n\n");
    
    //check wether there are still allocated vectors
    if(main_vect.next!=NULL && rank==0)
      {
	printf("Warning, there are still allocated vectors:\n");
	print_all_vect_content();
	printf("For a total of %d bytes\n",compute_vect_memory_usage());
      }
    
    tot_time+=take_time();
    master_printf("Total time: %lg s\n",tot_time);
#ifdef COMM_BENCH
    master_printf("Total communication time: %lg s\n",tot_comm_time);
#endif
    
    MPI_Barrier(MPI_COMM_WORLD);
    master_printf("   Ciao!\n\n");
    MPI_Finalize();
  }
}
