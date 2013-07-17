#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <mpi.h>
#include <stdio.h>

#include "debug.h"
#include "global_variables.h"
#include "random.h"
#include "vectors.h"
#include "../new_types/new_types_definitions.h"
#include "../geometry/geometry_eo.h"
#include "../geometry/geometry_lx.h"
#include "../geometry/geometry_Wsklx.h"
#ifdef BGQ
 #include "../bgq/geometry_bgq.h"
#endif
#include "../hmc/gauge/tree_level_Symanzik_force.h"
#include "../hmc/gauge/tree_level_Symanzik_action.h"
#include "../routines/ios.h"

void close_nissa()
{
  master_printf("Closing nissa\n");
  
  //unset lx geometry
  if(nissa_lx_geom_inited) unset_lx_geometry();
  if(nissa_Wsklx_order_inited) unset_Wsklx_order();
  
  //unset eo geometry
  if(nissa_eo_geom_inited) unset_eo_geometry();

#ifdef BGQ
  unset_bgq_geometry();
#endif
  
  //stop the random generator
  if(nissa_loc_rnd_gen_inited) stop_loc_rnd_gen();
  
  //print information over the maximum amount of memory used
  master_printf("Maximal memory used during the run: %d bytes (",nissa_max_required_memory);
  if(rank==0) fprintf_friendly_filesize(stdout,nissa_max_required_memory);
  master_printf(") per rank\n\n");
  
  //check wether there are still allocated vectors
  if(main_nissa_vect.next!=NULL && rank==0)
    {
      printf("Warning, there are still allocated vectors:\n");
      print_all_nissa_vect_content();
      printf("For a total of %d bytes\n",compute_nissa_vect_memory_usage());
    }

  tot_nissa_time+=take_time();
  master_printf("Total time: %lg s\n",tot_nissa_time);
#ifdef COMM_BENCH
  master_printf("Total communication time: %lg s\n",tot_nissa_comm_time);
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  master_printf("   Ciao!\n\n");
  MPI_Finalize();
}
