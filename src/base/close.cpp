#include <mpi.h>
#include <stdio.h>

#include "global_variables.h"
#include "random.h"
#include "routines.h"
#include "vectors.h"
#include "../new_types/new_types_definitions.h"
#include "../geometry/geometry_eo.h"
#include "../geometry/geometry_lx.h"
#include "../hmc/gauge/tree_level_Symanzik_force.h"
#include "../hmc/gauge/tree_level_Symanzik_action.h"

void close_nissa()
{
  master_printf("Closing nissa\n");
  
  //unset lx geometry
  if(nissa_lx_geom_inited) unset_lx_geometry();
  
  //unset eo geometry
  if(nissa_eo_geom_inited) unset_eo_geometry();
  
  //stop the random generator
  if(nissa_loc_rnd_gen_inited) stop_loc_rnd_gen();
  
  //stop the tree level Symanzik action calculation structures
  stop_Symanzik_action();
  stop_Symanzik_staples();
  
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
  master_printf("Total communication time: %lg s\n",tot_nissa_comm_time);

  MPI_Barrier(MPI_COMM_WORLD);
  master_printf("   Ciao!\n\n");
  MPI_Finalize();
}
