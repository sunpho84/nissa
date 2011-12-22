#pragma once

void close_nissa()
{
  master_printf("Closing nissa\n");
  
  //unset lx geometry
  if(nissa_lx_geom_inited) unset_lx_geometry();
  
  //unset eo geometry
  if(nissa_eo_geom_inited) unset_eo_geometry();
  
  //stop the random generator
  if(nissa_loc_rnd_gen_inited) stop_loc_rnd_gen();
  
  //print information over the maximum amount of memory used
  master_printf("Maximal memory used during the run: %d bytes (",nissa_max_required_memory);
  if(rank==0) fprintf_friendly_filesize(stdout,nissa_max_required_memory);
  master_printf(") per process\n\n");
  
  //check wether there are still allocated vectors
  if(main_nissa_vect.next!=NULL && rank==0)
    {
      printf("Warning, there are still allocated vectors:\n");
      print_all_nissa_vect_content();
      printf("For a total of %d bytes\n",compute_nissa_vect_memory_usage());
    }
  MPI_Barrier(MPI_COMM_WORLD);
  
  master_printf("   Ciao!\n\n");
  MPI_Finalize();
}
