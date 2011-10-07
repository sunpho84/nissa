#pragma once

void close_appretto()
{
  master_printf("Closing Appretto");
  
  //stop the random generator
  stop_loc_rnd_gen();
  
  //print information over the maximum amount of memory used
  master_printf("Maximal memory used during the run: %d bytes (",appretto_max_required_memory);
  if(rank==0) fprintf_friendly_filesize(stdout,appretto_max_required_memory);
  master_printf(") per process\n\n");
  
  //check wether there are still allocated vectors
  if(main_appretto_vect.next!=NULL && rank==0)
    {
      printf("Warning, there are still allocated vectors:\n");
      print_all_appretto_vect_content();
      printf("For a total of %d bytes\n",compute_appretto_vect_memory_usage());
    }
  MPI_Barrier(MPI_COMM_WORLD);
  
  master_printf("   Ciao!\n\n");
  MPI_Finalize();
}
