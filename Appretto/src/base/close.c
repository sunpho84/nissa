#pragma once

void close_appretto()
{
  //close the random generator
  if(random_initialized)
    {
      close_random();
      random_initialized=0;
    }
  
  //print information over the maximum amount of memory used
  if(rank==0)
    {
      printf("\nMaximal memory used during the run: %d bytes (",appretto_max_required_memory);
      fprintf_friendly_filesize(stdout,appretto_max_required_memory);
      printf(") per process\n\n");
    }
  
  //check wether there are still allocated vectors
  if(main_appretto_vect.next!=NULL && rank==0)
    {
      printf("Warning, there are still allocated vectors:\n");
      print_all_appretto_vect_content();
      printf("For a total of %d bytes\n",compute_appretto_vect_memory_usage());
    }
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Finalize();
}
