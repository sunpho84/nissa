#pragma once

void close_appretto()
{
  if(random_initialized)
    {
      close_random();
      random_initialized=0;
    }
  
  if(main_appretto_vect.next!=NULL && rank==0)
    {
      printf("Warning, there are still allocated vectors:\n");
      print_all_appretto_vect_content();
      printf("For a total of %d bytes\n",compute_appretto_vect_memory_usage());
    }
  MPI_Finalize();
}
