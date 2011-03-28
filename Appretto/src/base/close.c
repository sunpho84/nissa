#pragma once

#include "random.c"

void close_appretto()
{
  if(random_initialized)
    {
      close_random();
      random_initialized=0;
    }
  MPI_Finalize();
}
