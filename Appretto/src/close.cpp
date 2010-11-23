#pragma once

#include "random.cpp"

void close_appretto()
{
  if(random_initialized)
    {
      close_random();
      random_initialized=false;
    }
  MPI_Finalize();
}
