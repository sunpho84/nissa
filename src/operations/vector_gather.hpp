#ifndef _VECTOR_GATHER_HPP
#define _VECTOR_GATHER_HPP

#include <stdlib.h>

namespace nissa
{
  void average_list_of_gathered_vector_sites(double *vec,int *sites,int nsites,int dps);
  void gathered_vector_cubic_symmetrize(double *vec,int dps);
  void gathered_vector_mirrorize(double *vec,int dps);
  void vector_gather(char *glb,char *loc,size_t bps,int dest_rank);
}

#endif
