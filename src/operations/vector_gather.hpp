#ifndef _VECTOR_GATHER_H
#define _VECTOR_GATHER_H

namespace nissa
{
  void average_list_of_gathered_vector_sites(double *vec,int *sites,int nsites,int dps);
  void gathered_vector_cubic_symmetrize(double *vec,int dps);
  void gathered_vector_mirrorize(double *vec,int dps);
  void vector_gather(char *glb,char *loc,int bps,int dest_rank);
}

#endif
