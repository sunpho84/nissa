#ifndef _EDGES_HPP
#define _EDGES_HPP

#include "communicate.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"

namespace nissa
{
  void communicate_eo_edges(eo_ptr<void> data,comm_t &comm,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site);
  void communicate_eo_quad_su3_edges(eo_ptr<quad_su3> conf);
  void communicate_lx_edges(char *data,comm_t &comm,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site);
  void communicate_lx_quad_su3_edges(quad_su3 *conf);
  void communicate_lx_as2t_su3_edges(as2t_su3 *a);
  void communicate_lx_su3_edges(su3 *u);
}

#endif
