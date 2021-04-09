#ifndef _SQUARED_STAPLE_HPP
#define _SQUARED_STAPLE_HPP

#include <mpi.h>

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  typedef su3 squared_staples_t[NDIM][NDIM*(NDIM+1)/2];
  
  void compute_squared_staples_lx_conf(squared_staples_t *out,quad_su3 *conf);
  void compute_summed_squared_staples_eo_conf(eo_ptr<quad_su3> F,eo_ptr<quad_su3> eo_conf);
  void compute_summed_squared_staples_lx_conf(quad_su3 *out,quad_su3 *conf);
  CUDA_HOST_DEVICE void compute_point_summed_squared_staples_eo_conf(quad_su3 staple,eo_ptr<quad_su3> eo_conf,const LocLxSite& A);
  void compute_point_summed_squared_staples_lx_conf(quad_su3 staple,quad_su3 *lx_conf,const LocLxSite& A);
  CUDA_HOST_DEVICE void compute_point_summed_squared_staples_eo_conf_single_dir(su3 staple,eo_ptr<quad_su3> eo_conf,const LocLxSite& A,int mu);
  void compute_point_summed_squared_staples_lx_conf_single_dir(su3 staple,quad_su3 *lx_conf,const LocLxSite& A,int mu);
  void squared_staples_lx_conf_allocate_buffers(eo_ptr<quad_su3> send_buf,eo_ptr<quad_su3> recv_buf);
  void squared_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(quad_su3 *send_buf,quad_su3 *recv_buf,squared_staples_t *out,quad_su3 *conf,int (*nrequest),MPI_Request *request,int thread_id);
  void squared_staples_lx_conf_compute_fw_surf_fw_staples(squared_staples_t *out,quad_su3 *conf,int thread_id);
  void squared_staples_lx_conf_compute_non_fw_surf_bw_staples(squared_staples_t *out,quad_su3 *conf,int thread_id);
  void squared_staples_lx_conf_compute_non_fw_surf_fw_staples(squared_staples_t *out,quad_su3 *conf,int thread_id);
  void squared_staples_lx_conf_finish_communicating_fw_surf_bw_staples(squared_staples_t *out,quad_su3 *recv_buf,int (*nrequest),MPI_Request *request,int thread_id);
  void squared_staples_lx_conf_finish_communicating_lower_surface(quad_su3 *conf,quad_su3 *recv_buf,int (*nrequest),MPI_Request *request,int thread_id);
  void squared_staples_lx_conf_start_communicating_lower_surface(quad_su3 *send_buf,quad_su3 *recv_buf,quad_su3 *conf,int (*nrequest),MPI_Request *request,int thread_id);
}

#endif
