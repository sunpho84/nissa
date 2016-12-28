#ifndef _REMAP_VECTOR_HPP
#define _REMAP_VECTOR_HPP

#include "communicate/all_to_all.hpp"

namespace nissa
{
  struct vector_remap_t : all_to_all_comm_t
  {
    vector_remap_t(int nel_out,void (*index)(int &irank_to,int &iel_to,int iel_fr,void *pars),void *pars);
    void remap(void *out,void *in,size_t bps){communicate(out,in,bps);}
  };
  
  //local direction geometry
  extern vector_remap_t *remap_lx_to_locd[NDIM];
  extern vector_remap_t *remap_locd_to_lx[NDIM];
  extern int max_locd_perp_size_per_dir[NDIM],locd_perp_size_per_dir[NDIM];
  extern int max_locd_size;
  
  void remap_lx_vector_to_locd(void *out,void *in,int nbytes,int mu);
  void remap_locd_vector_to_lx(void *out,void *in,int nbytes,int mu);
}

#endif
