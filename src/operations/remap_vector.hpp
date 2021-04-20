#ifndef _REMAP_VECTOR_HPP
#define _REMAP_VECTOR_HPP

#ifndef EXTERN_REMAP
# define EXTERN_REMAP extern
#endif

#include <communicate/all_to_all.hpp>

namespace nissa
{
  struct vector_remap_t : all_to_all_comm_t
  {
    vector_remap_t(int nel_fr,void (*index)(Rank &irank_to,LocLxSite &iel_to,const LocLxSite& iel_fr,void *pars),void *pars);
    
    void remap(void *out,void *in,size_t bps)
    {
      communicate(out,in,bps);
    }
  };
  
  //local direction geometry
  EXTERN_REMAP vector_remap_t *remap_lx_to_locd[NDIM];
  EXTERN_REMAP vector_remap_t *remap_locd_to_lx[NDIM];
  EXTERN_REMAP Coords<int64_t> max_locd_perp_size_per_dir,locd_perp_size_per_dir;
  EXTERN_REMAP int64_t max_locd_size;
  EXTERN_REMAP Coords<int64_t> locd_size_per_dir;
  
  void remap_lx_vector_to_locd(void *out,void *in,int nbytes,const Dir& mu);
  void remap_locd_vector_to_lx(void *out,void *in,int nbytes,const Dir& mu);
}

#endif
