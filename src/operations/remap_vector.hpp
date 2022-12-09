#ifndef _REMAP_VECTOR_HPP
#define _REMAP_VECTOR_HPP

#ifndef EXTERN_REMAP
# define EXTERN_REMAP extern
#endif

#include <base/debug.hpp>
#include <communicate/all_to_all.hpp>

namespace nissa
{
  struct vector_remap_t :
    all_to_all_comm_t
  {
    /// Initializes the remap, with a lambda
    template <typename F>
    vector_remap_t(const int& nel_fr,
		   F&& index)
    {
      all_to_all_scattering_list_t sl;
      for(int iel_fr=0;iel_fr<nel_fr;iel_fr++)
	{
	  const auto [rank_to,iel_to]=index(iel_fr);
	  
	  if(rank_to>=nranks or rank_to<0)
	    crash("destination rank %d does not exist!",rank_to);
	  
	  sl.push_back(std::make_pair(iel_fr,iel_to*nranks+rank_to));
	}
      setup_knowing_where_to_send(sl);
    }
    
    // /// Initialize the remap, with an old-style lambda
    // vector_remap_t(int nel_fr,
    // 		   void (*index)(int& irank_to,
    // 				 int& iel_to,
    // 				 int iel_fr,
    // 				 void* pars),
    // 		   void* pars) :
    //   vector_remap_t(nel_fr,
    // 		     [pars,index](const int& iel_fr)
    // 		     {
    // 		       int rank_to,iel_to;
    // 		       index(rank_to,iel_to,iel_fr,pars);
		       
    // 		       return std::make_pair(rank_to,iel_to);
    // 		     })
    // {
    // }
    
    void remap(void *out,
	       void *in,
	       size_t bps) const
    {
      communicate(out,in,bps);
    }
    
    template <typename T>
    void remap(T* out,
	       const T* in) const
    {
      remap(out,in,sizeof(T));
    }
  };
  
  //local direction geometry
  EXTERN_REMAP vector_remap_t *remap_lx_to_locd[NDIM];
  EXTERN_REMAP vector_remap_t *remap_locd_to_lx[NDIM];
  EXTERN_REMAP int max_locd_perp_size_per_dir[NDIM],locd_perp_size_per_dir[NDIM];
  EXTERN_REMAP int max_locd_size,locd_size_per_dir[NDIM];
  
  void remap_lx_vector_to_locd(void *out,void *in,int nbytes,int mu);
  void remap_locd_vector_to_lx(void *out,void *in,int nbytes,int mu);
}

#endif
