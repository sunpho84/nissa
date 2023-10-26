#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#define EXTERN_REMAP
#include "remap_vector.hpp"

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/all_to_all.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "threads/threads.hpp"

#include "remap_vector.hpp"


namespace nissa
{
  //compute the remapping index to make dir mu local
  auto get_index_make_loc_dir(const int& mu,
			      const int& prp_max_vol)
  {
    return [mu,prp_max_vol](int iloc_lx)
    {
      int glb_perp_site=0;
      for(int nu=0;nu<NDIM;nu++)
	if(mu!=nu)
	  glb_perp_site=glb_perp_site*glbSizes[nu]+glbCoordOfLoclx[iloc_lx][nu];
      
      const int irank_locld=glb_perp_site/prp_max_vol;
      int iloc_locld=glb_perp_site-irank_locld*prp_max_vol;
      iloc_locld=iloc_locld*glbSizes[mu]+glbCoordOfLoclx[iloc_lx][mu];
      
      return std::make_pair(irank_locld,iloc_locld);
    };
  }
  
  //unmake
  auto get_index_unmake_loc_dir(const int mu,
				const int& prp_max_vol)
  {
    return [mu,prp_max_vol](int iloc_locld) // don't make constant
    {
      coords_t c;
      c[mu]=iloc_locld%glbSizes[mu];
      iloc_locld/=glbSizes[mu];
      
      int glb_perp_site=
	iloc_locld+rank*prp_max_vol;
      
      for(int nu=NDIM-1;nu>=0;nu--)
	if(mu!=nu)
	  {
	    c[nu]=glb_perp_site%glbSizes[nu];
	    glb_perp_site/=glbSizes[nu];
	  }
      //int &irank_lx,int &iloc_lx;
      return get_loclx_and_rank_of_coord(c);
    };
  }
  
  //remap to locd
  void remap_lx_vector_to_locd(void *out,void *in,int nbytes,int mu)
  {
    if(remap_lx_to_locd[mu]==NULL)
      remap_lx_to_locd[mu]=
	new vector_remap_t(locVol,get_index_make_loc_dir(mu,max_locd_perp_size_per_dir[mu]));
    remap_lx_to_locd[mu]->remap(out,in,nbytes);
  }
  
  void remap_locd_vector_to_lx(void *out,void *in,int nbytes,int mu)
  {
    if(remap_locd_to_lx[mu]==NULL)
      remap_locd_to_lx[mu]=new vector_remap_t(locd_size_per_dir[mu],get_index_unmake_loc_dir(mu,max_locd_perp_size_per_dir[mu]));
    remap_locd_to_lx[mu]->remap(out,in,nbytes);
  }
}
