#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#define EXTERN_REMAP
#include "remap_vector.hpp"

#include "base/debug.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/all_to_all.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"

#include "remap_vector.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //constructor
  vector_remap_t::vector_remap_t(int nel_fr,void (*index)(int &irank_to,int &iel_to,int iel_fr,void *pars),void *pars)
  {
    all_to_all_scattering_list_t sl;
    for(int iel_fr=0;iel_fr<nel_fr;iel_fr++)
      {
	int rank_to,iel_to;
	index(rank_to,iel_to,iel_fr,pars);
	if(rank_to>=nranks or rank_to<0) crash("destination rank %d does not exist!",rank_to);
	sl.push_back(std::make_pair(iel_fr,iel_to*nranks+rank_to));
      }
    setup_knowing_where_to_send(sl);
  }
  
  //compute the remapping index to make dir mu local
  void index_make_loc_dir(int &irank_locld,int &iloc_locld,int iloc_lx,void *pars)
  {
    int mu=((int*)pars)[0],prp_max_vol=((int*)pars)[1];
    int glb_perp_site=0;
    for(int nu=0;nu<NDIM;nu++) if(mu!=nu) glb_perp_site=glb_perp_site*glb_size[nu]+glb_coord_of_loclx[iloc_lx][nu];
    irank_locld=glb_perp_site/prp_max_vol;
    iloc_locld=glb_perp_site-irank_locld*prp_max_vol;
    iloc_locld=iloc_locld*glb_size[mu]+glb_coord_of_loclx[iloc_lx][mu];
  }
  
  //unmake
  void index_unmake_loc_dir(int &irank_lx,int &iloc_lx,int iloc_locld,void *pars)
  {
    int mu=((int*)pars)[0],prp_max_vol=((int*)pars)[1];
    coords c;
    c[mu]=iloc_locld%glb_size[mu];
    iloc_locld/=glb_size[mu];
    int glb_perp_site=iloc_locld+rank*prp_max_vol;
    for(int nu=NDIM-1;nu>=0;nu--)
      if(mu!=nu)
	{
	  c[nu]=glb_perp_site%glb_size[nu];
	  glb_perp_site/=glb_size[nu];
	}
    get_loclx_and_rank_of_coord(&iloc_lx,&irank_lx,c);
  }
  
  //remap to locd
  void remap_lx_vector_to_locd(void *out,void *in,int nbytes,int mu)
  {
    int pars[2]={mu,max_locd_perp_size_per_dir[mu]};
    if(remap_lx_to_locd[mu]==NULL) remap_lx_to_locd[mu]=new vector_remap_t(loc_vol,index_make_loc_dir,pars);
    remap_lx_to_locd[mu]->remap(out,in,nbytes);
  }
  void remap_locd_vector_to_lx(void *out,void *in,int nbytes,int mu)
  {
    int pars[2]={mu,max_locd_perp_size_per_dir[mu]};
    if(remap_locd_to_lx[mu]==NULL) remap_locd_to_lx[mu]=new vector_remap_t(locd_size_per_dir[mu],index_unmake_loc_dir,pars);
    remap_locd_to_lx[mu]->remap(out,in,nbytes);
  }
}
