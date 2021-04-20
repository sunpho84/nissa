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

namespace nissa
{
  //constructor
  vector_remap_t::vector_remap_t(int nel_fr,void (*index)(Rank &irank_to,LocLxSite &iel_to,const LocLxSite& iel_fr,void *pars),void *pars)
  {
    all_to_all_scattering_list_t sl;
    for(int iel_fr=0;iel_fr<nel_fr;iel_fr++)
      {
	Rank rank_to;
	LocLxSite iel_to;
	
	index(rank_to,iel_to,iel_fr,pars);
	if(rank_to>=nranks or rank_to<0)
	  crash("destination rank %d does not exist!",rank_to);
	
	sl.push_back(std::make_pair(iel_fr,iel_to()*nranks+rank_to()));
      }
    setup_knowing_where_to_send(sl);
  }
  
  //compute the remapping index to make dir mu local
  void index_make_loc_dir(Rank &irank_locld,LocLxSite &iloc_locld,const LocLxSite& iloc_lx,void *pars)
  {
    const Dir mu=((int*)pars)[0];
    const int prp_max_vol=((int*)pars)[1];
    GlbLxSite glb_perp_site=0;
    
    FOR_ALL_DIRS(nu)
      if(mu!=nu)
	glb_perp_site=glb_perp_site*glbSize(nu)+glbCoordOfLoclx(iloc_lx,nu);
    
    irank_locld=glb_perp_site()/prp_max_vol;
    iloc_locld=glb_perp_site()-irank_locld()*prp_max_vol;
    iloc_locld=iloc_locld*glbSize(mu)()+glbCoordOfLoclx(iloc_lx,mu)();
  }
  
  //unmake
  void index_unmake_loc_dir(Rank &irank_lx,LocLxSite &iloc_lx,const LocLxSite& ext_iloc_locld,void *pars)
  {
    LocLxSite iloc_locld=ext_iloc_locld;
    
    const Dir mu=((int*)pars)[0];
    const int prp_max_vol=((int*)pars)[1];
    
    GlbCoords c;
    c(mu)=iloc_locld()%glbSize(mu);
    
    iloc_locld/=glbSize(mu)();
    
    GlbLxSite glb_perp_site=iloc_locld()+rank*prp_max_vol;
    for(Dir nu=NDIM-1;nu>=0;nu--)
      if(mu!=nu)
	{
	  c(nu)=glb_perp_site%glbSize(nu);
	  glb_perp_site/=glbSize(nu);
	}
    get_loclx_and_rank_of_coord(iloc_lx,irank_lx,c);
  }
  
  //remap to locd
  void remap_lx_vector_to_locd(void *out,void *in,int nbytes,const Dir& mu)
  {
    const int64_t m=max_locd_perp_size_per_dir(mu);
    if((int64_t)(int32_t)(m)!=m)
      crash("integer overflow");
    
    int pars[2]={mu(),(int)m};
    if(remap_lx_to_locd[mu()]==NULL)
      remap_lx_to_locd[mu()]=new vector_remap_t(locVol(),index_make_loc_dir,pars);
    remap_lx_to_locd[mu()]->remap(out,in,nbytes);
  }
  
  void remap_locd_vector_to_lx(void *out,void *in,int nbytes,const Dir& mu)
  {
    const int64_t m=max_locd_perp_size_per_dir(mu);
    if((int64_t)(int32_t)(m)!=m)
      crash("integer overflow");
    
    int pars[2]={mu(),(int)m};
    if(remap_locd_to_lx[mu()]==NULL)
      remap_locd_to_lx[mu()]=new vector_remap_t(locd_size_per_dir(mu),index_unmake_loc_dir,pars);
    remap_locd_to_lx[mu()]->remap(out,in,nbytes);
  }
}
