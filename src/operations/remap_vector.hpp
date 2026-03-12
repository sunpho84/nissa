#ifndef _REMAP_VECTOR_HPP
#define _REMAP_VECTOR_HPP

#include <base/debug.hpp>
#include <communicate/all_to_all.hpp>

namespace nissa
{
  struct vector_remap_t;
  
  //local direction geometry
  inline int max_locd_perp_size_per_dir[NDIM];
  
  inline int locd_perp_size_per_dir[NDIM];
  
  inline int max_locd_size;
  
  inline int locd_size_per_dir[NDIM];
  
  inline vector_remap_t* _remapLxToLocd[NDIM];
  
  inline vector_remap_t* _remap_locd_to_lx[NDIM];
  
  struct vector_remap_t :
    all_to_all_comm_t
  {
    vector_remap_t(all_to_all_comm_t&& oth) :
      all_to_all_comm_t(std::move(oth))
    {
    }
    
    // /// Default initialization
    // vector_remap_t()=default;
    
    /// Initializes the remap, with a lambda
    template <typename F>
    vector_remap_t(const int64_t& nel_fr,
		   F&& index)
    {
      all_to_all_scattering_list_t sl;
      for(int64_t iel_fr=0;iel_fr<nel_fr;iel_fr++)
	{
	  const auto [rank_to,iel_to]=index(iel_fr);
	  
	  if(rank_to>=nranks or rank_to<0)
	    CRASH("destination rank %d does not exist!",rank_to);
	  
	  sl.push_back(std::make_pair(iel_fr,iel_to*nranks+rank_to));
	}
      
      setup_knowing_where_to_send(sl);
    }
    
    vector_remap_t inverse() const
    {
      return ((const all_to_all_comm_t*)this)->inverse();
    }
    
    void remap(void* out,
	       const void* in,
	       const size_t& bps) const
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
  
  /// Return an indexer to make dir mu local
  inline auto getLocDirIndexMaker(const int& mu,
				  const int64_t& prp_max_vol)
  {
    return
      [mu,prp_max_vol](const int64_t& iloc_lx)
      {
	int glb_perp_site=0;
	for(int nu=0;nu<NDIM;nu++)
	  if(mu!=nu)
	    glb_perp_site=glb_perp_site*glbSize[nu]+glbCoordOfLoclx[iloc_lx][nu];
	
	const int irank_locld=glb_perp_site/prp_max_vol;
	
	int64_t iloc_locld=glb_perp_site-irank_locld*prp_max_vol;
	iloc_locld=iloc_locld*glbSize[mu]+glbCoordOfLoclx[iloc_lx][mu];
	
	return std::make_pair(irank_locld,iloc_locld);
      };
  }
  
  /// Gets the lxToLocd remapper allocating it if not exists
  inline const vector_remap_t& lxToLocdRemapper(const int& mu)
  {
    if(_remapLxToLocd[mu]==nullptr)
      _remapLxToLocd[mu]=
	new vector_remap_t(locVol,getLocDirIndexMaker(mu,max_locd_perp_size_per_dir[mu]));
    
    return *_remapLxToLocd[mu];
  }
  
  /// Remap to locd
  inline void remapLxVectorToLocd(void* out,
				  const void* in,
				  const int& nbytes,
				  const int& mu)
  {
    lxToLocdRemapper(mu).remap(out,in,nbytes);
  }
  
  /// Provide the index to unmake a local dir
  inline auto getLocDirIndexUnmaker(const int& mu,
				    const int64_t& prp_max_vol)
  {
    return
      [mu,prp_max_vol](int64_t iloc_locld) // don't make constant
    {
      Coords c;
      c[mu]=iloc_locld%glbSize[mu];
      iloc_locld/=glbSize[mu];
      
      int64_t glb_perp_site=
	iloc_locld+rank*prp_max_vol;
      
      for(int nu=NDIM-1;nu>=0;nu--)
	if(mu!=nu)
	  {
	    c[nu]=glb_perp_site%glbSize[nu];
	    glb_perp_site/=glbSize[nu];
	  }
      
      //int &irank_lx,int &iloc_lx;
      return getLoclxAndRankOfCoords(c);
    };
  }
  
  /// Gets the locdToLx remapper allocating it if not exists
  inline vector_remap_t& locdToLxRemapper(const int& mu)
  {
    if(_remapLxToLocd[mu]==nullptr)
      _remapLxToLocd[mu]=
	new vector_remap_t(lxToLocdRemapper(mu).inverse());
    
    return *_remapLxToLocd[mu];
  }
  
  inline void remapLocdVectorToLx(void* out,
				  const void* in,
				  const int& nbytes,
				  const int& mu)
  {
    locdToLxRemapper(mu).remap(out,in,nbytes);
  }
}

#endif
