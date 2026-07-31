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
  
  inline vector_remap_t* _remapLocdToLx[NDIM];
  
  inline vector_remap_t* _topo_corr_rem;
  
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
  inline const vector_remap_t& locdToLxRemapper(const int& mu)
  {
    if(_remapLocdToLx[mu]==nullptr)
      _remapLocdToLx[mu]=
	new vector_remap_t(locd_size_per_dir[mu],getLocDirIndexUnmaker(mu,max_locd_perp_size_per_dir[mu]));//lxToLocdRemapper(mu).inverse());
    
    return *_remapLocdToLx[mu];
  }
  
  inline void remapLocdVectorToLx(void* out,
				  const void* in,
				  const int& nbytes,
				  const int& mu)
  {
    locdToLxRemapper(mu).remap(out,in,nbytes);
  }
  
  /// Finding the index to put only 1/16 of the data
  inline std::pair<int,int64_t> index_to_topo_corr_remapping(const int& iloc_lx)
  {
    int subcube=0,subcube_el=0;
    int subcube_size[NDIM][2],subcube_coord[NDIM],subcube_el_coord[NDIM];
    for(int mu=0;mu<NDIM;mu++)
      {
	subcube_size[mu][0]=glbSize[mu]/2+1;
	subcube_size[mu][1]=glbSize[mu]/2-1;
	
	//take global coord and identify subcube
	int glx_mu=glbCoordOfLoclx[iloc_lx][mu];
	subcube_coord[mu]=(glx_mu>=subcube_size[mu][0]);
	subcube=subcube*2+subcube_coord[mu];
	
	//identify also the local coord
	subcube_el_coord[mu]=glx_mu-subcube_coord[mu]*subcube_size[mu][0];
	subcube_el=subcube_el*subcube_size[mu][subcube_coord[mu]]+subcube_el_coord[mu];
      }
    
    //summ the smaller-index cubes
    Coords nsubcubes_per_dir;
    for(int mu=0;mu<NDIM;mu++) nsubcubes_per_dir[mu]=2;
    int64_t minind_cube_vol=0;
    for(int isubcube=0;isubcube<subcube;isubcube++)
      {
	//get coords
	Coords c=coordOfLx(isubcube,nsubcubes_per_dir);
	//compute vol
	int64_t subcube_vol=1;
	for(int mu=0;mu<NDIM;mu++)
	  subcube_vol*=subcube_size[mu][c[mu]];
	minind_cube_vol+=subcube_vol;
      }
    
    const int64_t tmp=
      subcube_el+minind_cube_vol;
    
    /// Destintation rank
    const int irank=
      tmp/locVol;
    
    /// Dest loclx
    const int64_t iloc=
      tmp%locVol;
    
    return std::pair{irank,iloc};
  }
  
  /// Gets the topo_corr_rem remapper
  inline const vector_remap_t& topoCorrRem()
  {
    if(_topo_corr_rem==nullptr)
      _topo_corr_rem=
	new vector_remap_t(locVol,index_to_topo_corr_remapping);
    
    return *_topo_corr_rem;
  }
}

#endif
