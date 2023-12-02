#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <math.h>
#include <string.h>

#define EXTERN_GEOMETRY_LX
# include "geometry_lx.hpp"

#include <base/debug.hpp>
#include <base/vectors.hpp>
#include <communicate/communicate.hpp>
#include <new_types/su3.hpp>
#include <operations/remap_vector.hpp>
#include <threads/threads.hpp>

namespace nissa
{
  //Return the index of site of coord x in the border mu,nu
  int edgelx_of_coord(const coords_t &x, const int &mu, const int &nu)
  {
    int ilx=0;
    
    for(int rho=0;rho<NDIM;rho++)
      if(rho!=mu && rho!=nu)
	ilx=ilx*locSize[rho]+x[rho];
    
    return ilx;
  }
  
  //Return the index of site of coord x in the border mu
  int bordlx_of_coord(const coords_t &x, const int &mu)
  {
    int ilx=0;
    for(int nu=0;nu<NDIM;nu++)
      if(nu!=mu)
	ilx=ilx*locSize[nu]+x[nu];
    
    return ilx;
  }
  
  //Return the index of site of coord x in a box of sides s
  CUDA_HOST_AND_DEVICE int lx_of_coord(const coords_t &x, const coords_t &s)
  {
    int ilx=0;
    
    for(int mu=0;mu<NDIM;mu++)
      ilx=ilx*s[mu]+x[mu];
    
    return ilx;
  }
  
  coords_t coord_of_lx(int ilx, const coords_t s)
  {
    coords_t x;
    
    for(int mu=NDIM-1;mu>=0;mu--)
      {
	x[mu]=ilx%s[mu];
	ilx/=s[mu];
      }
    
    return x;
  }
  
  //return the volume of a given box
  int vol_of_lx(const coords_t &size)
  {
    int vol=1;
    
    for(int mu=0;mu<NDIM;mu++)
      vol*=size[mu];
    
    return vol;
  }
  
  //wrappers
  CUDA_HOST_AND_DEVICE int loclx_of_coord(const coords_t& x)
  {
    return lx_of_coord(x,locSize);
  }
  
  //wrappers
  int glblx_of_coord(const coords_t& x)
  {
    return lx_of_coord(x,glbSizes);
  }
  
  int glblx_of_coord_list(int a,int b,int c,int d)
  {
    coords_t co={a,b,c,d};
    
    return glblx_of_coord(co);
  }
  
  //combine two points
  int glblx_of_comb(int b,int wb,int c,int wc)
  {
    coords_t co;
    
    for(int mu=0;mu<NDIM;mu++)
      {
	co[mu]=glbCoordOfLoclx[b][mu]*wb+glbCoordOfLoclx[c][mu]*wc;
	while(co[mu]<0) co[mu]+=glbSizes[mu];
	co[mu]%=glbSizes[mu];
      }
    
    return glblx_of_coord(co);
  }
  
  coords_t glb_coord_of_glblx(int gx)
  {
    coords_t x;
    
    for(int mu=NDIM-1;mu>=0;mu--)
      {
	int next=gx/glbSizes[mu];
	x[mu]=gx-next*glbSizes[mu];
	gx=next;
      }
    
    return x;
  }
  
  int glblx_of_diff(const int& b,const int& c)
  {
    return glblx_of_comb(b,+1,c,-1);
  }
  
  int glblx_of_summ(const int& b,const int& c)
  {
    return glblx_of_comb(b,+1,c,+1);
  }
  
  int glblx_opp(const int& b)
  {
    return glblx_of_diff(0,b);
  }
  
  //Return the coordinate of the rank containing the global coord
  coords_t rank_coord_of_site_of_coord(const coords_t& glb_coord)
  {
    coords_t rank_coord;
    
    for(int mu=0;mu<NDIM;mu++)
      rank_coord[mu]=glb_coord[mu]/locSize[mu];
    
    return rank_coord;
  }
  
  //Return the rank of passed coord
  int rank_of_coord(const coords_t& x)
  {
    return lx_of_coord(x,nrank_dir);
  }
  
  coords_t coord_of_rank(const int& x)
  {
    return coord_of_lx(x,nrank_dir);
  }
  
  //Return the rank containing the global coordinates
  int rank_hosting_site_of_coord(const coords_t& x)
  {
    const coords_t p=rank_coord_of_site_of_coord(x);
    
    return rank_of_coord(p);
  }
  //Return the rank containing the glblx passed
  int rank_hosting_glblx(const int& gx)
  {
    const coords_t c=glb_coord_of_glblx(gx);
    
    return rank_hosting_site_of_coord(c);
  }
  
  //return the index of the site of passed "pseudolocal" coordinate
  //if the coordinates are local, return the index according to the function loclx_of_coord
  //if exactly one of the coordinate is just out return its index according to bordlx_of_coord, incremented of previous border and loc_vol
  //if exactly two coordinates are outside, return its index according to edgelx_of_coord, incremented as before stated
  int full_lx_of_coords(const coords_t ext_x)
  {
    //pseudo-localize it
    coords_t x;
    for(int mu=0;mu<NDIM;mu++)
      {
	x[mu]=ext_x[mu];
	while(x[mu]<0) x[mu]+=glbSizes[mu];
	while(x[mu]>=glbSizes[mu]) x[mu]-=glbSizes[mu];
      }
    
    //check locality
    int isloc=1;
    for(int mu=0;mu<NDIM;mu++)
      {
	isloc&=(x[mu]>=0);
	isloc&=(x[mu]<locSize[mu]);
      }
    
    if(isloc) return loclx_of_coord(x);
    
    //check borderity
    coords_t is_bord;
    for(int mu=0;mu<NDIM;mu++)
      {
	is_bord[mu]=0;
	if(is_dir_parallel[mu])
	  {
	    if(x[mu]==glbSizes[mu]-1) is_bord[mu]=-1;
	    if(x[mu]==locSize[mu]) is_bord[mu]=+1;
	  }
      }
    
    //check if it is in one of the NDIM forward or backward borders
    for(int mu=0;mu<NDIM;mu++)
      {
	bool is=is_bord[mu];
	for(int inu=0;inu<NDIM-1;inu++) is&=(is_bord[perp_dir[mu][inu]]==0);
	
	if(is)
	  {
	    if(is_bord[mu]==-1) return locVol+bord_offset[mu]+bordlx_of_coord(x,mu);             //backward border comes first
	    if(is_bord[mu]==+1) return locVol+bord_vol/2+bord_offset[mu]+bordlx_of_coord(x,mu);  //forward border comes after
	    crash("if is bord should not arrive here %d %d %d %d",ext_x[0],ext_x[1],ext_x[2],ext_x[3]);
	  }
      }
    
    //check if it is in one of the NDIM*(NDIM-1)/2 --,-+,+-,++ edges
    for(int mu=0;mu<NDIM;mu++)
      for(int inu=0;inu<NDIM-1;inu++)
	{
	  int nu=perp_dir[mu][inu];
	  
	  //order mu,nu
	  int al=(mu<nu)?mu:nu;
	  int be=(mu>nu)?mu:nu;
	  
	  bool is=is_bord[mu]&&is_bord[nu];
#if NDIM>=3
	  for(int irho=0;irho<NDIM-2;irho++) is&=(is_bord[perp2_dir[mu][inu][irho]]==0);
#endif
	  
	  if(is)
	    {
	      int iedge=edge_numb[mu][nu];
	      if((is_bord[al]==-1)&&(is_bord[be]==-1)) return locVol+bord_vol+edge_offset[iedge]+0*edge_vol/4+edgelx_of_coord(x,mu,nu);
	      if((is_bord[al]==-1)&&(is_bord[be]==+1)) return locVol+bord_vol+edge_offset[iedge]+1*edge_vol/4+edgelx_of_coord(x,mu,nu);
	      if((is_bord[al]==+1)&&(is_bord[be]==-1)) return locVol+bord_vol+edge_offset[iedge]+2*edge_vol/4+edgelx_of_coord(x,mu,nu);
	      if((is_bord[al]==+1)&&(is_bord[be]==+1)) return locVol+bord_vol+edge_offset[iedge]+3*edge_vol/4+edgelx_of_coord(x,mu,nu);
	      crash("Edge: %d, mu=%d, nu=%d %d %d %d %d",iedge,mu,nu,ext_x[0],ext_x[1],ext_x[2],ext_x[3]);
	    }
	}
    
    return -1;
  }
  
  //return the border site adiacent at surface
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  int bordlx_of_surflx(const int& loclx,
		       const int& mu)
  {
    if(not is_dir_parallel[mu])
      return -1;
    
#ifndef COMPILING_FOR_DEVICE
    if(locSize[mu]<2)
      crash("not working if one dir is smaller than 2");
#endif
    
    if(locCoordOfLoclx[loclx][mu]==0)
      return loclxNeighdw[loclx][mu]-locVol;
    
    if(locCoordOfLoclx[loclx][mu]==locSize[mu]-1)
      return loclxNeighup[loclx][mu]-locVol;
    
    return -1;
  }
  
  /// Return the edge site adiacent at a border
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  int edgelx_of_surflx(const int& loclx,
		       const int& iEdge)
  {
    const auto [mu,nu]=edge_dirs[iEdge];
    
    if(not (is_dir_parallel[mu] and is_dir_parallel[nu]))
	return -1;
    
#ifndef COMPILING_FOR_DEVICE
    if(locSize[mu]<2 or locSize[nu]<2)
      crash("not working if one dir is smaller than 2");
#endif
    
    auto iter=
      [&loclx](auto& iter,
	       const int& s,
	       const int& dir,
	       const auto&...tail)
      {
	const int c=
	  locCoordOfLoclx[loclx][dir];
	
	if(c%(locSize[dir]-1)==0)
	  {
	    const int& m=
	      loclx_neigh[c!=0][s][dir];
	    
	    if constexpr(sizeof...(tail)==0)
	      return m-locVol-bord_vol;
	    else
	      return iter(iter,m,tail...);
	  }
	else
	  return -1;
      };
    
    return iter(iter,loclx,mu,nu);
  }
  
  //label all the sites: bulk, border and edge
  void label_all_sites()
  {
    //defined a box extending over borders/edges
    int extended_box_vol=1;
    coords_t extended_box_size;
    for(int mu=0;mu<NDIM;mu++)
      {
	extended_box_size[mu]=is_dir_parallel[mu]*2+locSize[mu];
	extended_box_vol*=extended_box_size[mu];
      }
    
    for(int ivol=0;ivol<extended_box_vol;ivol++)
      {
	//subtract by one if dir is parallelized
	coords_t x=coord_of_lx(ivol,extended_box_size);
	for(int mu=0;mu<NDIM;mu++) if(is_dir_parallel[mu]) x[mu]--;
	
	//check if it is defined
	int iloc=full_lx_of_coords(x);
	if(iloc!=-1)
	  {
	    //compute global coordinates, assigning
	    for(int nu=0;nu<NDIM;nu++)
	      glbCoordOfLoclx[iloc][nu]=(x[nu]+rank_coord[nu]*locSize[nu]+glbSizes[nu])%glbSizes[nu];
	    
	    //find the global index
	    int iglb=glblx_of_coord(glbCoordOfLoclx[iloc]);
	    
	    //if it is on the bulk store it
	    if(iloc<locVol)
	      {
		for(int nu=0;nu<NDIM;nu++) locCoordOfLoclx[iloc][nu]=x[nu];
		glblxOfLoclx[iloc]=iglb;
	      }
	    
	    //if it is on the border store it
	    if(iloc>=locVol&&iloc<locVol+bord_vol)
	      {
		int ibord=iloc-locVol;
		glblxOfBordlx[ibord]=iglb;
		loclxOfBordlx[ibord]=iloc;
	      }
	    
	    //if it is on the edge store it
	    if(iloc>=locVol+bord_vol)
	      {
		int iedge=iloc-locVol-bord_vol;
		glblxOfEdgelx[iedge]=iglb;
		//loclx_of_edgelx[iedge]=iloc;
	      }
	  }
      }
  }
  
  //find the neighbours
  void find_neighbouring_sites()
  {
    //loop over the four directions
    for(int ivol=0;ivol<locVol+bord_vol+edge_vol;ivol++)
      for(int mu=0;mu<NDIM;mu++)
	{
	  //copy the coords
	  coords_t n;
	  for(int nu=0;nu<NDIM;nu++) n[nu]=glbCoordOfLoclx[ivol][nu]-locSize[nu]*rank_coord[nu];
	  
	  //move forward
	  n[mu]++;
	  int nup=full_lx_of_coords(n);
	  //move backward
	  n[mu]-=2;
	  int ndw=full_lx_of_coords(n);
	  
	  //if "local" assign it (automatically -1 otherwise)
	  loclxNeighup[ivol][mu]=nup;
	  loclxNeighdw[ivol][mu]=ndw;
	}
  }
  
  //finds how to fill the borders with opposite surface (up b->dw s)
  void findSurfOfBord()
  {
    PAR(0,locVol,
	CAPTURE(),loclx,
	{
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      const int bordlx=bordlx_of_surflx(loclx,mu);
	      if(bordlx!=-1) surflxOfBordlx[bordlx]=loclx;
	    }
	});
  }
  
  //finds how to fill the borders with opposite surface (up b->dw s)
  void findSurfOfEdges()
  {
    PAR(0,locVol,
	CAPTURE(),loclx,
	{
	  for(int iEdge=0;iEdge<nEdges;iEdge++)
	    {
	      const int edgeLx=edgelx_of_surflx(loclx,iEdge);
	      
	      if(edgeLx!=-1)
		surflxOfEdgelx[edgeLx]=loclx;
	    }
	});
  }
  
  //index all the sites on bulk
  void find_bulk_sites()
  {
    //check surfacity
    int ibulk=0,inon_fw_surf=0,inon_bw_surf=0;
    int isurf=0,ifw_surf=0,ibw_surf=0;
    NISSA_LOC_VOL_LOOP(ivol)
    {
      //find if it is on bulk or non_fw or non_bw surf
      int is_bulk=true,is_non_fw_surf=true,is_non_bw_surf=true;
      for(int mu=0;mu<NDIM;mu++)
	if(is_dir_parallel[mu])
	  {
	    if(locCoordOfLoclx[ivol][mu]==locSize[mu]-1) is_bulk=is_non_fw_surf=false;
	    if(locCoordOfLoclx[ivol][mu]==0)              is_bulk=is_non_bw_surf=false;
	  }
      
      //mark it
      if(is_bulk) loclxOfBulklx[ibulk++]=ivol;
      else        loclxOfSurflx[isurf++]=ivol;
      if(is_non_fw_surf) loclxOfNonFwSurflx[inon_fw_surf++]=ivol;
      else               loclxOfFwSurflx[ifw_surf++]=ivol;
      if(is_non_bw_surf) loclxOfNonBwSurflx[inon_bw_surf++]=ivol;
      else               loclxOfBwSurflx[ibw_surf++]=ivol;
    }
    
    if(ibulk!=bulkVol) crash("mismatch in bulk id");
    if(isurf!=surfVol) crash("mismatch in surf id");
    if(inon_fw_surf!=nonFwSurfVol) crash("mismatch in non_fw_surf id");
    if(inon_bw_surf!=nonBwSurfVol) crash("mismatch in non_bw_surf id");
    if(ifw_surf!=fwSurfVol) crash("mismatch in fw_surf id");
    if(ibw_surf!=bwSurfVol) crash("mismatch in bw_surf id");
  }
  
  //indexes run as t,x,y,z (faster:z)
  void set_lx_geometry()
  {
    if(lxGeomInited==1) crash("cartesian geometry already intialized!");
    lxGeomInited=1;
    
    if(grid_inited!=1) crash("grid not initialized!");
    
    //find the rank of the neighbour in the various dir
    for(int mu=0;mu<NDIM;mu++)
      {
	do
	  {
	    coords_t temp=rank_coord;
	    temp[mu]=(temp[mu]+nrank_dir[mu]-1)%nrank_dir[mu];
	    rank_neighdw[mu]=rank_of_coord(temp);
	  }
	while(0);
	
	do
	  {
	    coords_t temp=rank_coord;
	    temp[mu]=(temp[mu]+1)%nrank_dir[mu];
	    rank_neighup[mu]=rank_of_coord(temp);
	  }
	while(0);
      }
    rank_neigh[0]=rank_neighdw;
    rank_neigh[1]=rank_neighup;
    
    locCoordOfLoclx=nissa_malloc("loc_coord_of_loclx",locVol,coords_t);
    glbCoordOfLoclx=nissa_malloc("glb_coord_of_loclx",locVol+bord_vol+edge_vol,coords_t);
    loclx_neigh[0]=loclxNeighdw=nissa_malloc("loclx_neighdw",locVol+bord_vol+edge_vol,coords_t);
    loclx_neigh[1]=loclxNeighup=nissa_malloc("loclx_neighup",locVol+bord_vol+edge_vol,coords_t);
    
    //local to global
    glblxOfLoclx=nissa_malloc("glblx_of_loclx",locVol,int);
    
    //borders
    glblxOfBordlx=nissa_malloc("glblx_of_bordlx",bord_vol,int);
    loclxOfBordlx=nissa_malloc("loclx_of_bordlx",bord_vol,int);
    surflxOfBordlx=nissa_malloc("surflx_of_bordlx",bord_vol,int);
    
    //bulk and surfs
    loclxOfBulklx=nissa_malloc("loclx_of_bulklx",bulkVol,int);
    loclxOfSurflx=nissa_malloc("loclx_of_surflx",surfVol,int);
    loclxOfNonBwSurflx=nissa_malloc("loclx_of_non_bw_surflx",nonBwSurfVol,int);
    loclxOfNonFwSurflx=nissa_malloc("loclx_of_non_fw_surflx",nonFwSurfVol,int);
    loclxOfBwSurflx=nissa_malloc("loclx_of_bw_surflx",bwSurfVol,int);
    loclxOfFwSurflx=nissa_malloc("loclx_of_fw_surflx",fwSurfVol,int);
    
    //edges
    glblxOfEdgelx=nissa_malloc("glblx_of_edgelx",edge_vol,int);
    surflxOfEdgelx=nissa_malloc("bordlx_of_edgelx",edge_vol,int);
    
    //label the sites and neighbours
    label_all_sites();
    find_neighbouring_sites();
    
    //matches surface and opposite border and edges
    findSurfOfBord();
    findSurfOfEdges();
    
    //find bulk sites
    find_bulk_sites();
    
    //allocate a buffer large enough to allow communications of su3spinspin lx border
    recv_buf_size=std::max(recv_buf_size,bord_vol*sizeof(su3spinspin));
    send_buf_size=std::max(send_buf_size,bord_vol*sizeof(su3spinspin));
    
    master_printf("Cartesian geometry intialized\n");
  }
  
  //global movements
  int glblx_neighup(const int& gx,const int& mu)
  {
    coords_t c=glb_coord_of_glblx(gx);
    c[mu]=(c[mu]+1)%glbSizes[mu];
    
    return glblx_of_coord(c);
  }
  
  int glblx_neighdw(const int& gx,const int& mu)
  {
    coords_t c=glb_coord_of_glblx(gx);
    c[mu]=(c[mu]+glbSizes[mu]-1)%glbSizes[mu];
    
    return glblx_of_coord(c);
  }
  
  //unset cartesian geometry
  void unset_lx_geometry()
  {
    if(lxGeomInited!=1) crash("cartesian geometry not initialized!");
    
    master_printf("Unsetting cartesian geometry\n");
    lxGeomInited=0;
    
    nissa_free(recv_buf);
    nissa_free(send_buf);
    
    nissa_free(locCoordOfLoclx);
    nissa_free(glbCoordOfLoclx);
    nissa_free(loclxNeighup);
    nissa_free(loclxNeighdw);
    
    nissa_free(glblxOfLoclx);
    nissa_free(glblxOfBordlx);
    nissa_free(loclxOfBordlx);
    nissa_free(surflxOfBordlx);
    nissa_free(glblxOfEdgelx);
    nissa_free(surflxOfEdgelx);
    
    nissa_free(loclxOfBulklx);
    nissa_free(loclxOfSurflx);
    nissa_free(loclxOfNonFwSurflx);
    nissa_free(loclxOfFwSurflx);
    nissa_free(loclxOfNonBwSurflx);
    nissa_free(loclxOfBwSurflx);
  }
}
