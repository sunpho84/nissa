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
#include <operations/su3_paths/gauge_sweeper.hpp>
#include <routines/ios.hpp>
#include <routines/mpi_routines.hpp>
#include <threads/threads.hpp>

namespace nissa
{
  //Return the index of site of coord x in the border mu,nu
  int64_t edgelxOfCoord(const Coords &x,
			const int &mu,
			const int &nu)
  {
    int64_t ilx=0;
    
    for(int rho=0;rho<NDIM;rho++)
      if(rho!=mu and rho!=nu)
	ilx=ilx*locSize[rho]+x[rho];
    
    return ilx;
  }
  
  //Return the index of site of coord x in the border mu
  int64_t bordlxOfCoord(const Coords &x,
			const int &mu)
  {
    int64_t ilx=0;
    
    for(int nu=0;nu<NDIM;nu++)
      if(nu!=mu)
	ilx=ilx*locSize[nu]+x[nu];
    
    return ilx;
  }
  
  //Return the index of site of coord x in a box of sides s
  CUDA_HOST_AND_DEVICE int64_t lxOfCoord(const Coords &x,
					 const Coords &s)
  {
    int64_t ilx=0;
    
    for(int mu=0;mu<NDIM;mu++)
      ilx=ilx*s[mu]+x[mu];
    
    return ilx;
  }
  
  Coords coordOfLx(const int64_t& _ilx,
		   const Coords& s)
  {
    int64_t ilx=
      _ilx;
    
    Coords x;
    for(int mu=NDIM-1;mu>=0;mu--)
      {
	x[mu]=ilx%s[mu];
	ilx/=s[mu];
      }
    
    return x;
  }
  
  //return the volume of a given box
  int64_t volOfLx(const Coords& size)
  {
    int64_t vol=1;
    
    for(int mu=0;mu<NDIM;mu++)
      vol*=size[mu];
    
    return vol;
  }
  
  //wrappers
  CUDA_HOST_AND_DEVICE int64_t loclxOfCoord(const Coords& x)
  {
    return lxOfCoord(x,locSize);
  }
  
  //wrappers
  int64_t glblxOfCoord(const Coords& x)
  {
    return lxOfCoord(x,glbSize);
  }
  
  int64_t glblxOfCoordList(const int& a,
			   const int& b,
			   const int& c,
			   const int& d)
    
  {
    return glblxOfCoord(Coords{a,b,c,d});
  }
  
  //combine two points
  int64_t glblxOfComb(const int64_t& b,
		      const int& wb,
		      const int64_t& c,
		      const int& wc)
  {
    Coords co;
    
    for(int mu=0;mu<NDIM;mu++)
      {
	co[mu]=glbCoordOfLoclx[b][mu]*wb+glbCoordOfLoclx[c][mu]*wc;
	while(co[mu]<0)
	  co[mu]+=glbSize[mu];
	co[mu]%=glbSize[mu];
      }
    
    return glblxOfCoord(co);
  }
  
  Coords glbCoordOfGlblx(const int64_t& _gx)
  {
    int64_t gx=
      _gx;
    
    Coords x;
    for(int mu=NDIM-1;mu>=0;mu--)
      {
	int next=gx/glbSize[mu];
	x[mu]=gx-next*glbSize[mu];
	gx=next;
      }
    
    return x;
  }
  
  int64_t glblxOfDiff(const int64_t& b,
		      const int64_t& c)
  {
    return glblxOfComb(b,+1,c,-1);
  }
  
  int64_t glblxOfSum(const int64_t& b,
		     const int64_t& c)
  {
    return glblxOfComb(b,+1,c,+1);
  }
  
  int64_t glblxOpp(const int64_t& b)
  {
    return glblxOfDiff(0l,b);
  }
  
  //Return the coordinate of the rank containing the global coord
  Coords rankCoordsOfSiteOfCoord(const Coords& glb_coord)
  {
    Coords rankCoord;
    
    for(int mu=0;mu<NDIM;mu++)
      rankCoord[mu]=glb_coord[mu]/locSize[mu];
    
    return rankCoord;
  }
  
  //Return the rank of passed coord
  int rankOfCoords(const Coords& x)
  {
    return lxOfCoord(x,nRanksDir);
  }
  
  Coords coordOfRank(const int& x)
  {
    return coordOfLx(x,nRanksDir);
  }
  
  //Return the rank containing the global coordinates
  int rankHostingSiteOfCoords(const Coords& x)
  {
    const Coords p=
      rankCoordsOfSiteOfCoord(x);
    
    return rankOfCoords(p);
  }
  //Return the rank containing the glblx passed
  int rankHostingGlblx(const int64_t& gx)
  {
    const Coords c=
      glbCoordOfGlblx(gx);
    
    return rankHostingSiteOfCoords(c);
  }
  
  //return the index of the site of passed "pseudolocal" coordinate
  //if the coordinates are local, return the index according to the function loclx_of_coord
  //if exactly one of the coordinate is just out return its index according to bordlx_of_coord, incremented of previous border and loc_vol
  //if exactly two coordinates are outside, return its index according to edgelx_of_coord, incremented as before stated
  int64_t fullLxOfCoords(const Coords& ext_x)
  {
    //pseudo-localize it
    Coords x;
    for(int mu=0;mu<NDIM;mu++)
      {
	x[mu]=ext_x[mu];
	while(x[mu]<0)
	  x[mu]+=glbSize[mu];
	while(x[mu]>=glbSize[mu])
	  x[mu]-=glbSize[mu];
      }
    
    //check locality
    bool isloc=true;
    for(int mu=0;mu<NDIM;mu++)
      {
	isloc&=(x[mu]>=0);
	isloc&=(x[mu]<locSize[mu]);
      }
    
    if(isloc) return loclxOfCoord(x);
    
    //check borderity
    Coords is_bord;
    for(int mu=0;mu<NDIM;mu++)
      {
	is_bord[mu]=0;
	if(isDirParallel[mu])
	  {
	    if(x[mu]==glbSize[mu]-1) is_bord[mu]=-1;
	    if(x[mu]==locSize[mu]) is_bord[mu]=+1;
	  }
      }
    
    //check if it is in one of the NDIM forward or backward borders
    for(int mu=0;mu<NDIM;mu++)
      {
	bool is=is_bord[mu];
	for(int inu=0;inu<NDIM-1;inu++)
	  is&=(is_bord[perpDirs[mu][inu]]==0);
	
	if(is)
	  {
	    if(is_bord[mu]==-1) return locVol+bordOffset[mu]+bordlxOfCoord(x,mu);             //backward border comes first
	    if(is_bord[mu]==+1) return locVol+bordVol/2+bordOffset[mu]+bordlxOfCoord(x,mu);  //forward border comes after
	    crash("if is bord should not arrive here %d %d %d %d",ext_x[0],ext_x[1],ext_x[2],ext_x[3]);
	  }
      }
    
    //check if it is in one of the NDIM*(NDIM-1)/2 --,-+,+-,++ edges
    for(int mu=0;mu<NDIM;mu++)
      for(int inu=0;inu<NDIM-1;inu++)
	{
	  int nu=perpDirs[mu][inu];
	  
	  //order mu,nu
	  int al=(mu<nu)?mu:nu;
	  int be=(mu>nu)?mu:nu;
	  
	  bool is=is_bord[mu]&&is_bord[nu];
#if NDIM>=3
	  for(int irho=0;irho<NDIM-2;irho++) is&=(is_bord[perp2Dirs[mu][inu][irho]]==0);
#endif
	  
	  if(is)
	    {
	      int iedge=edge_numb[mu][nu];
	      if((is_bord[al]==-1)&&(is_bord[be]==-1)) return locVol+bordVol+edge_offset[iedge]+0*edgeVol/4+edgelxOfCoord(x,mu,nu);
	      if((is_bord[al]==-1)&&(is_bord[be]==+1)) return locVol+bordVol+edge_offset[iedge]+1*edgeVol/4+edgelxOfCoord(x,mu,nu);
	      if((is_bord[al]==+1)&&(is_bord[be]==-1)) return locVol+bordVol+edge_offset[iedge]+2*edgeVol/4+edgelxOfCoord(x,mu,nu);
	      if((is_bord[al]==+1)&&(is_bord[be]==+1)) return locVol+bordVol+edge_offset[iedge]+3*edgeVol/4+edgelxOfCoord(x,mu,nu);
	      crash("Edge: %d, mu=%d, nu=%d %d %d %d %d",iedge,mu,nu,ext_x[0],ext_x[1],ext_x[2],ext_x[3]);
	    }
	}
    
    return -1;
  }
  
  //return the border site adiacent at surface
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  int64_t bordlx_of_surflx(const int64_t& loclx,
			   const int& mu)
  {
    if(not isDirParallel[mu])
      return -1;
    
#ifndef COMPILING_FOR_DEVICE
    if(locSize[mu]<2) //nasty
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
  int64_t edgelx_of_surflx(const int64_t& loclx,
			   const int& iEdge)
  {
    const auto [mu,nu]=edge_dirs[iEdge];
    
    if(not (isDirParallel[mu] and isDirParallel[nu]))
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
	const int64_t c=
	  locCoordOfLoclx[loclx][dir];
	
	if(c%(locSize[dir]-1)==0)
	  {
	    const int64_t& m=
	      loclx_neigh[c!=0][s][dir];
	    
	    if constexpr(sizeof...(tail)==0)
	      return m-locVol-bordVol;
	    else
	      return iter(iter,m,tail...);
	  }
	else
	  return (int64_t)-1;
      };
    
    return iter(iter,loclx,mu,nu);
  }
  
  //label all the sites: bulk, border and edge
  void label_all_sites()
  {
    //defined a box extending over borders/edges
    int extended_box_vol=1;
    Coords extended_box_size;
    for(int mu=0;mu<NDIM;mu++)
      {
	extended_box_size[mu]=isDirParallel[mu]*2+locSize[mu];
	extended_box_vol*=extended_box_size[mu];
      }
    
    for(int64_t ivol=0;ivol<extended_box_vol;ivol++)
      {
	//subtract by one if dir is parallelized
	Coords x=coordOfLx(ivol,extended_box_size);
	for(int mu=0;mu<NDIM;mu++) if(isDirParallel[mu]) x[mu]--;
	
	//check if it is defined
	int64_t iloc=fullLxOfCoords(x);
	if(iloc!=-1)
	  {
	    //compute global coordinates, assigning
	    for(int nu=0;nu<NDIM;nu++)
	      glbCoordOfLoclx[iloc][nu]=(x[nu]+rankCoord[nu]*locSize[nu]+glbSize[nu])%glbSize[nu];
	    
	    //find the global index
	    int64_t iglb=glblxOfCoord(glbCoordOfLoclx[iloc]);
	    
	    //if it is on the bulk store it
	    if(iloc<locVol)
	      {
		for(int nu=0;nu<NDIM;nu++) locCoordOfLoclx[iloc][nu]=x[nu];
		glblxOfLoclx[iloc]=iglb;
	      }
	    
	    //if it is on the border store it
	    if(iloc>=locVol&&iloc<locVol+bordVol)
	      {
		int64_t ibord=iloc-locVol;
		glblxOfBordlx[ibord]=iglb;
		loclxOfBordlx[ibord]=iloc;
	      }
	    
	    //if it is on the edge store it
	    if(iloc>=locVol+bordVol)
	      {
		int64_t iedge=iloc-locVol-bordVol;
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
    for(int64_t ivol=0;ivol<locVol+bordVol+edgeVol;ivol++)
      for(int mu=0;mu<NDIM;mu++)
	{
	  //copy the coords
	  Coords n;
	  for(int nu=0;nu<NDIM;nu++) n[nu]=glbCoordOfLoclx[ivol][nu]-locSize[nu]*rankCoord[nu];
	  
	  //move forward
	  n[mu]++;
	  int64_t nup=fullLxOfCoords(n);
	  //move backward
	  n[mu]-=2;
	  int64_t ndw=fullLxOfCoords(n);
	  
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
	      const int64_t bordlx=bordlx_of_surflx(loclx,mu);
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
	      const int64_t edgeLx=edgelx_of_surflx(loclx,iEdge);
	      
	      if(edgeLx!=-1)
		surflxOfEdgelx[edgeLx]=loclx;
	    }
	});
  }
  
  //index all the sites on bulk
  void find_bulk_sites()
  {
    //check surfacity
    int64_t ibulk=0,inon_fw_surf=0,inon_bw_surf=0;
    int64_t isurf=0,ifw_surf=0,ibw_surf=0;
    NISSA_LOC_VOL_LOOP(ivol)
      {
	//find if it is on bulk or non_fw or non_bw surf
	bool is_bulk=true,is_non_fw_surf=true,is_non_bw_surf=true;
	for(int mu=0;mu<NDIM;mu++)
	  if(isDirParallel[mu])
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
    
    if(gridInited!=1) crash("grid not initialized!");
    
    //find the rank of the neighbour in the various dir
    for(int mu=0;mu<NDIM;mu++)
      {
	do
	  {
	    Coords temp=rankCoord;
	    temp[mu]=(temp[mu]+nRanksDir[mu]-1)%nRanksDir[mu];
	    rankNeighdw[mu]=rankOfCoords(temp);
	  }
	while(0);
	
	do
	  {
	    Coords temp=rankCoord;
	    temp[mu]=(temp[mu]+1)%nRanksDir[mu];
	    rankNeighup[mu]=rankOfCoords(temp);
	  }
	while(0);
      }
    rankNeigh[0]=rankNeighdw;
    rankNeigh[1]=rankNeighup;
    
    locCoordOfLoclx=nissa_malloc("loc_coord_of_loclx",locVol,Coords);
    glbCoordOfLoclx=nissa_malloc("glb_coord_of_loclx",locVol+bordVol+edgeVol,Coords);
    loclx_neigh[0]=loclxNeighdw=nissa_malloc("loclx_neighdw",locVol+bordVol+edgeVol,Coords);
    loclx_neigh[1]=loclxNeighup=nissa_malloc("loclx_neighup",locVol+bordVol+edgeVol,Coords);
    
    //local to global
    glblxOfLoclx=nissa_malloc("glblx_of_loclx",locVol,int64_t);
    
    //borders
    glblxOfBordlx=nissa_malloc("glblx_of_bordlx",bordVol,int64_t);
    loclxOfBordlx=nissa_malloc("loclx_of_bordlx",bordVol,int64_t);
    surflxOfBordlx=nissa_malloc("surflx_of_bordlx",bordVol,int64_t);
    
    //bulk and surfs
    loclxOfBulklx=nissa_malloc("loclx_of_bulklx",bulkVol,int64_t);
    loclxOfSurflx=nissa_malloc("loclx_of_surflx",surfVol,int64_t);
    loclxOfNonBwSurflx=nissa_malloc("loclx_of_non_bw_surflx",nonBwSurfVol,int64_t);
    loclxOfNonFwSurflx=nissa_malloc("loclx_of_non_fw_surflx",nonFwSurfVol,int64_t);
    loclxOfBwSurflx=nissa_malloc("loclx_of_bw_surflx",bwSurfVol,int64_t);
    loclxOfFwSurflx=nissa_malloc("loclx_of_fw_surflx",fwSurfVol,int64_t);
    
    //edges
    glblxOfEdgelx=nissa_malloc("glblx_of_edgelx",edgeVol,int64_t);
    surflxOfEdgelx=nissa_malloc("bordlx_of_edgelx",edgeVol,int64_t);
    
    //label the sites and neighbours
    label_all_sites();
    find_neighbouring_sites();
    
    //matches surface and opposite border and edges
    findSurfOfBord();
    findSurfOfEdges();
    
    //find bulk sites
    find_bulk_sites();
    
    //allocate a buffer large enough to allow communications of su3spinspin lx border
    recv_buf_size=std::max(recv_buf_size,bordVol*sizeof(su3spinspin));
    send_buf_size=std::max(send_buf_size,bordVol*sizeof(su3spinspin));
    
    //create the sweepers but do not fully initialize
    Wilson_sweeper=new gauge_sweeper_t;
    Symanzik_sweeper=new gauge_sweeper_t;
    
    //set locd geom (one of the dimension local and fastest running, the other as usual)
    max_locd_size=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	remap_lx_to_locd[mu]=remap_locd_to_lx[mu]=NULL;
	max_locd_perp_size_per_dir[mu]=(glbVol/glbSize[mu]+nranks-1)/nranks;
	locd_perp_size_per_dir[mu]=(int)std::min((int64_t)max_locd_perp_size_per_dir[mu],glbVol/glbSize[mu]-max_locd_perp_size_per_dir[mu]*rank);
	verbosity_lv3_master_printf("rank %d locd_perp_size_per_dir[%d]: %d\n",rank,mu,locd_perp_size_per_dir[mu]);
	locd_size_per_dir[mu]=locd_perp_size_per_dir[mu]*glbSize[mu];
	max_locd_size=std::max(max_locd_size,locd_size_per_dir[mu]);
      }
    
    master_printf("Cartesian geometry intialized\n");
  }
  
  //global movements
  int64_t glblxNeighup(const int64_t& gx,
			const int& mu)
  {
    Coords c=glbCoordOfGlblx(gx);
    c[mu]=(c[mu]+1)%glbSize[mu];
    
    return glblxOfCoord(c);
  }
  
  int64_t glblxNeighdw(const int64_t& gx,
			const int& mu)
  {
    Coords c=glbCoordOfGlblx(gx);
    c[mu]=(c[mu]+glbSize[mu]-1)%glbSize[mu];
    
    return glblxOfCoord(c);
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
    
    delete Wilson_sweeper;
    delete Symanzik_sweeper;
  }
  
  //definitions of lexical ordered senders for edges
  void initialize_lx_edge_senders_of_kind(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *base)
  {
    //Various type useful for edges and sub-borders
    MPI_Datatype MPI_3_SLICE;
    MPI_Type_contiguous(locSize[3],*base,&MPI_3_SLICE);
    
    ///////////define the sender for the 6 kinds of edges////////////
    //the 01 sender, that is simply a vector of L[2] vector of L[3] 
    MPI_Type_contiguous(locSize[2]*locSize[3],*base,&(MPI_EDGE_SEND[0]));
    //the 02 sender is a vector of L[1] segment of length L[3] (already defined) separated by L[2] of them
    MPI_Type_vector(locSize[1],1,locSize[2],MPI_3_SLICE,&(MPI_EDGE_SEND[1]));
    //the 03 sender is a vector of length L[1]xL[2] of single elements, separated by L[3] of them
    MPI_Type_vector(locSize[1]*locSize[2],1,locSize[3],*base,&(MPI_EDGE_SEND[2]));
    //the 12 sender should be equal to the 02 sender, with 1->0
    MPI_Type_vector(locSize[0],1,locSize[2],MPI_3_SLICE,&(MPI_EDGE_SEND[3]));
    //the 13 sender should be equal to the 03 sender, with 1->0
    MPI_Type_vector(locSize[0]*locSize[2],1,locSize[3],*base,&(MPI_EDGE_SEND[4]));
    //the 23 sender should be equal to the 03 sender with 1<->2
    MPI_Type_vector(locSize[0]*locSize[1],1,locSize[3],*base,&(MPI_EDGE_SEND[5]));
    //Commit
    for(int iedge=0;iedge<6;iedge++) MPI_Type_commit(&(MPI_EDGE_SEND[iedge]));
  }
  
  //definitions of lexical ordered receivers for edges
  void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base)
  {
    //define the NDIM*(NDIM-1)/2 edges receivers, which are contiguous in memory
    int iedge=0;
    for(int mu=0;mu<NDIM;mu++)
      for(int nu=mu+1;nu<NDIM;nu++)
	{
	  MPI_Type_contiguous(locVol/locSize[mu]/locSize[nu],*base,&(MPI_EDGE_RECE[iedge]));
	  MPI_Type_commit(&(MPI_EDGE_RECE[iedge]));
	  iedge++;
	}
  }
  
  //initalize senders and receivers for edges of lexically ordered vectors
  void set_lx_edge_senders_and_receivers(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base)
  {
    initialize_lx_edge_senders_of_kind(MPI_EDGE_SEND,base);
    initialize_lx_edge_receivers_of_kind(MPI_EDGE_RECE,base);
  }
  
  //define all the local lattice momenta
  void define_local_momenta(Momentum* k,double *k2,Momentum* ktilde,double *ktilde2,const Momentum& bc)
  {
    if(!lxGeomInited) set_lx_geometry();
    
    //first of all, defines the local momenta for the various directions
    NISSA_LOC_VOL_LOOP(imom)
    {
      k2[imom]=ktilde2[imom]=0;
      for(int mu=0;mu<NDIM;mu++)
	{
	  k[imom][mu]=M_PI*(2*glbCoordOfLoclx[imom][mu]+bc[mu])/glbSize[mu];
	  ktilde[imom][mu]=sin(k[imom][mu]);
	  
	  k2[imom]+=k[imom][mu]*k[imom][mu];
	  ktilde2[imom]+=ktilde[imom][mu]*ktilde[imom][mu];
	}
    }
  }
  
  //return the staggered phases for a given site
  CUDA_HOST_AND_DEVICE Coords getStagphaseOfLx(const int64_t& ivol)
  {
    Coords ph;
    
    ph[0]=1;
    for(int mu=1;mu<NDIM;mu++)
      ph[mu]=ph[mu-1]*(1-2*(glbCoordOfLoclx[ivol][mu-1]%2));
    
    return ph;
  }
  
  //return the staggered phases for a given site
  CUDA_HOST_AND_DEVICE int getStagphaseOfLx(const int64_t& ivol,
					    const int& mu)
  {
    int ph=1;
    
    for(int nu=1;nu<=mu;nu++)
      ph*=(1-2*(glbCoordOfLoclx[ivol][nu-1]%2));
    
    return ph;
  }
}
