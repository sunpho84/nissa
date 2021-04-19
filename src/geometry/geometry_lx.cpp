
#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <math.h>
#include <string.h>

#define EXTERN_GEOMETRY_LX
# include <geometry/geometry_lx.hpp>

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
  //wrappers
  CUDA_HOST_DEVICE LocLxSite loclx_of_coord(const LocCoords& x)
  {
    return lx_of_coord<LocLxSite>(x,locSize);
  }
  
  //Return the coordinate of the rank containing the global coord
  void rank_coord_of_site_of_coord(RankCoords& rank_coord,const GlbCoords& glb_coord)
  {
    FOR_ALL_DIRECTIONS(mu)
      rank_coord(mu)=glb_coord(mu)()/locSize(mu)();
  }
  
  //Return the rank of passed coord
  Rank rank_of_coord(const RankCoords& x)
  {
    return lx_of_coord<Rank>(x,nrank_dir);
  }
  
  void coord_of_rank(RankCoords& c,const Rank& x)
  {
    coord_of_lx(c,x,nrank_dir);
  }
  
  //Return the rank containing the global coordinates
  Rank rank_hosting_site_of_coord(const GlbCoords& x)
  {
    RankCoords p;
    rank_coord_of_site_of_coord(p,x);
    
    return rank_of_coord(p);
  }
  
  //Return the rank containing the glblx passed
  Rank rank_hosting_glblx(const GlbLxSite& gx)
  {
    GlbCoords c;
    
    glb_coord_of_glblx(c,gx);
    
    return rank_hosting_site_of_coord(c);
  }
  
  //Return the local site and rank containing the global coordinates
  void get_loclx_and_rank_of_coord(LocLxSite& ivol,Rank& rank,const GlbCoords& g)
  {
    LocCoords l;
    RankCoords p;
    FOR_ALL_DIRECTIONS(mu)
      {
	p(mu)=g(mu)()/locSize(mu)();
	l(mu)=g(mu)()%locSize(mu)();
      }
    
    rank=rank_of_coord(p);
    ivol=loclx_of_coord(l);
  }
  
  GlbLxSite get_glblx_of_rank_and_loclx(const int irank,const LocLxSite& loclx)
  {
    RankCoords p;
    coord_of_rank(p,irank);
    crash("I believe it's wrong");
    GlbLxSite iglblx=0;
    FOR_ALL_DIRECTIONS(mu)
      iglblx=iglblx*glbSize(mu)()+locCoordOfLoclx(loclx,mu)();
    
    return iglblx;
  }
  
  /// Return the index of site of coord x in the 3d space obtained projecting away mu
  LocLxSite spatLxOfProjectedCoords(const LocCoords& x,const Direction& mu)
  {
    LocLxSite ilx=0;
    
    FOR_ALL_DIRECTIONS(nu)
      if(nu!=mu)
	ilx=ilx*locSize(nu)+x(nu);
    
    return ilx;
  }
  
  //return the index of the site of passed "pseudolocal" coordinate
  //if the coordinates are local, return the index according to the function loclx_of_coord
  //if exactly one of the coordinate is just out return its index according to bordlx_of_coord, incremented of previous border and loc_vol
  //if exactly two coordinates are outside, return its index according to edgelx_of_coord, incremented as before stated
  LocLxSite full_lx_of_coords(const LocCoords& ext_x)
  {
    //pseudo-localize it
    LocCoords x;
    FOR_ALL_DIRECTIONS(mu)
      {
	x(mu)=ext_x(mu);
	while(x(mu)<0) x(mu)+=glbSize(mu)();
	while(x(mu)>=glbSize(mu)()) x(mu)-=glbSize(mu)();
      }
    
    //check locality
    bool isloc=true;
    FOR_ALL_DIRECTIONS(mu)
      {
	isloc&=(x(mu)>=0);
	isloc&=(x(mu)<locSize(mu));
      }
    
    if(isloc)
      return loclx_of_coord(x);
    
    //check borderity
    Coords<int> is_bord;
    FOR_ALL_DIRECTIONS(mu)
      {
	is_bord(mu)=0;
	if(paral_dir(mu))
	  {
	    if(x(mu)==glbSize(mu)()-1) is_bord(mu)=-1;
	    if(x(mu)==locSize(mu)) is_bord(mu)=+1;
	  }
      }
    
    //check if it is in one of the NDIM forward or backward borders
    FOR_ALL_DIRECTIONS(mu)
      {
	int is=is_bord(mu);
	for(int inu=0;inu<NDIM-1;inu++)
	  is&=(is_bord(Direction(perp_dir[mu.nastyConvert()][inu]))==0);
	
	if(is)
	  {
	    /// Projected index
	    const auto iInsideDirBoord=+spatLxOfProjectedCoords(x,mu)();
	    const BordLxSite iInsideOrieBord=bord_offset(mu)+iInsideDirBoord;
	    
	    if(is_bord(mu)==-1)
	      return extenedLocLxSiteOfBordLxSite(iInsideOrieBord);             //backward border comes first
	    if(is_bord(mu)==+1)
	      return extenedLocLxSiteOfBordLxSite(bordVol()/2+iInsideOrieBord);  //forward border comes after
	    crash("if is bord should not arrive here");//nasty %d %d %d %d",ext_x[0],ext_x[1],ext_x[2],ext_x[3]);
	  }
      }
    
    //check if it is in one of the NDIM*(NDIM-1)/2 --,-+,+-,++ edges
    FOR_ALL_DIRECTIONS(mu)
      for(int inu=0;inu<NDIM-1;inu++)
	{
	  const Direction nu=perp_dir[mu.nastyConvert()][inu];
	  
	  //order mu,nu
	  const Direction al=(mu<nu)?mu:nu;
	  const Direction be=(mu>nu)?mu:nu;
	  
	  bool is=is_bord(mu) and is_bord(nu);
#if NDIM>=3
	  for(int irho=0;irho<NDIM-2;irho++)
	    is&=(is_bord(Direction(perp2_dir[mu.nastyConvert()][inu][irho]))==0);
#endif
	  
	  if(is)
	    {
	      int iedge=edge_numb[mu.nastyConvert()][nu.nastyConvert()];
	      if((is_bord(al)==-1) and (is_bord(be)==-1)) return (locVolWithBord+edge_offset[iedge].nastyConvert()+0*edge_vol.nastyConvert()/4+lineLxOfDoublyProjectedCoords(x,mu,nu)).nastyConvert();
	      if((is_bord(al)==-1) and (is_bord(be)==+1)) return (locVolWithBord+edge_offset[iedge].nastyConvert()+1*edge_vol.nastyConvert()/4+lineLxOfDoublyProjectedCoords(x,mu,nu)).nastyConvert();
	      if((is_bord(al)==+1) and (is_bord(be)==-1)) return (locVolWithBord+edge_offset[iedge].nastyConvert()+2*edge_vol.nastyConvert()/4+lineLxOfDoublyProjectedCoords(x,mu,nu)).nastyConvert();
	      if((is_bord(al)==+1) and (is_bord(be)==+1)) return (locVolWithBord+edge_offset[iedge].nastyConvert()+3*edge_vol.nastyConvert()/4+lineLxOfDoublyProjectedCoords(x,mu,nu)).nastyConvert();
	      crash("Edge: %d, mu=%d, nu=%d %d %d %d %d",iedge,mu(),nu(),ext_x(Direction(0))(),ext_x(xDirection)(),ext_x(yDirection)(),ext_x(zDirection)());
	    }
	}
    
    return -1;
  }
  
  /// Returns the border site adiacent at surface
  BordLxSite bordlx_of_surflx(const LocLxSite& loclx,const Direction& mu)
  {
    if(not paral_dir(mu))
      return -1;
    if(locSize(mu)<2)
      crash("not working if one dir is smaller than 2");
    
    if(locCoordOfLoclx(loclx,mu)==0)
      return bordLxSiteOfExtendedLocLxSize(loclxNeighdw(loclx,mu));
    if(locCoordOfLoclx(loclx,mu)==locSize(mu)-1)
      return bordLxSiteOfExtendedLocLxSize(loclxNeighup(loclx,mu));
    
    return -1;
  }
  
  //label all the sites: bulk, border and edge
  void label_all_sites()
  {
    //defined a box extending over borders/edges
    LocLxSite extended_box_vol=1;
    LocCoords extended_box_size;
    FOR_ALL_DIRECTIONS(mu)
      {
	extended_box_size(mu)=paral_dir(mu)*2+locSize(mu);
	extended_box_vol*=extended_box_size(mu)();
      }
    
    for(LocLxSite ivol=0;ivol<extended_box_vol;ivol++)
      {
	//subtract by one if dir is parallelized
	LocCoords x;
	coord_of_lx(x,ivol,extended_box_size);
	FOR_ALL_DIRECTIONS(mu)
	  if(paral_dir(mu))
	    x(mu)--;
	
	//check if it is defined
	const LocLxSite iloc=full_lx_of_coords(x);
	
	  if(iloc!=-1)
	  {
	    //compute global coordinates, assigning
	    FOR_ALL_DIRECTIONS(nu)
	      glbCoordOfLoclx(iloc,nu)=(x(nu)()+rank_coord(nu)()*locSize(nu)()+glbSize(nu)())%glbSize(nu)();
	    
	    /// Global index
	    GlbCoords g;
	    FOR_ALL_DIRECTIONS(mu)
	      g(mu)=glbCoordOfLoclx(iloc,mu);
	    const GlbLxSite iglb=glblx_of_coord(g);
	    
	    //if it is on the bulk store it
	    if(iloc<locVol)
	      {
		FOR_ALL_DIRECTIONS(nu)
		  locCoordOfLoclx(iloc,nu)=x(nu);
		glblxOfLoclx(iloc)=iglb;
	      }
	    
	    //if it is on the border store it
	    if(iloc>=locVol and iloc<locVolWithBord)
	      {
		const BordLxSite& ibord=bordLxSiteOfExtendedLocLxSize(iloc);
		glblxOfBordlx(ibord)=iglb;
		loclxOfBordlx(ibord)=iloc;
	      }
	    
	    //if it is on the edge store it
	    if(iloc>=locVolWithBord)
	      {
		const EdgeLxSite iedge=edgeLxSiteOfExtendedLocLxSize(iloc);
		glblxOfEdgelx(iedge)=iglb;
	      }
	  }
      }
  }
  
  //find the neighbours
  void find_neighbouring_sites()
  {
    //loop over the four directions
    for(LocLxSite ivol=0;ivol<locVolWithBordAndEdge;ivol++)
      FOR_ALL_DIRECTIONS(mu)
	{
	  //copy the coords
	  LocCoords n;
	  FOR_ALL_DIRECTIONS(nu)
	    n(nu)=glbCoordOfLoclx(ivol,nu)()-locSize(nu)()*rank_coord(nu)();
	  
	  //move forward
	  n(mu)++;
	  const LocLxSite nup=full_lx_of_coords(n); //nasty
	  
	  //move backward
	  n(mu)-=2;
	  const LocLxSite ndw=full_lx_of_coords(n);
	  
	  //if "local" assign it (automatically -1 otherwise)
	  loclxNeighup(ivol,mu)=nup;
	  loclxNeighdw(ivol,mu)=ndw;
	}
  }
  
  GlbLxSite glblxNeighup(const GlbLxSite& gx,const Direction& mu)
  {
    GlbCoords c;
    glb_coord_of_glblx(c,gx);
    c(mu)=(c(mu)+1)%glbSize(mu);
    
    return glblx_of_coord(c);
  }
  
  GlbLxSite glblxNeighdw(const GlbLxSite& gx,const Direction& mu)
  {
    GlbCoords c;
    glb_coord_of_glblx(c,gx);
    c(mu)=(c(mu)+glbSize(mu)-1)%glbSize(mu);
    
    return glblx_of_coord(c);
  }
  
  //finds how to fill the borders with opposite surface (up b->dw s)
  void find_surf_of_bord()
  {
    NISSA_LOC_VOL_LOOP(loclx)
      FOR_ALL_DIRECTIONS(mu)
	{
	  const BordLxSite bordlx=bordlx_of_surflx(loclx,mu);
	  if(bordlx!=-1) loclxSiteAdjacentToBordLx(bordlx)=loclx;
	}
  }
  
  //index all the sites on bulk
  void find_bulk_sites()
  {
    //check surfacity
    BulkLxSite ibulk=0;
    LocLxSite inon_fw_surf=0,inon_bw_surf=0;
    LocLxSite isurf=0,ifw_surf=0,ibw_surf=0;
    NISSA_LOC_VOL_LOOP(ivol)
      {
	//find if it is on bulk or non_fw or non_bw surf
	int is_bulk=true,is_non_fw_surf=true,is_non_bw_surf=true;
	FOR_ALL_DIRECTIONS(mu)
	  if(paral_dir(mu))
	    {
	      if(locCoordOfLoclx(ivol,mu)==locSize(mu)-1) is_bulk=is_non_fw_surf=false;
	      if(locCoordOfLoclx(ivol,mu)==0)              is_bulk=is_non_bw_surf=false;
	    }
	
	//mark it
	if(is_bulk) loclxOfBulklx(ibulk++)=ivol;
	if(is_non_fw_surf) loclxOfNonFwSurflx(inon_fw_surf++)=ivol;
	else               loclxOfFwSurflx(ifw_surf++)=ivol;
	if(is_non_bw_surf) loclxOfNonBwSurflx(inon_bw_surf++)=ivol;
	else               loclxOfBwSurflx(ibw_surf++)=ivol;
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
    FOR_ALL_DIRECTIONS(mu)
      MPI_Cart_shift(cart_comm,mu(),1,&(rank_neighdw(mu)()),&(rank_neighup(mu)()));
    
    FOR_ALL_DIRECTIONS(mu)
      {
	rank_neigh[0](mu)=rank_neighdw(mu);
	rank_neigh[1](mu)=rank_neighup(mu);
      }
    
    //borders
    glblxOfBordlx.allocate(bordVol);
    loclxOfBordlx.allocate(bordVol);
    loclxSiteAdjacentToBordLx.allocate(bordVol);
    
    //bulk and surfs
    loclxOfBulklx.allocate(bulkVol);
    loclxOfNonBwSurflx.allocate(nonBwSurfVol);
    loclxOfNonFwSurflx.allocate(nonFwSurfVol);
    loclxOfBwSurflx.allocate(bwSurfVol);
    loclxOfFwSurflx.allocate(fwSurfVol);
    
    //edges
    glblxOfEdgelx.allocate(edge_vol);
    
    locCoordOfLoclx.allocate(locVol);
    glbCoordOfLoclx.allocate(locVolWithBordAndEdge);
    loclxNeighdw.allocate(locVolWithBordAndEdge);
    loclxNeighup.allocate(locVolWithBordAndEdge);
    
    //local to global
    glblxOfLoclx.allocate(locVol);
    
    //borders
    glblxOfBordlx.allocate(bordVol);
    loclxOfBordlx.allocate(bordVol);
    loclxSiteAdjacentToBordLx.allocate(bordVol);
    
    //bulk and surfs
    loclxOfBulklx.allocate(bulkVol);
    loclxOfNonBwSurflx.allocate(nonBwSurfVol);
    loclxOfNonFwSurflx.allocate(nonFwSurfVol);
    loclxOfBwSurflx.allocate(bwSurfVol);
    loclxOfFwSurflx.allocate(fwSurfVol);
    
    //edges
    glblxOfEdgelx.allocate(edge_vol);
    
    //label the sites and neighbours
    label_all_sites();
    find_neighbouring_sites();
    
    //matches surface and opposite border
    find_surf_of_bord();
    
    //find bulk sites
    find_bulk_sites();
    
    //allocate a buffer large enough to allow communications of su3spinspin lx border
    recv_buf_size=std::max(recv_buf_size,(int64_t)(bordVol()*sizeof(su3spinspin)));
    send_buf_size=std::max(send_buf_size,(int64_t)(bordVol()*sizeof(su3spinspin)));
    
    //create the sweepers but do not fully initialize
    Wilson_sweeper=new gauge_sweeper_t;
    Symanzik_sweeper=new gauge_sweeper_t;
    
    //set locd geom (one of the dimension local and fastest running, the other as usual)
    max_locd_size=0;
    FOR_ALL_DIRECTIONS(mu)
      {//nasty
	remap_lx_to_locd[mu()]=remap_locd_to_lx[mu()]=NULL;
	max_locd_perp_size_per_dir(mu)=(glbVol()/glbSize(mu)()+nranks-1)/nranks;
	locd_perp_size_per_dir(mu)=(int)std::min((int64_t)max_locd_perp_size_per_dir(mu),glbVol()/glbSize(mu)()-max_locd_perp_size_per_dir(mu)*rank);
	verbosity_lv3_master_printf("rank %d locd_perp_size_per_dir[%d]: %d\n",rank,mu,locd_perp_size_per_dir(mu));
	locd_size_per_dir(mu)=locd_perp_size_per_dir(mu)*glbSize(mu).nastyConvert();
	max_locd_size=std::max(max_locd_size,locd_size_per_dir(mu));
      }
    
    master_printf("Cartesian geometry intialized\n");
  }
  
  //unset cartesian geometry
  void unset_lx_geometry()
  {
    if(lxGeomInited!=1) crash("cartesian geometry not initialized!");
    
    master_printf("Unsetting cartesian geometry\n");
    lxGeomInited=0;
    
    nissa_free(recv_buf);
    nissa_free(send_buf);
    
    delete Wilson_sweeper;
    delete Symanzik_sweeper;
  }
  
  //definitions of lexical ordered senders for edges
  void initialize_lx_edge_senders_of_kind(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *base)
  {
    const Direction T=0;
    const Direction X=1;
    const Direction Y=2;
    const Direction Z=3;
    
    //Various type useful for edges and sub-borders
    MPI_Datatype MPI_3_SLICE;
    MPI_Type_contiguous(locSize(Z)(),*base,&MPI_3_SLICE);
    
    ///////////define the sender for the 6 kinds of edges////////////
    //the 01 sender, that is simply a vector of L[Y] vector of L[Z] 
      MPI_Type_contiguous(locSize(Y)()*locSize(Z)(),*base,&(MPI_EDGE_SEND[0]));
      //the 0Y sender is a vector of L[1] segment of length L[Z] (already defined) separated by L[Y] of them
      MPI_Type_vector(locSize(X)(),1,locSize(Y)(),MPI_3_SLICE,&(MPI_EDGE_SEND[1]));
      //the 0Z sender is a vector of length L[1]xL[Y] of single elements, separated by L[Z] of them
      MPI_Type_vector(locSize(X)()*locSize(Y)(),1,locSize(Z)(),*base,&(MPI_EDGE_SEND[2]));
      //the 1Y sender should be equal to the 0Y sender, with 1->0
      MPI_Type_vector(locSize(T)(),1,locSize(Y)(),MPI_3_SLICE,&(MPI_EDGE_SEND[3]));
      //the 1Z sender should be equal to the 0Z sender, with 1->0
      MPI_Type_vector(locSize(T)()*locSize(Y)(),1,locSize(Z)(),*base,&(MPI_EDGE_SEND[4]));
      //the YZ sender should be equal to the 0Z sender with 1<->Y
      MPI_Type_vector(locSize(T)()*locSize(X)(),1,locSize(Z)(),*base,&(MPI_EDGE_SEND[5]));
      //Commit
      for(int iedge=0;iedge<6;iedge++) MPI_Type_commit(&(MPI_EDGE_SEND[iedge]));
  }
  
  //definitions of lexical ordered receivers for edges
  void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base)
  {
    //define the NDIM*(NDIM-1)/2 edges receivers, which are contiguous in memory
    int iedge=0;
    FOR_ALL_DIRECTIONS(mu)
      for(Direction nu=mu+1;nu<NDIM;nu++)
	{
	  MPI_Type_contiguous(locVol()/locSize(mu)()/locSize(nu)(),*base,&(MPI_EDGE_RECE[iedge]));
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
  
  // //define all the local lattice momenta
  // void define_local_momenta(momentum_t *k,double *k2,momentum_t *ktilde,double *ktilde2,momentum_t bc)
  // {
  //   if(!lxGeomInited) set_lx_geometry();
    
  //   //first of all, defines the local momenta for the various directions
  //   NISSA_LOC_VOL_LOOP(_imom)
  //     {
  // 	auto imom=_imom.nastyConvert();
  // 	k2[imom]=ktilde2[imom]=0;
  // 	FOR_ALL_DIRECTIONS(mu)
  // 	  {
  // 	    k[imom](mu)=M_PI*(2*glbCoordOfLoclx[imom](mu)+bc(mu))/glbSize(mu);
  // 	    ktilde[imom](mu)=sin(k[imom](mu));
	    
  // 	    k2[imom]+=k[imom](mu)*k[imom](mu);
  // 	    ktilde2[imom]+=ktilde[imom](mu)*ktilde[imom](mu);
  // 	  }
  //     }
  // }
  
  CUDA_HOST_DEVICE void get_stagphase_of_lx(Coords<int>& ph,const LocLxSite& ivol)
  {
    ph(timeDirection)=1;
    for(Direction mu=1;mu<NDIM;mu++)
      ph(mu)=ph(mu-1)*(1-2*(glbCoordOfLoclx(ivol,mu-1)()%2));
  }
  
  CUDA_HOST_DEVICE int get_stagphase_of_lx(const LocLxSite& ivol,const Direction& mu)
  {
    int ph=1;
    for(Direction nu=1;nu<=mu;nu++)
      ph*=(1-2*(glbCoordOfLoclx(ivol,nu-1)%2))();
    return ph;
  }
}
