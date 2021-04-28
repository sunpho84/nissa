#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <string.h>

#define EXTERN_GEOMETRY_EO
# include "geometry_eo.hpp"

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //compute the parity of a global site
  Parity glb_coord_parity(const GlbCoords& c)
  {
    GlbCoord sum=0;
    FOR_ALL_DIRS(mu)
      sum+=c(mu);
    
    return static_cast<int>(sum()%2);
  }
  
  Parity glblx_parity(const GlbLxSite& glx)
  {
    GlbCoords c;
    glb_coord_of_glblx(c,glx);
    
    return glb_coord_parity(c);
  }
  
  //set the eo geometry
  void set_eo_geometry()
  {
    if(not use_eo_geom)
      crash("E/O Geometry was not to be used!");
    if(eo_geom_inited)
      crash("E/O Geometry already initialized!");
    
    //check that all local sizes are multiples of 2
    bool ok=true;
   FOR_ALL_DIRS(mu)
     ok&=(locSize(mu)%2==0);
    
    if(not ok)
      crash("local lattice size odd!");
    
    //set half the vol, bord and edge size
    glbVolh=glbVol()/2;
    locVolh=locVol()/2;
    
    //set the parity
    loclx_parity.allocate(locVolWithBordAndEdge);
    
    loceo_of_loclx.allocate(locVolWithBordAndEdge);
    
    for(int par=0;par<2;par++)
      {
	loclx_of_loceo.allocate(locVolhWithBordAndEdge);
	loceo_neighup.allocate(locVolhWithBordAndEdge);
	loceo_neighdw.allocate(locVolhWithBordAndEdge);
	surfeo_of_bordeo.allocate(bordVolh);
      }
    
    //Label the sites
    Tensor<OfComps<Parity>,LocEoSite> iloc_eo;
    FOR_BOTH_PARITIES(par)
      iloc_eo(par)=0;
    
    for(LocLxSite loclx=0;loclx<locVolWithBordAndEdge;loclx++)
      {
	//fix parity of local index
	const Parity par=
	  loclx_parity(loclx)=glblx_parity(glblxOfLoclx(loclx));
	
	//associate the e/o index to lx sites and vice-versa
	loceo_of_loclx(loclx)=iloc_eo(par);
	loclx_of_loceo(par,iloc_eo(par))=loclx;
	iloc_eo(par)++;
      }
    
    //Fix the movements among e/o ordered sites
    for(LocLxSite loclx=0;loclx<locVolWithBordAndEdge;loclx++)
      FOR_ALL_DIRS(mu)
	{
	  //take parity and e/o corresponding site
	  const Parity& par=loclx_parity(loclx);
	  const LocEoSite& loceo=loceo_of_loclx(loclx);
	  
	  //up movements
	  const LocLxSite& loclxUp=loclxNeighup(loclx,mu);
	  if(loclxUp>=0 and loclxUp<locVolWithBordAndEdge)
	    loceo_neighup(par,loceo,mu)=loceo_of_loclx(loclxUp);
	  
	  //dw movements
	  const LocLxSite& loclxDw=loclxNeighdw(loclx,mu);
	  if(loclxDw>=0 and loclxDw<locVolWithBordAndEdge)
	    loceo_neighdw(par,loceo,mu)=loceo_of_loclx(loclxDw);
	}
    
    //finds how to fill the borders with surface
    for(BordLxSite bordlx=0;bordlx<bordVol;bordlx++)
      {
	const LocLxSite& surflx=loclxSiteAdjacentToBordLx(bordlx);
	const LocLxSite& loclx=extenedLocLxSiteOfBordLxSite(bordlx);
	surfeo_of_bordeo(loclx_parity(surflx),bordEoSiteOfExtendedLocEoSize(loceo_of_loclx(loclx)))=loceo_of_loclx(surflx);
      }
    
    master_printf("E/O Geometry intialized\n");
    
    loclx_parity.consolidate();
    loceo_of_loclx.consolidate();
    loclx_of_loceo.consolidate();
    surfeo_of_bordeo.consolidate();
    loceo_neighup.consolidate();
    loceo_neighdw.consolidate();
    
    eo_geom_inited=1;
  }
  
  //definitions of e/o split sender edges
  void initialize_eo_edge_senders_of_kind(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *base)
  {
    FOR_BOTH_PARITIES(par)
      for(int vmu=0;vmu<2;vmu++)
	FOR_ALL_DIRS(mu)
	  for(int vnu=0;vnu<2;vnu++)
	    for(Dir nu=mu+1;nu<NDIM;nu++)
	      if(paral_dir(mu) and paral_dir(nu))
		{
		  int iedge=edge_numb[mu.nastyConvert()][nu.nastyConvert()];
		  int icomm=((par()*2+vmu)*2+vnu)*NDIM*(NDIM-1)/2+iedge;
		  
		  //the sending edge might be a mess
		  const LocEoSite eo_edge_size=locVolh/locSize(mu)()/locSize(nu)();
		  int *edge_pos_disp=nissa_malloc("edge_disp",eo_edge_size.nastyConvert(),int);
		  int *single=nissa_malloc("single",eo_edge_size.nastyConvert(),int);
		  for(int iedge_eo=0;iedge_eo<eo_edge_size;iedge_eo++) single[iedge_eo]=1;
		  
		  int iedge_site=0;
		  for(int b_eo=0;b_eo<bordVolh;b_eo++)
		    {
		      const LocLxSite ivol=loclx_of_loceo(par,extenedLocEoSiteOfBordEoSite(b_eo));
		      if(loclxNeigh(!vmu)(ivol,mu)>=0 and loclxNeigh(!vmu)(ivol,mu)<locVol and loclxNeigh(vnu)(ivol,nu)>=locVolWithBord)
			edge_pos_disp[iedge_site++]=b_eo;
		    }
		  if(iedge_site!=eo_edge_size)
		    crash("iedge_site=%d did not arrive to eo_edge_size=%d",iedge_site,eo_edge_size);
		  
		  MPI_Type_indexed(eo_edge_size(),single,edge_pos_disp,*base,&(MPI_EO_EDGES_SEND[icomm]));
		  //commit the mess
		  MPI_Type_commit(&(MPI_EO_EDGES_SEND[icomm]));
		  
		  nissa_free(single);
		  nissa_free(edge_pos_disp);
		}
  }
  
  //definitions of e/o split receivers for edges
  void initialize_eo_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base)
  {
    //define the NDIM*(NDIM-1)/2 edges receivers, which are contiguous in memory
    int iedge=0;
   FOR_ALL_DIRS(mu)
      for(Dir nu=mu+1;nu<NDIM;nu++)
	{
	  MPI_Type_contiguous(locVol()/locSize(mu)()/locSize(nu)()/2,*base,&(MPI_EDGE_RECE[iedge]));
	  MPI_Type_commit(&(MPI_EDGE_RECE[iedge]));
	  iedge++;
	}
  }
  
  //initalize senders and receivers for edges of lexically ordered vectors
  void set_eo_edge_senders_and_receivers(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *MPI_EO_EDGES_RECE,MPI_Datatype *base)
  {
    initialize_eo_edge_senders_of_kind(MPI_EO_EDGES_SEND,base);
    initialize_eo_edge_receivers_of_kind(MPI_EO_EDGES_RECE,base);
  }
  
  //unset the eo geometry
  void unset_eo_geometry()
  {
    if(not eo_geom_inited)
      crash("asking to unset never initialized E/O Geometry!");
    
    master_printf("Unsetting E/O Geometry\n");
    
    loclx_parity.dealloc();
    loceo_of_loclx.dealloc();
    loclx_of_loceo.dealloc();
    surfeo_of_bordeo.dealloc();
    loceo_neighup.dealloc();
    loceo_neighdw.dealloc();
    
    eo_geom_inited=0;
  }
  
  //filter the points retaining only those having all even coord
  void filter_hypercube_origin_sites(color **vec)
  {
    NISSA_LOC_VOL_LOOP(ivol)
    {
      int save=1;
      FOR_ALL_DIRS(mu)
	save=save and glbCoordOfLoclx(ivol,mu)%2==0;
      
      if(not save)
	memset(vec[loclx_parity(ivol).nastyConvert()][loceo_of_loclx(ivol).nastyConvert()],0,sizeof(color));
    }
  }
}
