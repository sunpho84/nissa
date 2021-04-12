#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#define EXTERN_GEOMETRY_EO
 #include "geometry_eo.hpp"

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

//add a check that loc_vol is a multiple of 2
//#define REM_2 if(0)
#define REM_2

namespace nissa
{
  //compute the parity of a global site
  int glb_coord_parity(coords c)
  {
    int par=0;
    for(int mu=0;mu<NDIM;mu++) par+=c[mu];
    par%=2;
    
    return par;
  }
  int glblx_parity(int glx)
  {
    coords c;
    glb_coord_of_glblx(c,glx);
    
    return glb_coord_parity(c);
  }
  
  //set the eo geometry
  void set_eo_geometry()
  {
    if(not use_eo_geom) crash("E/O Geometry was not to be used!");
    if(eo_geom_inited) crash("E/O Geometry already initialized!");
    
    //check that all local sizes are multiples of 2
    int ok=1;
    REM_2 for(int mu=0;mu<NDIM;mu++) ok&=(locSize[mu]%2==0);
    REM_2 if(!ok) crash("local lattice size odd!");
    
    //set half the vol, bord and edge size
    glbVolh=(glbVol/2).nastyConvert();
    locVolh=(locVol/2).nastyConvert();
    
    //set the parity
    loclx_parity=nissa_malloc("loclx_parity",locVolWithBordAndEdge.nastyConvert(),int);
    ignore_borders_communications_warning(loclx_parity);
    
    loceo_of_loclx=nissa_malloc("loceo_of_loclx",locVolWithBordAndEdge.nastyConvert(),int);
    ignore_borders_communications_warning(loceo_of_loclx);
    
    for(int par=0;par<2;par++)
      {
	loclx_of_loceo[par]=nissa_malloc("loclx_of_loceo",(locVolh+bord_volh).nastyConvert()+edge_volh,int);
	loceo_neighup[par]=nissa_malloc("loceo_neighup",(locVolh+bord_volh).nastyConvert()+edge_volh,coords);
	loceo_neighdw[par]=nissa_malloc("loceo_neighdw",(locVolh+bord_volh).nastyConvert()+edge_volh,coords);
	surfeo_of_bordeo[par]=nissa_malloc("surfeo_of_bordeo",bord_volh,int);
	ignore_borders_communications_warning(loclx_of_loceo[par]);
	ignore_borders_communications_warning(loceo_neighup[par]);
	ignore_borders_communications_warning(loceo_neighdw[par]);
      }
    
    //Label the sites
    int iloc_eo[2]={0,0};
    for(int loclx=0;loclx<locVolWithBordAndEdge;loclx++)
      {
	//fix parity of local index
	int par=loclx_parity[loclx]=glb_coord_parity(glbCoordOfLoclx[loclx]);
	
	//associate the e/o index to lx sites and vice-versa
	loceo_of_loclx[loclx]=iloc_eo[par];
	loclx_of_loceo[par][iloc_eo[par]]=loclx;
	iloc_eo[par]++;
      }
    
    //Fix the movements among e/o ordered sites
    for(LocLxSite loclx=0;loclx<locVolWithBordAndEdge;loclx++)
      FOR_ALL_DIRECTIONS(mu)
	{
	  //take parity and e/o corresponding site
	  int par=loclx_parity[loclx.nastyConvert()];
	  int loceo=loceo_of_loclx[loclx.nastyConvert()];
	  
	  //up movements
	  const LocLxSite& loclxUp=loclxNeighup(loclx,mu);
	  if(loclxUp>=0 and loclxUp<locVolWithBordAndEdge)
	    loceo_neighup[par][loceo][mu.nastyConvert()]=loceo_of_loclx[loclxUp.nastyConvert()];
	  
	  //dw movements
	  const LocLxSite& loclxDw=loclxNeighdw(loclx,mu);
	  if(loclxDw>=0 and loclxDw<locVolWithBordAndEdge)
	    loceo_neighdw[par][loceo][mu.nastyConvert()]=loceo_of_loclx[loclxDw.nastyConvert()];
	}
    
    //finds how to fill the borders with surface
    for(BordLxSite bordlx=0;bordlx<bordVol;bordlx++)
      {
	const LocLxSite& surflx=loclxSiteAdjacentToBordLx(bordlx);
	const LocLxSite& loclx=extenedLocLxSiteOfBordLxSite(bordlx);
	surfeo_of_bordeo[loclx_parity[surflx.nastyConvert()]][loceo_of_loclx[loclx.nastyConvert()]-locVolh.nastyConvert()]=loceo_of_loclx[surflx.nastyConvert()];
      }
    
    master_printf("E/O Geometry intialized\n");
    
    eo_geom_inited=1;
  }
  
  //definitions of e/o split sender edges
  void initialize_eo_edge_senders_of_kind(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *base)
  {
    for(int par=0;par<2;par++)
      for(int vmu=0;vmu<2;vmu++)
	FOR_ALL_DIRECTIONS(mu)
	  for(int vnu=0;vnu<2;vnu++)
	    for(Direction nu=mu+1;nu<NDIM;nu++)
	      if(paral_dir[mu.nastyConvert()] and paral_dir[nu.nastyConvert()])
		{
		  int iedge=edge_numb[mu.nastyConvert()][nu.nastyConvert()];
		  int icomm=((par*2+vmu)*2+vnu)*NDIM*(NDIM-1)/2+iedge;
		  
		  //the sending edge might be a mess
		  const LocEoSite eo_edge_size=locVolh/locSize[mu.nastyConvert()]/locSize[nu.nastyConvert()];
		  int *edge_pos_disp=nissa_malloc("edge_disp",eo_edge_size.nastyConvert(),int);
		  int *single=nissa_malloc("single",eo_edge_size.nastyConvert(),int);
		  for(int iedge_eo=0;iedge_eo<eo_edge_size;iedge_eo++) single[iedge_eo]=1;
		  
		  int iedge_site=0;
		  for(int b_eo=0;b_eo<bord_volh;b_eo++)
		    {
		      const LocLxSite ivol=loclx_of_loceo[par][(locVolh+b_eo).nastyConvert()];
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
    for(int mu=0;mu<NDIM;mu++)
      for(int nu=mu+1;nu<NDIM;nu++)
	{
	  MPI_Type_contiguous(locVol.nastyConvert()/locSize[mu]/locSize[nu]/2,*base,&(MPI_EDGE_RECE[iedge]));
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
    
    for(int par=0;par<2;par++)
      {
	nissa_free(loclx_of_loceo[par]);
	nissa_free(surfeo_of_bordeo[par]);
	nissa_free(loceo_neighup[par]);
	nissa_free(loceo_neighdw[par]);
      }
    nissa_free(loclx_parity);
    nissa_free(loceo_of_loclx);
    
    eo_geom_inited=0;
  }
  
  //filter the points retaining only those having all even coord
  void filter_hypercube_origin_sites(color **vec)
  {
    NISSA_LOC_VOL_LOOP(ivol)
    {
      int save=1;
      for(int mu=0;mu<NDIM;mu++)
	save=save and glbCoordOfLoclx[ivol.nastyConvert()][mu]%2==0;
      
      if(!save)
	memset(vec[loclx_parity[ivol.nastyConvert()]][loceo_of_loclx[ivol.nastyConvert()]],0,sizeof(color));
    }
  }
}
