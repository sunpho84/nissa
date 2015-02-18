#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/global_variables.hpp"
#include "base/debug.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

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
    if(!use_eo_geom) crash("E/O Geometry was not to be used!");
    if(eo_geom_inited) crash("E/O Geometry already initialized!");
    
    //check that all local sizes are multiples of 2
    int ok=1;
    for(int mu=0;mu<NDIM;mu++) ok&=(loc_size[mu]%2==0);
    if(!ok) crash("local lattice size odd!");
    
    //set half the vol, bord and edge size
    glb_volh=glb_vol/2;
    loc_volh=loc_vol/2;
    
    //set the parity
    loclx_parity=nissa_malloc("loclx_parity",loc_vol+bord_vol+edge_vol,int);
    ignore_borders_communications_warning(loclx_parity);
    
    loceo_of_loclx=nissa_malloc("loceo_of_loclx",loc_vol+bord_vol+edge_vol,int);
    ignore_borders_communications_warning(loceo_of_loclx);
    
    for(int par=0;par<2;par++) loclx_of_loceo[par]=nissa_malloc("loclx_of_loceo",loc_volh+bord_volh+edge_volh,int);
    for(int par=0;par<2;par++) loceo_neighup[par]=nissa_malloc("loceo_neighup",loc_volh+bord_volh+edge_volh,coords);
    for(int par=0;par<2;par++) loceo_neighdw[par]=nissa_malloc("loceo_neighdw",loc_volh+bord_volh+edge_volh,coords);
    for(int par=0;par<2;par++) surfeo_of_bordeo[par]=nissa_malloc("surfeo_of_bordeo",bord_volh,int);
    for(int par=0;par<2;par++) ignore_borders_communications_warning(loclx_of_loceo[par]);
    for(int par=0;par<2;par++) ignore_borders_communications_warning(loceo_neighup[par]);
    for(int par=0;par<2;par++) ignore_borders_communications_warning(loceo_neighdw[par]);
    
    //Label the sites
    int iloc_eo[2]={0,0};
    for(int loclx=0;loclx<loc_vol+bord_vol+edge_vol;loclx++)
      {
	//fix parity of local index
	int par=loclx_parity[loclx]=glb_coord_parity(glb_coord_of_loclx[loclx]);
	
	//associate the e/o index to lx sites and vice-versa
	loceo_of_loclx[loclx]=iloc_eo[par];
	loclx_of_loceo[par][iloc_eo[par]]=loclx;
	iloc_eo[par]++;
      }
    
    //Fix the movements among e/o ordered sites
    for(int loclx=0;loclx<loc_vol+bord_vol+edge_vol;loclx++)
      for(int mu=0;mu<NDIM;mu++)
	{
	  //take parity and e/o corresponding site
	  int par=loclx_parity[loclx];
	  int loceo=loceo_of_loclx[loclx];
	  
	  //up movements
	  int loclx_up=loclx_neighup[loclx][mu];
	  if(loclx_up>=0 && loclx_up<loc_vol+bord_vol+edge_vol)
	    loceo_neighup[par][loceo][mu]=loceo_of_loclx[loclx_up];
	  
	  //dw movements
	  int loclx_dw=loclx_neighdw[loclx][mu];
	  if(loclx_dw>=0 && loclx_dw<loc_vol+bord_vol+edge_vol)
	    loceo_neighdw[par][loceo][mu]=loceo_of_loclx[loclx_dw];
	}
    
    //finds how to fill the borders with surface
    for(int bordlx=0;bordlx<bord_vol;bordlx++)
      {
	int surflx=surflx_of_bordlx[bordlx];
	surfeo_of_bordeo[loclx_parity[surflx]][loceo_of_loclx[bordlx+loc_vol]-loc_volh]=loceo_of_loclx[surflx];
      }

    master_printf("E/O Geometry intialized\n");
    
    eo_geom_inited=1;
  }
  
  //definitions of e/o split sender edges
  void initialize_eo_edge_senders_of_kind(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *base)
  {
    for(int par=0;par<2;par++)
      for(int vmu=0;vmu<2;vmu++)
	for(int mu=0;mu<NDIM;mu++)
	  for(int vnu=0;vnu<2;vnu++)
	    for(int nu=mu+1;nu<NDIM;nu++)
	      if(paral_dir[mu] && paral_dir[nu])
		{
		  int iedge=edge_numb[mu][nu];
		  int icomm=((par*2+vmu)*2+vnu)*NDIM*(NDIM-1)/2+iedge;
		  
		  //the sending edge might be a mess
		  int eo_edge_size=loc_volh/loc_size[mu]/loc_size[nu];
		  int *edge_pos_disp=nissa_malloc("edge_disp",eo_edge_size,int);
		  int *single=nissa_malloc("single",eo_edge_size,int);
		  for(int iedge_eo=0;iedge_eo<eo_edge_size;iedge_eo++) single[iedge_eo]=1;
		  
		  int iedge_site=0;
		  for(int b_eo=0;b_eo<bord_volh;b_eo++)
		    {
		      int ivol=loclx_of_loceo[par][loc_volh+b_eo];
		      if(loclx_neigh[!vmu][ivol][mu]>=0 && loclx_neigh[!vmu][ivol][mu]<loc_vol && loclx_neigh[vnu][ivol][nu]>=loc_vol+bord_vol) edge_pos_disp[iedge_site++]=b_eo;
		    }
		  if(iedge_site!=eo_edge_size) crash("iedge_site=%d did not arrive to eo_edge_size=%d",iedge_site,eo_edge_size);
		  
		  MPI_Type_indexed(eo_edge_size,single,edge_pos_disp,*base,&(MPI_EO_EDGES_SEND[icomm]));
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
	  MPI_Type_contiguous(loc_vol/loc_size[mu]/loc_size[nu]/2,*base,&(MPI_EDGE_RECE[iedge]));
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
    if(!eo_geom_inited)
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
  
  //add or remove staggered phases to/from a conf
  THREADABLE_FUNCTION_1ARG(addrem_stagphases_to_eo_conf, quad_su3**,eo_conf)
  {
    //we must ensure that nobody is using the conf
    THREAD_BARRIER();
    
    //work also on borders and edges if allocated and valid
    int ending=loc_volh;
    if(check_borders_allocated(eo_conf[0]) && check_borders_allocated(eo_conf[1]) && check_borders_valid(eo_conf[0]) && check_borders_valid(eo_conf[1])) ending+=bord_volh;
    if(check_edges_allocated(eo_conf[0])   && check_edges_allocated(eo_conf[1])   && check_edges_valid(eo_conf[0])   && check_edges_valid(eo_conf[1])) ending+=edge_volh;
    
    GET_THREAD_ID();
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol_eo,0,ending)
	  {
	    int ivol_lx=loclx_of_loceo[par][ivol_eo];
	    
	    int d=0;
	    
	    //phase in direction 1 is always 0 so nothing has to be done in that dir
	    //if(d%2==1) su3_prod_double(eo_conf[par][ivol_eo][1],eo_conf[par][ivol_eo][1],-1);
	    
	    //direction 2
	    d+=glb_coord_of_loclx[ivol_lx][1];
	    if(d%2==1) su3_prod_double(eo_conf[par][ivol_eo][2],eo_conf[par][ivol_eo][2],-1);
	    
	    //direction 3
	    d+=glb_coord_of_loclx[ivol_lx][2];
	    if(d%2==1) su3_prod_double(eo_conf[par][ivol_eo][3],eo_conf[par][ivol_eo][3],-1);
	    
	    //direction 0
	    d+=glb_coord_of_loclx[ivol_lx][3];
	    //debug: putting the anti-periodic condition on the temporal border
	    //in future remove it!!!
	    if(glb_coord_of_loclx[ivol_lx][0]==glb_size[0]-1) d+=1;
	    if(d%2==1) su3_prod_double(eo_conf[par][ivol_eo][0],eo_conf[par][ivol_eo][0],-1);
	  }
      }
    
    //now the conf is again available
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END

  //filter the points retaining only those having all even coord
  void filter_hypercube_origin_sites(color **vec)
  {
    NISSA_LOC_VOL_LOOP(ivol)
    {
      int save=1;
      for(int mu=0;mu<NDIM;mu++)
	save=save and glb_coord_of_loclx[ivol][mu]%2==0;
      
      if(!save)
	memset(vec[loclx_parity[ivol]][loceo_of_loclx[ivol]],0,sizeof(color));
    }
  }
}
