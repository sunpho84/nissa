#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>

#include "communicate.hpp"

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_lx.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

/* and this is the shape of the edges
_______________________________________________________________________________________________________________________
|______________t-x-_edge______________|||______________t-y-_edge______________|||______________t-z-_edge______________|
|____y=0____||____y=1____||____y=2____|||____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|
| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z ||| y | y | y || y | y | y || y | y | y |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
-----------------------------------------------------------------------------------------------------------------------

_______________________________________________________________________________________________________________________
|______________x-y-_edge______________|||______________x-z-_edge______________|||______________y-z-_edge______________|
|____t=0____||____t=1____||____t=2____|||____t=0____||____t=1____||____t=2____|||____t=0____||____t=1____||____t=2____|
| z | z | z || z | z | z || z | z | z ||| y | y | z || y | y | y || y | y | y ||| x | x | x || x | x | x || x | x | x |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
-----------------------------------------------------------------------------------------------------------------------

then follows all the i-j+ edges, the i+j- and after the i+j+
*/

namespace nissa
{
  //Send the edges of lx vector
  void communicate_lx_edges(char *data,comm_t &bord_comm,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site)
  {
    int min_size=nbytes_per_site*(edge_vol+bord_vol+loc_vol);
    crash_if_edges_not_allocated(data,min_size);
    communicate_lx_borders(data,bord_comm);
    
    if(!check_edges_valid(data))
      {
	GET_THREAD_ID();
	
	if(IS_MASTER_THREAD)
	  {
	    int nrequest=0;
	    MPI_Request request[NDIM*(NDIM-1)*4];
	    MPI_Status status[NDIM*(NDIM-1)*4];
	    int send,rece;
	    int imessage=0;
	    coords x;
	    memset(x,0,sizeof(coords));
	    
	    for(int idir=0;idir<NDIM;idir++)
	      for(int jdir=idir+1;jdir<NDIM;jdir++)
		if(paral_dir[idir] and paral_dir[jdir])
		  {
		    int iedge=edge_numb[idir][jdir];
		    int pos_edge_offset;
		    
		    //take the starting point of the border
		    x[jdir]=loc_size[jdir]-1;
		    pos_edge_offset=bordlx_of_coord(x,idir);
		    x[jdir]=0;
		    
		    //Send the i-j- internal edge to the j- rank as i-j+ external edge
		    send=(loc_vol+bord_offset[idir])*nbytes_per_site;
		    rece=(loc_vol+bord_vol+edge_offset[iedge]+edge_vol/4)*nbytes_per_site;
		    MPI_Irecv(data+rece,1,MPI_EDGES_RECE[iedge],rank_neighup[jdir],imessage,cart_comm,request+(nrequest++));
		    MPI_Isend(data+send,1,MPI_EDGES_SEND[iedge],rank_neighdw[jdir],imessage++,cart_comm,request+(nrequest++));
		    
		    //Send the i-j+ internal edge to the j+ rank as i-j- external edge
		    send=(loc_vol+bord_offset[idir]+pos_edge_offset)*nbytes_per_site;
		    rece=(loc_vol+bord_vol+edge_offset[iedge])*nbytes_per_site;
		    MPI_Irecv(data+rece,1,MPI_EDGES_RECE[iedge],rank_neighdw[jdir],imessage,cart_comm,request+(nrequest++));
		    MPI_Isend(data+send,1,MPI_EDGES_SEND[iedge],rank_neighup[jdir],imessage++,cart_comm,request+(nrequest++));
		    
		    //Send the i+j- internal edge to the j- rank as i+j+ external edge
		    send=(loc_vol+bord_offset[idir]+bord_vol/2)*nbytes_per_site;
		    rece=(loc_vol+bord_vol+edge_offset[iedge]+3*edge_vol/4)*nbytes_per_site;
		    MPI_Irecv(data+rece,1,MPI_EDGES_RECE[iedge],rank_neighup[jdir],imessage,cart_comm,request+(nrequest++));
		    MPI_Isend(data+send,1,MPI_EDGES_SEND[iedge],rank_neighdw[jdir],imessage++,cart_comm,request+(nrequest++));
		    
		    //Send the i+j+ internal edge to the j+ rank as i+j- external edge
		    send=(loc_vol+bord_offset[idir]+bord_vol/2+pos_edge_offset)*nbytes_per_site;
		    rece=(loc_vol+bord_vol+edge_offset[iedge]+edge_vol/2)*nbytes_per_site;
		    MPI_Irecv(data+rece,1,MPI_EDGES_RECE[iedge],rank_neighdw[jdir],imessage,cart_comm,request+(nrequest++));
		    MPI_Isend(data+send,1,MPI_EDGES_SEND[iedge],rank_neighup[jdir],imessage++,cart_comm,request+(nrequest++));
		    imessage++;
		  }
	    
	    if(nrequest>0) MPI_Waitall(nrequest,request,status);
	  } //end IS_MASTER_THREAD
	set_edges_valid(data);
      }
  }
  
  //Send the edges u su3: usefuls for hyp
  void communicate_lx_su3_edges(su3 *u)
  {communicate_lx_edges((char*)u,lx_su3_comm,MPI_LX_SU3_EDGES_SEND,MPI_LX_SU3_EDGES_RECE,sizeof(su3));}
  
  //Send the edges of the gauge configuration
  void communicate_lx_quad_su3_edges(quad_su3 *conf)
  {communicate_lx_edges((char*)conf,lx_quad_su3_comm,MPI_LX_QUAD_SU3_EDGES_SEND,MPI_LX_QUAD_SU3_EDGES_RECE,sizeof(quad_su3));}
  
  //Send the edges of the as2t_su3 field (for clover derivative)
  void communicate_lx_as2t_su3_edges(as2t_su3 *a)
  {communicate_lx_edges((char*)a,lx_as2t_su3_comm,MPI_LX_AS2T_SU3_EDGES_SEND,MPI_LX_AS2T_SU3_EDGES_RECE,sizeof(as2t_su3));}
  
  ///////////////////////////////////////////// e/o geometry /////////////////////////////////////////
  
  //Send the edges of eo vector
  void communicate_eo_edges(char **data,comm_t &bord_comm,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site)
  {
    if(!check_edges_valid(data[EVN])||!check_edges_valid(data[ODD]))
      {
	//first make sure that borders are communicated
	communicate_ev_and_od_borders((void**)data,bord_comm);
	
	GET_THREAD_ID();
	
	if(IS_MASTER_THREAD)
	  {
	    int nrequest=0;
	    int nrequest_tot=0;
	    for(int mu=0;mu<NDIM;mu++)
	      for(int nu=mu+1;nu<NDIM;nu++)
		if(paral_dir[mu] && paral_dir[nu])
		  nrequest_tot+=16;
	    
	    MPI_Request request[nrequest_tot];
	    MPI_Status status[nrequest_tot];
	    int imessage=0;
	    
	    verbosity_lv3_master_printf("Communicating edges of %s\n",get_vect_name(data[EVN]));
	    
	    for(int par=0;par<2;par++)
	      {
		//check_edges_allocated(data);
		int min_size=nbytes_per_site*(edge_volh+bord_volh+loc_volh);
		crash_if_edges_not_allocated(data[par],min_size);
		
		//"v" refer to the verse of the dir
		for(int vmu=0;vmu<2;vmu++)
		  for(int mu=0;mu<NDIM;mu++)
		    for(int vnu=0;vnu<2;vnu++)
		      for(int nu=mu+1;nu<NDIM;nu++)
			if(paral_dir[mu] && paral_dir[nu])
			  {
			    int iedge=edge_numb[mu][nu];
			    
			    //communicators verse refers to the internal edge
			    int icomm_send=((par*2+vmu)*2+vnu)*6+iedge;
			    int icomm_recv=iedge;
			    
			    //Send the (mu,nu) internal edge to the nu rank as (mu,-nu) external edge
			    //Receive the (mu,-nu) external edge from the -nu rank (mu,nu) internal edge
			    int ext_edge_recv_start=(loc_volh+bord_volh+edge_volh/4*(vmu*2+!vnu)+edge_offset[iedge]/2)*nbytes_per_site;
			    int int_edge_send_start=loc_volh*nbytes_per_site;
			    MPI_Irecv(data[par]+ext_edge_recv_start,1,MPI_EDGES_RECE[icomm_recv],rank_neigh[!vnu][nu],imessage,cart_comm,request+(nrequest++));
			    MPI_Isend(data[par]+int_edge_send_start,1,MPI_EDGES_SEND[icomm_send],rank_neigh[vnu][nu],imessage,cart_comm,request+(nrequest++));
			    
			    imessage++;
			  }
	      }
	    
	    if(nrequest!=nrequest_tot) crash("something went wrong");
	    
	    if(nrequest>0) MPI_Waitall(nrequest,request,status);
	  }
	
	for(int par=0;par<2;par++) set_edges_valid(data[par]);
      }
  }
  
  //Send the edges of the gauge configuration
  void communicate_eo_quad_su3_edges(quad_su3 **conf)
  {communicate_eo_edges((char**)conf,lx_quad_su3_comm,MPI_EO_QUAD_SU3_EDGES_SEND,MPI_EO_QUAD_SU3_EDGES_RECE,sizeof(quad_su3));}
  
  void communicate_eo_as2t_su3_edges(as2t_su3 **a)
  {communicate_lx_edges((char*)a,lx_as2t_su3_comm,MPI_EO_AS2T_SU3_EDGES_SEND,MPI_EO_AS2T_SU3_EDGES_RECE,sizeof(as2t_su3));}
}
