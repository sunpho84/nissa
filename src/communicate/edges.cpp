#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <mpi.h>

#include "communicate.h"

#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/thread_macros.h"
#include "../base/vectors.h"
#include "../geometry/geometry_lx.h"
#include "../new_types/new_types_definitions.h"
#include "../routines/ios.h"
#ifdef USE_THREADS
 #include "../routines/thread.h"
#endif

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

//Send the edges of lx vector
void communicate_lx_edges(char *data,comm_t &bord_comm,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site)
{
  crash_if_edges_not_allocated(data);
  communicate_lx_borders(data,bord_comm);
  
  if(!check_edges_valid(data))
    {
      GET_THREAD_ID();
      
      if(IS_MASTER_THREAD)
	{
	  int nrequest=0;
	  MPI_Request request[48];
	  MPI_Status status[48];
	  int send,rece;
	  int imessage=4532543;
	  int x[4]={0,0,0,0};
	  
	  for(int idir=0;idir<4;idir++)
	    for(int jdir=idir+1;jdir<4;jdir++)
	      if(paral_dir[idir] && paral_dir[jdir])
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
	  for(int mu=0;mu<4;mu++)
	    for(int nu=mu+1;nu<4;nu++)
	      if(paral_dir[mu] && paral_dir[nu])
		nrequest_tot+=16;
	  
	  MPI_Request request[nrequest_tot];
	  MPI_Status status[nrequest_tot];
	  int imessage=8534240;
	  
	  verbosity_lv3_master_printf("Communicating edges of %s\n",get_vec_name(data[EVN]));
	  
	  //check_edges_allocated(data);
	  for(int par=0;par<2;par++)
	    {
	      crash_if_edges_not_allocated(data[par]);
	      
	      //"v" refer to the verse of the dir
	      for(int vmu=0;vmu<2;vmu++)
		for(int mu=0;mu<4;mu++)
		  for(int vnu=0;vnu<2;vnu++)
		    for(int nu=mu+1;nu<4;nu++)
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
