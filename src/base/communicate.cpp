#include <mpi.h>

#include "debug.h"
#include "../geometry/geometry_lx.h"
#include "../new_types/new_types_definitions.h"
#include "global_variables.h"
#include "routines.h"
#include "vectors.h"

/* This is the shape and ordering of the border in the memory, for a 3^4 lattice
_______________________________________________________________________________________________________________________
|___________________________________________________dir____=_____0____________________________________________________|
|_______________x__=__0_______________|||_______________x__=__1_______________|||_______________x__=__2_______________|
|____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|
| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
-----------------------------------------------------------------------------------------------------------------------

_______________________________________________________________________________________________________________________
|___________________________________________________dir____=_____1____________________________________________________|
|_______________t__=__0_______________|||_______________t__=__1_______________|||_______________t__=__2_______________|
|____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|||____y=0____||____y=1____||____y=2____|
| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
-----------------------------------------------------------------------------------------------------------------------

_______________________________________________________________________________________________________________________
|___________________________________________________dir____=_____2____________________________________________________|
|_______________t__=__0_______________|||_______________t__=__1_______________|||_______________t__=__2_______________|
|____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|
| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z ||| z | z | z || z | z | z || z | z | z |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
-----------------------------------------------------------------------------------------------------------------------

_______________________________________________________________________________________________________________________
|___________________________________________________dir____=_____3____________________________________________________|
|_______________t__=__0_______________|||_______________t__=__1_______________|||_______________t__=__2_______________|
|____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|||____x=0____||____x=1____||____x=2____|
| y | y | y || y | y | y || y | y | y ||| y | y | y || y | y | y || y | y | y ||| y | y | y || y | y | y || y | y | y |
|012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|||012|012|012||012|012|012||012|012|012|
-----------------------------------------------------------------------------------------------------------------------

and following, the positive (forward sending) direction borders
*/

//Send the borders of the data
void start_communicating_lx_borders(int &nrequest,MPI_Request *request,char *data,MPI_Datatype *MPI_BORDS_SEND,MPI_Datatype *MPI_BORDS_RECE,int nbytes_per_site)
{
  int imessage=654325;
  nrequest=0;
  
  crash_if_borders_not_allocated(data);  
  
  if(!check_borders_valid(data))
    {
      verbosity_lv3_master_printf("Start communicating borders of %s\n",get_vec_name((void*)data));
      for(int mu=0;mu<4;mu++)
	if(paral_dir[mu]!=0)
	  {
	    //sending the upper border to the lower node
	    MPI_Irecv(data+start_lx_bord_rece_up[mu]*nbytes_per_site,1,MPI_BORDS_RECE[mu],rank_neighup[mu],imessage,cart_comm,request+(nrequest++));
	    MPI_Isend(data+start_lx_bord_send_up[mu]*nbytes_per_site,1,MPI_BORDS_SEND[mu],rank_neighdw[mu],imessage++,cart_comm,request+(nrequest++));
	    
	    //sending the lower border to the upper node
	    MPI_Irecv(data+start_lx_bord_rece_dw[mu]*nbytes_per_site,1,MPI_BORDS_RECE[mu],rank_neighdw[mu],imessage,cart_comm,request+(nrequest++));
	    MPI_Isend(data+start_lx_bord_send_dw[mu]*nbytes_per_site,1,MPI_BORDS_SEND[mu],rank_neighup[mu],imessage++,cart_comm,request+(nrequest++));
	  }
    }
}

//wait to finish communications
void finish_communicating_lx_borders(int &nrequest,MPI_Request *request,char *data)
{
  tot_nissa_comm_time-=take_time();
  if(nrequest>0)
    {
      verbosity_lv3_master_printf("Waiting to finish %d communication of borders of vector %s\n",nrequest,get_vec_name(data));
      MPI_Status status[nrequest];
      MPI_Waitall(nrequest,request,status);
      nrequest=0;
    }
  set_borders_valid(data);
  set_edges_invalid(data);
  
  tot_nissa_comm_time+=take_time();
}

void communicate_lx_borders(char *data,MPI_Datatype *MPI_BORDS_SEND,MPI_Datatype *MPI_BORDS_RECE,int nbytes_per_site)
{
  MPI_Request request[16];
  int nrequest;
  start_communicating_lx_borders(nrequest,request,data,MPI_BORDS_SEND,MPI_BORDS_RECE,nbytes_per_site);
  finish_communicating_lx_borders(nrequest,request,data);
}

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
void communicate_lx_edges(char *data,MPI_Datatype *MPI_BORDS_SEND,MPI_Datatype *MPI_BORDS_RECE,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site)
{
  check_edges_allocated(data);
  communicate_lx_borders(data,MPI_BORDS_SEND,MPI_BORDS_RECE,nbytes_per_site);
  
  if(!check_edges_valid(data))
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
      set_edges_valid(data);
    }  
}

//Useful for gauge fixing and hyp
void communicate_lx_su3_borders(su3 *u)
{communicate_lx_borders((char*)u,MPI_LX_SU3_BORDS_SEND,MPI_LX_SU3_BORDS_RECE,sizeof(su3));}

//Send the borders of the gauge configuration
void communicate_lx_quad_su3_borders(quad_su3 *conf)
{communicate_lx_borders((char*)conf,MPI_LX_QUAD_SU3_BORDS_SEND,MPI_LX_QUAD_SU3_BORDS_RECE,sizeof(quad_su3));}

//Send the edges: usefuls for hyp
void communicate_lx_su3_edges(su3 *u)
{communicate_lx_edges((char*)u,MPI_LX_SU3_BORDS_SEND,MPI_LX_SU3_BORDS_RECE,MPI_LX_SU3_EDGES_SEND,MPI_LX_SU3_EDGES_RECE,sizeof(su3));}

//Send the edges of the gauge configuration
void communicate_lx_quad_su3_edges(quad_su3 *conf)
{communicate_lx_edges((char*)conf,MPI_LX_QUAD_SU3_BORDS_SEND,MPI_LX_QUAD_SU3_BORDS_RECE,MPI_LX_QUAD_SU3_EDGES_SEND,MPI_LX_QUAD_SU3_EDGES_RECE,sizeof(quad_su3));}

//Send the borders of a spincolor vector
void communicate_lx_spincolor_borders(spincolor *s)
{communicate_lx_borders((char*)s,MPI_LX_SPINCOLOR_BORDS_SEND,MPI_LX_SPINCOLOR_BORDS_RECE,sizeof(spincolor));}

//Send the borders of a spin vector
void communicate_lx_spin_borders(spin *s)
{communicate_lx_borders((char*)s,MPI_LX_SPIN_BORDS_SEND,MPI_LX_SPIN_BORDS_RECE,sizeof(spin));}

//Send the borders of a color vector
void communicate_lx_color_borders(color *s)
{communicate_lx_borders((char*)s,MPI_LX_COLOR_BORDS_SEND,MPI_LX_COLOR_BORDS_RECE,sizeof(color));}

//Send the borders of a spinspin vector
void communicate_lx_spinspin_borders(spinspin *s)
{communicate_lx_borders((char*)s,MPI_LX_SPINSPIN_BORDS_SEND,MPI_LX_SPINSPIN_BORDS_RECE,sizeof(spinspin));}

//Send the borders of a spincolor_128 vector
void communicate_lx_spincolor_128_borders(spincolor_128 *s)
{communicate_lx_borders((char*)s,MPI_LX_SPINCOLOR_128_BORDS_SEND,MPI_LX_SPINCOLOR_128_BORDS_RECE,sizeof(spincolor_128));}

//Separate spincolor start/stop communication
void start_communicating_lx_spincolor_borders(int &nrequest,MPI_Request *request,spincolor *s)
{start_communicating_lx_borders(nrequest,request,(char*)s,MPI_LX_SPINCOLOR_BORDS_SEND,MPI_LX_SPINCOLOR_BORDS_RECE,sizeof(spincolor));}
void finish_communicating_lx_spincolor_borders(int &nrequest,MPI_Request *request,spincolor *s)
{finish_communicating_lx_borders(nrequest,request,(char*)s);}


///////////////////////////////////////////////// even/odd split vectors communicators ///////////////////////////////////

void finish_communicating_ev_borders(int &nrequest,MPI_Request *request,char *ev_data)
{
  tot_nissa_comm_time-=take_time();
  if(nrequest>0)
    {
      verbosity_lv3_master_printf("Waiting to finish %d communication of ev borders of vector %s\n",nrequest,get_vec_name(ev_data));
      MPI_Status status[nrequest];
      MPI_Waitall(nrequest,request,status);
      nrequest=0;
    }
  set_borders_valid(ev_data);
  set_edges_invalid(ev_data);
  tot_nissa_comm_time+=take_time();
}

//Send the borders of the data
void start_communicating_ev_borders(int &nrequest,MPI_Request *request,char *ev_data,MPI_Datatype *MPI_EV_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_EV_BORDS_RECE,int nbytes_per_site)
{
  nrequest=0;
  
  crash_if_borders_not_allocated(ev_data);
  
  if(!check_borders_valid(ev_data))
    {
      int imessage=534245;
      
      for(int mu=0;mu<3;mu++)
	if(paral_dir[mu]!=0)
	  {
	    //sending the upper border to the lower node
	    MPI_Irecv(ev_data+start_eo_bord_rece_up[mu]*nbytes_per_site,1,MPI_EV_BORDS_RECE[mu],rank_neighup[mu],imessage,cart_comm,request+(nrequest++));
	    MPI_Isend(ev_data+start_eo_bord_send_up[mu]*nbytes_per_site,1,MPI_EV_BORDS_SEND_TXY[mu],rank_neighdw[mu],imessage++,cart_comm,request+(nrequest++));
	    
	    //sending the lower border to the upper node
	    MPI_Irecv(ev_data+start_eo_bord_rece_dw[mu]*nbytes_per_site,1,MPI_EV_BORDS_RECE[mu],rank_neighdw[mu],imessage,cart_comm,request+(nrequest++));
	    MPI_Isend(ev_data+start_eo_bord_send_dw[mu]*nbytes_per_site,1,MPI_EV_BORDS_SEND_TXY[mu],rank_neighup[mu],imessage++,cart_comm,request+(nrequest++));
	  }
      
      if(paral_dir[3]!=0)
	{
	  //sending the upper border to the lower node
	  MPI_Irecv(ev_data+start_eo_bord_rece_up[3]*nbytes_per_site,1,MPI_EV_BORDS_RECE[3],rank_neighup[3],imessage,cart_comm,request+(nrequest++));
	  MPI_Isend(ev_data,1,MPI_EV_BORDS_SEND_Z[0],rank_neighdw[3],imessage++,cart_comm,request+(nrequest++));
	  
	  //sending the lower border to the upper node
	  MPI_Irecv(ev_data+start_eo_bord_rece_dw[3]*nbytes_per_site,1,MPI_EV_BORDS_RECE[3],rank_neighdw[3],imessage,cart_comm,request+(nrequest++));
	  MPI_Isend(ev_data,1,MPI_EV_BORDS_SEND_Z[1],rank_neighup[3],imessage++,cart_comm,request+(nrequest++));
	}
      
      if(nrequest>0) verbosity_lv3_master_printf("Starting communication of ev borders of vector %s\n",get_vec_name(ev_data));
    }
}

//Send the borders of the data
void communicate_ev_borders(char *ev_data,MPI_Datatype *MPI_EV_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_EV_BORDS_RECE,int nbytes_per_site)
{
  MPI_Request request[16];
  int nrequest;
  start_communicating_ev_borders(nrequest,request,ev_data,MPI_EV_BORDS_SEND_TXY,MPI_EV_BORDS_SEND_Z,MPI_EV_BORDS_RECE,nbytes_per_site);
  finish_communicating_ev_borders(nrequest,request,ev_data);
}

//Send the borders of the data
void communicate_od_borders(char *od_data,MPI_Datatype *MPI_EV_BORDS_SEND_TXY,MPI_Datatype *MPI_OD_BORDS_SEND_Z,MPI_Datatype *MPI_EV_BORDS_RECE,int nbytes_per_site)
{
  crash_if_borders_not_allocated(od_data);
  
  if(!check_borders_valid(od_data))
    {
      tot_nissa_comm_time-=take_time();

      int nrequest=0;
      int imessage=349234;
      MPI_Request request[16];
      MPI_Status status[16];
      
      for(int mu=0;mu<3;mu++)
	if(paral_dir[mu]!=0)
	  {
	    //sending the upper border to the lower node
	    MPI_Irecv(od_data+start_eo_bord_rece_up[mu]*nbytes_per_site,1,MPI_EV_BORDS_RECE[mu],rank_neighup[mu],imessage,cart_comm,request+(nrequest++));
	    MPI_Isend(od_data+start_eo_bord_send_up[mu]*nbytes_per_site,1,MPI_EV_BORDS_SEND_TXY[mu],rank_neighdw[mu],imessage++,cart_comm,request+(nrequest++));
	    
	    //sending the lower border to the upper node
	    MPI_Irecv(od_data+start_eo_bord_rece_dw[mu]*nbytes_per_site,1,MPI_EV_BORDS_RECE[mu],rank_neighdw[mu],imessage,cart_comm,request+(nrequest++));
	    MPI_Isend(od_data+start_eo_bord_send_dw[mu]*nbytes_per_site,1,MPI_EV_BORDS_SEND_TXY[mu],rank_neighup[mu],imessage++,cart_comm,request+(nrequest++));
	  }
      
      if(paral_dir[3]!=0)
	{
	  //sending the upper border to the lower node
	  MPI_Irecv(od_data+start_eo_bord_rece_up[3]*nbytes_per_site,1,MPI_EV_BORDS_RECE[3],rank_neighup[3],imessage,cart_comm,request+(nrequest++));
	  MPI_Isend(od_data,1,MPI_OD_BORDS_SEND_Z[0],rank_neighdw[3],imessage++,cart_comm,request+(nrequest++));
	  
	  //sending the lower border to the upper node
	  MPI_Irecv(od_data+start_eo_bord_rece_dw[3]*nbytes_per_site,1,MPI_EV_BORDS_RECE[3],rank_neighdw[3],imessage,cart_comm,request+(nrequest++));
	  MPI_Isend(od_data,1,MPI_OD_BORDS_SEND_Z[1],rank_neighup[3],imessage++,cart_comm,request+(nrequest++));
	}
      
      if(nrequest>0)
	{
	  verbosity_lv3_master_printf("Communicating od borders of vector %s, nrequest: %d\n",get_vec_name(od_data),nrequest);
	  MPI_Waitall(nrequest,request,status);
	}
      
      set_borders_valid(od_data);
      set_edges_invalid(od_data);
      tot_nissa_comm_time+=take_time();
    }
}

//Send the borders of the data
void communicate_eo_borders(char **data,MPI_Datatype *MPI_EO_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_OD_BORDS_SEND_Z,MPI_Datatype *MPI_EO_BORDS_RECE,int nbytes_per_site)
{
  int nrequest=0;
  MPI_Request request[32];
  int imessage=63458729;
  for(int par=0;par<2;par++)
    {
      crash_if_borders_not_allocated(data[par]);
      
      if(!check_borders_valid(data[par]))
	{
	  for(int mu=0;mu<3;mu++) //0,1,2
	    if(paral_dir[mu]!=0)
	      {
		//sending the upper border to the lower node
		imessage++;
		MPI_Irecv(data[par]+start_eo_bord_rece_up[mu]*nbytes_per_site,1,MPI_EO_BORDS_RECE[mu],rank_neighup[mu],imessage,
			  cart_comm,request+(nrequest++));
		MPI_Isend(data[par]+start_eo_bord_send_up[mu]*nbytes_per_site,1,MPI_EO_BORDS_SEND_TXY[mu],rank_neighdw[mu],imessage,
			  cart_comm,request+(nrequest++));

		//sending the lower border to the upper node
		imessage++;
		MPI_Irecv(data[par]+start_eo_bord_rece_dw[mu]*nbytes_per_site,1,MPI_EO_BORDS_RECE[mu],rank_neighdw[mu],imessage, 
			  cart_comm,request+(nrequest++));
		MPI_Isend(data[par]+start_eo_bord_send_dw[mu]*nbytes_per_site,1,MPI_EO_BORDS_SEND_TXY[mu],rank_neighup[mu],imessage,
			  cart_comm,request+(nrequest++));
	      }
	}
    }

  if(paral_dir[3]!=0)
    {
      if(!check_borders_valid(data[EVN]))
	{
	  //sending the upper border to the lower node
	  imessage++;
	  MPI_Irecv(data[EVN]+start_eo_bord_rece_up[3]*nbytes_per_site,1,MPI_EO_BORDS_RECE[3],rank_neighup[3],imessage,cart_comm,request+(nrequest++));
	  MPI_Isend(data[EVN],1,MPI_EV_BORDS_SEND_Z[0],rank_neighdw[3],imessage,cart_comm,request+(nrequest++));
	  
	  //sending the lower border to the upper node
	  imessage++;
	  MPI_Irecv(data[EVN]+start_eo_bord_rece_dw[3]*nbytes_per_site,1,MPI_EO_BORDS_RECE[3],rank_neighdw[3],imessage,cart_comm,request+(nrequest++));
	  MPI_Isend(data[EVN],1,MPI_EV_BORDS_SEND_Z[1],rank_neighup[3],imessage,cart_comm,request+(nrequest++));
	}
      
      if(!check_borders_valid(data[ODD]))
	{
	  //sending the upper border to the lower node
	  imessage++;
	  MPI_Irecv(data[ODD]+start_eo_bord_rece_up[3]*nbytes_per_site,1,MPI_EO_BORDS_RECE[3],rank_neighup[3],imessage,cart_comm,request+(nrequest++));
	  MPI_Isend(data[ODD],1,MPI_OD_BORDS_SEND_Z[0],rank_neighdw[3],imessage,cart_comm,request+(nrequest++));
	  
	  //sending the lower border to the upper node
	  imessage++;
	  MPI_Irecv(data[ODD]+start_eo_bord_rece_dw[3]*nbytes_per_site,1,MPI_EO_BORDS_RECE[3],rank_neighdw[3],imessage,cart_comm,request+(nrequest++));
	  MPI_Isend(data[ODD],1,MPI_OD_BORDS_SEND_Z[1],rank_neighup[3],imessage,cart_comm,request+(nrequest++));
	}
    }
  
  tot_nissa_comm_time-=take_time();
  if(nrequest>0)
    {
      if(nrequest>0) verbosity_lv3_master_printf("Starting communication of borders of vector %s\n",get_vec_name(data[0]));
      MPI_Status status[nrequest];
      MPI_Waitall(nrequest,request,status);
      nrequest=0;
    }
  tot_nissa_comm_time+=take_time();

  set_borders_valid(data[EVN]);
  set_borders_valid(data[ODD]);
}

//Send the borders of the gauge configuration
void communicate_eo_quad_su3_borders(quad_su3 **eo_conf)
{communicate_eo_borders((char**)eo_conf,MPI_EO_QUAD_SU3_BORDS_SEND_TXY,MPI_EV_QUAD_SU3_BORDS_SEND_Z,MPI_OD_QUAD_SU3_BORDS_SEND_Z,MPI_EO_QUAD_SU3_BORDS_RECE,sizeof(quad_su3));}

//Send the borders of an even color
void communicate_ev_color_borders(color *ev)
{communicate_ev_borders((char*)ev,MPI_EO_COLOR_BORDS_SEND_TXY,MPI_EV_COLOR_BORDS_SEND_Z,MPI_EO_COLOR_BORDS_RECE,sizeof(color));}
void communicate_od_color_borders(color *od)
{communicate_od_borders((char*)od,MPI_EO_COLOR_BORDS_SEND_TXY,MPI_OD_COLOR_BORDS_SEND_Z,MPI_EO_COLOR_BORDS_RECE,sizeof(color));}
void communicate_eo_color_borders(color **eos)
{communicate_eo_borders((char**)eos,MPI_EO_COLOR_BORDS_SEND_TXY,MPI_EV_COLOR_BORDS_SEND_Z,MPI_OD_COLOR_BORDS_SEND_Z,MPI_EO_COLOR_BORDS_RECE,sizeof(color));}
void start_communicating_ev_color_borders(int &nrequest,MPI_Request *request,color *ev)
{start_communicating_ev_borders(nrequest,request,(char*)ev,MPI_EO_COLOR_BORDS_SEND_TXY,MPI_EV_COLOR_BORDS_SEND_Z,MPI_EO_COLOR_BORDS_RECE,sizeof(color));}
void finish_communicating_ev_color_borders(int &nrequest,MPI_Request *request,color *ev)
{finish_communicating_ev_borders(nrequest,request,(char*)ev);}

//Send the borders of an even spincolor
void communicate_ev_spincolor_borders(spincolor *ev)
{communicate_ev_borders((char*)ev,MPI_EO_SPINCOLOR_BORDS_SEND_TXY,MPI_EV_SPINCOLOR_BORDS_SEND_Z,MPI_EO_SPINCOLOR_BORDS_RECE,sizeof(spincolor));}
void communicate_od_spincolor_borders(spincolor *od)
{communicate_od_borders((char*)od,MPI_EO_SPINCOLOR_BORDS_SEND_TXY,MPI_OD_SPINCOLOR_BORDS_SEND_Z,MPI_EO_SPINCOLOR_BORDS_RECE,sizeof(spincolor));}
void start_communicating_ev_spincolor_borders(int &nrequest,MPI_Request *request,spincolor *ev)
{start_communicating_ev_borders(nrequest,request,(char*)ev,MPI_EO_SPINCOLOR_BORDS_SEND_TXY,MPI_EV_SPINCOLOR_BORDS_SEND_Z,MPI_EO_SPINCOLOR_BORDS_RECE,sizeof(spincolor));}
void finish_communicating_ev_spincolor_borders(int &nrequest,MPI_Request *request,spincolor *ev)
{finish_communicating_ev_borders(nrequest,request,(char*)ev);}
void communicate_ev_spin_borders(spin *ev)
{communicate_ev_borders((char*)ev,MPI_EO_SPIN_BORDS_SEND_TXY,MPI_EV_SPIN_BORDS_SEND_Z,MPI_EO_SPIN_BORDS_RECE,sizeof(spin));}
void communicate_od_spin_borders(spin *od)
{communicate_od_borders((char*)od,MPI_EO_SPIN_BORDS_SEND_TXY,MPI_OD_SPIN_BORDS_SEND_Z,MPI_EO_SPIN_BORDS_RECE,sizeof(spin));}

//128 bit version
void communicate_ev_spincolor_128_borders(spincolor_128 *ev)
{communicate_ev_borders((char*)ev,MPI_EO_SPINCOLOR_128_BORDS_SEND_TXY,MPI_EV_SPINCOLOR_128_BORDS_SEND_Z,MPI_EO_SPINCOLOR_128_BORDS_RECE,sizeof(spincolor_128));}
void communicate_od_spincolor_128_borders(spincolor_128 *od)
{communicate_od_borders((char*)od,MPI_EO_SPINCOLOR_128_BORDS_SEND_TXY,MPI_OD_SPINCOLOR_128_BORDS_SEND_Z,MPI_EO_SPINCOLOR_128_BORDS_RECE,sizeof(spincolor_128));}

//Send the edges of eo vector
void communicate_eo_edges(char **data,MPI_Datatype *MPI_EO_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_OD_BORDS_SEND_Z,MPI_Datatype *MPI_EO_BORDS_RECE,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site)
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
  
  if(!check_edges_valid(data[EVN])||!check_edges_valid(data[ODD]))
    {
      communicate_eo_borders(data,MPI_EO_BORDS_SEND_TXY,MPI_EV_BORDS_SEND_Z,MPI_OD_BORDS_SEND_Z,MPI_EO_BORDS_RECE,nbytes_per_site);
      
      verbosity_lv3_master_printf("Communicating edges\n");
      
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
	  
	  set_edges_valid(data[par]);
	}
      
      if(nrequest!=nrequest_tot) crash("something went wrong");
      
      if(nrequest>0) MPI_Waitall(nrequest,request,status);
    }
}

//Send the edges of the gauge configuration
void communicate_eo_quad_su3_edges(quad_su3 **conf)
{communicate_eo_edges((char**)conf,MPI_EO_QUAD_SU3_BORDS_SEND_TXY,MPI_EV_QUAD_SU3_BORDS_SEND_Z,MPI_OD_QUAD_SU3_BORDS_SEND_Z,MPI_EO_QUAD_SU3_BORDS_RECE,MPI_EO_QUAD_SU3_EDGES_SEND,MPI_EO_QUAD_SU3_EDGES_RECE,sizeof(quad_su3));}
