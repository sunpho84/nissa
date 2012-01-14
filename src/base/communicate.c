#pragma once

#include "bgp_instructions.c"

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
_______________t__=__0_______________|||_______________t__=__1_______________|||_______________t__=__2_______________|
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

and following, the negative direction borders
*/

//Send the borders of the data
void communicate_lx_borders(char *data,MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *MPI_BORD_RECE,int nbytes_per_site)
{
  if(debug_lvl>1) check_borders_allocated(data);
  
  int nrequest=0;
  MPI_Request request[16];
  MPI_Status status[16];

  for(int i=0;i<4;i++)
    if(paral_dir[i]!=0)
      {
	//sending the upper border to the lower node
	MPI_Irecv((void*)(data+start_lx_bord_rece_up[i]*nbytes_per_site),1,MPI_BORD_RECE[i],rank_neighup[i],83+i, 
		  cart_comm,&request[nrequest++]);
	MPI_Isend((void*)(data+start_lx_bord_send_up[i]*nbytes_per_site),1,MPI_BORD_SEND[i],rank_neighdw[i],83+i,
		  cart_comm,&request[nrequest++]);
	
	//sending the lower border to the upper node
	MPI_Irecv((void*)(data+start_lx_bord_rece_dw[i]*nbytes_per_site),1,MPI_BORD_RECE[i],rank_neighdw[i],87+i, 
		  cart_comm,&request[nrequest++]);
	MPI_Isend((void*)(data+start_lx_bord_send_dw[i]*nbytes_per_site),1,MPI_BORD_SEND[i],rank_neighup[i],87+i,
		  cart_comm,&request[nrequest++]);
      }
  
  if(nrequest>0) MPI_Waitall(nrequest,request,status);
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
void communicate_lx_edges(char *data,MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,int nbytes_per_site)
{
  if(debug_lvl>1) check_edges_allocated(data);
  
  int nrequest=0;
  MPI_Request request[48];
  MPI_Status status[48];
  int send,rece;
  int imessage=0;
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
	  
	  //Send the i-j- internal edge to the j- site as i-j+ external edge
	  send=(loc_vol+bord_offset[idir])*nbytes_per_site;
	  rece=(loc_vol+loc_bord+edge_offset[iedge]+loc_edge/4)*nbytes_per_site;
	  MPI_Irecv((void*)(data+rece),1,MPI_EDGE_RECE[iedge],rank_neighup[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  MPI_Isend((void*)(data+send),1,MPI_EDGE_SEND[iedge],rank_neighdw[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  imessage++;
	  
	  //Send the i-j+ internal edge to the j+ site as i-j- external edge
	  send=(loc_vol+bord_offset[idir]+pos_edge_offset)*nbytes_per_site;
	  rece=(loc_vol+loc_bord+edge_offset[iedge])*nbytes_per_site;
	  MPI_Irecv((void*)(data+rece),1,MPI_EDGE_RECE[iedge],rank_neighdw[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  MPI_Isend((void*)(data+send),1,MPI_EDGE_SEND[iedge],rank_neighup[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  imessage++;
	  
	  //Send the i+j- internal edge to the j- site as i+j+ external edge
	  send=(loc_vol+bord_offset[idir]+loc_bord/2)*nbytes_per_site;
	  rece=(loc_vol+loc_bord+edge_offset[iedge]+3*loc_edge/4)*nbytes_per_site;
	  MPI_Irecv((void*)(data+rece),1,MPI_EDGE_RECE[iedge],rank_neighup[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  MPI_Isend((void*)(data+send),1,MPI_EDGE_SEND[iedge],rank_neighdw[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  imessage++;
	  
	  //Send the i+j+ internal edge to the j+ site as i+j- external edge
	  send=(loc_vol+bord_offset[idir]+loc_bord/2+pos_edge_offset)*nbytes_per_site;
	  rece=(loc_vol+loc_bord+edge_offset[iedge]+loc_edge/2)*nbytes_per_site;
	  MPI_Irecv((void*)(data+rece),1,MPI_EDGE_RECE[iedge],rank_neighdw[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  MPI_Isend((void*)(data+send),1,MPI_EDGE_SEND[iedge],rank_neighup[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  imessage++;
	}
  
  if(nrequest>0) MPI_Waitall(nrequest,request,status);
}

//Useful for gauge fixing
void communicate_lx_su3_borders(su3 *u)
{communicate_lx_borders((char*)u,MPI_LX_SU3_BORD_SEND,MPI_LX_SU3_BORD_RECE,sizeof(su3));}

//Send the borders of the gauge configuration
void communicate_lx_gauge_borders(quad_su3 *conf)
{communicate_lx_borders((char*)conf,MPI_LX_GAUGE_BORD_SEND,MPI_LX_GAUGE_BORD_RECE,sizeof(quad_su3));}

//Send the edges of the gauge configuration
void communicate_lx_gauge_edges(quad_su3 *conf)
{communicate_lx_edges((char*)conf,MPI_LX_GAUGE_EDGE_SEND,MPI_LX_GAUGE_EDGE_RECE,sizeof(quad_su3));}

//Send the borders of a spincolor vector
void communicate_lx_spincolor_borders(spincolor *s)
{communicate_lx_borders((char*)s,MPI_LX_SPINCOLOR_BORD_SEND,MPI_LX_SPINCOLOR_BORD_RECE,sizeof(spincolor));}


///////////////////////////////////////////////// even/odd split vectors communicators ///////////////////////////////////

//Send the borders of the data
void communicate_ev_borders(char *ev_data,MPI_Datatype *MPI_EV_BORD_SEND,MPI_Datatype *MPI_EV_BORD_RECE,int nbytes_per_site)
{
  if(debug_lvl>1) check_borders_allocated(ev_data);
  
  int nrequest=0;
  MPI_Request request[16];
  MPI_Status status[16];

  for(int mu=0;mu<4;mu++)
    if(paral_dir[mu]!=0)
      {
	//sending the upper border to the lower node
	MPI_Irecv((void*)(ev_data+start_eo_bord_rece_up[mu]*nbytes_per_site),1,MPI_EV_BORD_RECE[mu],rank_neighup[mu],83+mu,
		  cart_comm,&request[nrequest++]);
	MPI_Isend((void*)(ev_data+start_eo_bord_send_up[mu]*nbytes_per_site),1,MPI_EV_BORD_SEND[mu],rank_neighdw[mu],83+mu,
		  cart_comm,&request[nrequest++]);
	
	//sending the lower border to the upper node
	MPI_Irecv((void*)(ev_data+start_eo_bord_rece_dw[mu]*nbytes_per_site),1,MPI_EV_BORD_RECE[mu],rank_neighdw[mu],91+mu, 
		  cart_comm,&request[nrequest++]);
	MPI_Isend((void*)(ev_data+start_eo_bord_send_dw[mu]*nbytes_per_site),1,MPI_EV_BORD_SEND[mu],rank_neighup[mu],91+mu,
		  cart_comm,&request[nrequest++]);
      }
  
  if(nrequest>0) MPI_Waitall(nrequest,request,status);
}

//Send the borders of the data
void communicate_eo_borders(char *ev_data,char *od_data,MPI_Datatype *MPI_EO_BORD_SEND,MPI_Datatype *MPI_EO_BORD_RECE,int nbytes_per_site)
{
  char *data[2]={ev_data,od_data};
  
  if(debug_lvl>1) check_borders_allocated(data);
  
  int nrequest=0;
  MPI_Request request[32];
  MPI_Status status[32];

  for(int par=0;par<2;par++)
    for(int mu=0;mu<4;mu++)
      if(paral_dir[mu]!=0)
	{
	  //sending the upper border to the lower node
	  MPI_Irecv((void*)(data[par]+start_eo_bord_rece_up[mu]*nbytes_per_site),1,MPI_EO_BORD_RECE[mu],rank_neighup[mu],83+par*4+mu,
		    cart_comm,&request[nrequest++]);
	  MPI_Isend((void*)(data[par]+start_eo_bord_send_up[mu]*nbytes_per_site),1,MPI_EO_BORD_SEND[mu],rank_neighdw[mu],83+par*4+mu,
		    cart_comm,&request[nrequest++]);
	  
	  //sending the lower border to the upper node
	  MPI_Irecv((void*)(data[par]+start_eo_bord_rece_dw[mu]*nbytes_per_site),1,MPI_EO_BORD_RECE[mu],rank_neighdw[mu],91+par*4+mu, 
		    cart_comm,&request[nrequest++]);
	  MPI_Isend((void*)(data[par]+start_eo_bord_send_dw[mu]*nbytes_per_site),1,MPI_EO_BORD_SEND[mu],rank_neighup[mu],91+par*4+mu,
		    cart_comm,&request[nrequest++]);
	}
  
  if(nrequest>0) MPI_Waitall(nrequest,request,status);
}

//Send the borders of the gauge configuration
void communicate_eo_gauge_borders(quad_su3 *ev_conf,quad_su3 *od_conf)
{communicate_eo_borders((char*)ev_conf,(char*)od_conf,MPI_EO_GAUGE_BORD_SEND,MPI_EO_GAUGE_BORD_RECE,sizeof(quad_su3));}

//Send the borders of an even color
void communicate_ev_color_borders(color *ev)
{communicate_ev_borders((char*)ev,MPI_EO_COLOR_BORD_SEND,MPI_EO_COLOR_BORD_RECE,sizeof(color));}
void communicate_od_color_borders(color *ev)
{communicate_ev_borders((char*)ev,MPI_EO_COLOR_BORD_SEND,MPI_EO_COLOR_BORD_RECE,sizeof(color));}
