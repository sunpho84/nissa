#pragma once

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

and following, the negative direction borders
*/

//Send the borders of the gauge configuration
void communicate_gauge_borders(quad_su3 *conf)
{
  int nrequest=0;
  MPI_Request request[16];
  MPI_Status status[16];
  int send,rece;

  for(int i=0;i<4;i++)
    {
      //sending the upper border to the lower node
      send=loclx_of_coord_list(0,0,0,0);
      rece=loc_vol+bord_offset[i]+loc_bord/2;
      MPI_Isend(conf[send],1,MPI_GAUGE_SLICE_SEND[i],rank_neighdw[i],83+i,
		cart_comm,&request[nrequest++]);
      MPI_Irecv(conf[rece],1,MPI_GAUGE_SLICE_RECE[i],rank_neighup[i],83+i, 
		cart_comm,&request[nrequest++]);
      
      //sending the lower border to the upper node
      int x[4];
      for(int jdir=0;jdir<4;jdir++)
	if(jdir==i) x[jdir]=loc_size[i]-1;
	else x[jdir]=0;
      send=loclx_of_coord(x);
      rece=loc_vol+bord_offset[i];
      MPI_Isend(conf[send],1,MPI_GAUGE_SLICE_SEND[i],rank_neighup[i],87+i,
		cart_comm,&request[nrequest++]);
      MPI_Irecv(conf[rece],1,MPI_GAUGE_SLICE_RECE[i],rank_neighdw[i],87+i, 
		  cart_comm,&request[nrequest++]);
    }
  
  MPI_Waitall(nrequest,request,status);
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

//Send the edges of the gauge configuration
void communicate_gauge_edges(quad_su3 *conf)
{
  int nrequest=0;
  MPI_Request request[48];
  MPI_Status status[48];
  int send,rece;
  int imessage=0;
  int x[4]={0,0,0,0};
  
  for(int idir=0;idir<4;idir++)
    for(int jdir=idir+1;jdir<4;jdir++)
      {
	int iedge=edge_numb[idir][jdir];
	int pos_edge_offset;

	//take the starting point of the border
	x[jdir]=loc_size[jdir]-1;
	pos_edge_offset=bordlx_of_coord(x,idir);
	x[jdir]=0;

	//Send the i-j- internal edge to the j- site as i-j+ external edge
	send=loc_vol+bord_offset[idir];
	rece=loc_vol+loc_bord+edge_offset[iedge]+loc_edge/4;
	MPI_Isend(conf[send],1,MPI_GAUGE_EDGE_SEND[iedge],rank_neighdw[jdir],83+imessage,
		  cart_comm,&request[nrequest++]);
	MPI_Irecv(conf[rece],1,MPI_GAUGE_EDGE_RECE[iedge],rank_neighup[jdir],83+imessage,
		  cart_comm,&request[nrequest++]);
	imessage++;
	
	//Send the i-j+ internal edge to the j+ site as i-j- external edge
	send=loc_vol+bord_offset[idir]+pos_edge_offset;
	rece=loc_vol+loc_bord+edge_offset[iedge];
	MPI_Isend(conf[send],1,MPI_GAUGE_EDGE_SEND[iedge],rank_neighup[jdir],83+imessage,
		  cart_comm,&request[nrequest++]);
	MPI_Irecv(conf[rece],1,MPI_GAUGE_EDGE_RECE[iedge],rank_neighdw[jdir],83+imessage,
		  cart_comm,&request[nrequest++]);
	imessage++;

	//Send the i+j- internal edge to the j- site as i+j+ external edge
	send=loc_vol+bord_offset[idir]+loc_bord/2;
	rece=loc_vol+loc_bord+edge_offset[iedge]+3*loc_edge/4;
	MPI_Isend(conf[send],1,MPI_GAUGE_EDGE_SEND[iedge],rank_neighdw[jdir],83+imessage,
		  cart_comm,&request[nrequest++]);
	MPI_Irecv(conf[rece],1,MPI_GAUGE_EDGE_RECE[iedge],rank_neighup[jdir],83+imessage,
		  cart_comm,&request[nrequest++]);
	imessage++;
	
	//Send the i+j+ internal edge to the j+ site as i+j- external edge
	send=loc_vol+bord_offset[idir]+loc_bord/2+pos_edge_offset;
	rece=loc_vol+loc_bord+edge_offset[iedge]+loc_edge/2;
	MPI_Isend(conf[send],1,MPI_GAUGE_EDGE_SEND[iedge],rank_neighup[jdir],83+imessage,
		  cart_comm,&request[nrequest++]);
	MPI_Irecv(conf[rece],1,MPI_GAUGE_EDGE_RECE[iedge],rank_neighdw[jdir],83+imessage,
		  cart_comm,&request[nrequest++]);
	imessage++;
      }

  MPI_Waitall(nrequest,request,status);
}
	 

