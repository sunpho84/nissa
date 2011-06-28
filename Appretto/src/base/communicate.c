#pragma once

#include "types/su3.c"
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
void communicate_su3_borders(su3 *u)
{communicate_lx_borders((char*)u,MPI_SU3_BORD_SEND,MPI_SU3_BORD_RECE,sizeof(su3));}
	 
//Send the borders of the gauge configuration
void communicate_gauge_borders(quad_su3 *conf)
{communicate_lx_borders((char*)conf,MPI_GAUGE_BORD_SEND,MPI_GAUGE_BORD_RECE,sizeof(quad_su3));}

//Send the edges of the gauge configuration
void communicate_gauge_edges(quad_su3 *conf)
{communicate_lx_edges((char*)conf,MPI_GAUGE_EDGE_SEND,MPI_GAUGE_EDGE_RECE,sizeof(quad_su3));}

//Send the borders of a spincolor vector
void communicate_lx_spincolor_borders(spincolor *s)
{communicate_lx_borders((char*)s,MPI_LXSPINCOLOR_BORD_SEND,MPI_LXSPINCOLOR_BORD_RECE,sizeof(spincolor));}

//apply the projector over the internal border
void apply_lxspincolor_border_projector(redspincolor *p,spincolor *s)
{
#ifdef BGP
  bgp_complex S00,S01,S02;
  bgp_complex S10,S11,S12;
  bgp_complex S20,S21,S22;
  bgp_complex S30,S31,S32;

  for(int ibord=0;ibord<loc_bord;ibord++)
    {
      int idir=dir_of_bordlx[ibord];
      int ivol=loclx_of_bordlx[ibord];

      bgp_load_color(S00,S01,S02,s[ivol][0]);
      bgp_load_color(S10,S11,S12,s[ivol][1]);
      bgp_load_color(S20,S21,S22,s[ivol][2]);
      bgp_load_color(S30,S31,S32,s[ivol][3]);
  
      switch(idir)
	{
	case 0:
	  bgp_summassign_color(S00,S01,S02, S20,S21,S22);
	  bgp_summassign_color(S10,S11,S12, S30,S31,S32);
	  break;
	case 1:
	  bgp_subtassign_color(S00,S01,S02, S20,S21,S22);
	  bgp_subtassign_color(S10,S11,S12, S30,S31,S32);
	  break;
	case 2:
	  bgp_summassign_icolor(S00,S01,S02, S30,S31,S32);
	  bgp_summassign_icolor(S10,S11,S12, S20,S21,S22);
	  break;
	case 3:
	  bgp_subtassign_icolor(S00,S01,S02, S30,S31,S32);
	  bgp_subtassign_icolor(S10,S11,S12, S20,S21,S22);
	  break;
	case 4:
	  bgp_summassign_color(S00,S01,S02, S30,S31,S32);
	  bgp_subtassign_color(S10,S11,S12, S20,S21,S22);
	  break;
	case 5:
	  bgp_subtassign_color(S00,S01,S02, S30,S31,S32);
	  bgp_summassign_color(S10,S11,S12, S20,S21,S22);
	  break;
	case 6:
	  bgp_summassign_icolor(S00,S01,S02, S20,S21,S22);
	  bgp_subtassign_icolor(S10,S11,S12, S30,S31,S32);
	  break;
	case 7:
	  bgp_subtassign_icolor(S00,S01,S02, S20,S21,S22);
	  bgp_summassign_icolor(S10,S11,S12, S30,S31,S32);
	  break;
	}
      
      bgp_save_color(p[ibord][0],S00,S01,S02);
      bgp_save_color(p[ibord][1],S10,S11,S12);
    }
#else
  int ib=0;
  
  /////////// border for the backward derivate, to be sent backward //////////

  if(paral_dir[0]) //Backward 0
    for(int X=loc_vol;X<loc_vol+bord_offset[1];X++)
      {
	int Xup0=loclx_neighup[X][0];
	color_summ(p[ib][0],s[Xup0][0],s[Xup0][2]);
	color_summ(p[ib++][1],s[Xup0][1],s[Xup0][3]);
      }
  if(paral_dir[1]) //Backward 1
    for(int X=loc_vol+bord_offset[1];X<loc_vol+bord_offset[2];X++)
      {
	int Xup1=loclx_neighup[X][1];
	color_isumm(p[ib][0],s[Xup1][0],s[Xup1][3]);
	color_isumm(p[ib++][1],s[Xup1][1],s[Xup1][2]);
      }
  if(paral_dir[2]) //Backward 2
    for(int X=loc_vol+bord_offset[2];X<loc_vol+bord_offset[3];X++)
      {
	int Xup2=loclx_neighup[X][2];
	color_summ(p[ib][0],s[Xup2][0],s[Xup2][3]);
	color_subt(p[ib++][1],s[Xup2][1],s[Xup2][2]);
      }
  if(paral_dir[3]) //Backward 3
    for(int X=loc_vol+bord_offset[3];X<loc_vol+loc_bord/2;X++)
      {
	int Xup3=loclx_neighup[X][3];
	color_isumm(p[ib][0],s[Xup3][0],s[Xup3][2]);
	color_isubt(p[ib++][1],s[Xup3][1],s[Xup3][3]);
      }
  
  /////////// border for the forward derivate, to be sent backward //////////
  
  if(paral_dir[0]) //Forward 0
    for(int X=loc_vol+loc_bord/2;X<loc_vol+loc_bord/2+bord_offset[1];X++)
      {
	int Xdw0=loclx_neighdw[X][0];
	color_subt(p[ib][0],s[Xdw0][0],s[Xdw0][2]);
	color_subt(p[ib++][1],s[Xdw0][1],s[Xdw0][3]);
      }
  if(paral_dir[1]) //Forward 1
    for(int X=loc_vol+loc_bord/2+bord_offset[1];X<loc_vol+loc_bord/2+bord_offset[2];X++)
      {
	int Xdw1=loclx_neighdw[X][1];
	color_isubt(p[ib][0],s[Xdw1][0],s[Xdw1][3]);
	color_isubt(p[ib++][1],s[Xdw1][1],s[Xdw1][2]);
      }
  if(paral_dir[2]) //Forward 2
    for(int X=loc_vol+loc_bord/2+bord_offset[2];X<loc_vol+loc_bord/2+bord_offset[3];X++)
      {
	int Xdw2=loclx_neighdw[X][2];
	color_subt(p[ib][0],s[Xdw2][0],s[Xdw2][3]);
	color_summ(p[ib++][1],s[Xdw2][1],s[Xdw2][2]);
      }
  if(paral_dir[3]) //Forward 3
    for(int X=loc_vol+loc_bord/2+bord_offset[3];X<loc_vol+loc_bord;X++)
      {
	int Xdw3=loclx_neighdw[X][3];
	color_isubt(p[ib][0],s[Xdw3][0],s[Xdw3][2]);
	color_isumm(p[ib++][1],s[Xdw3][1],s[Xdw3][3]);
      }
#endif
}

//Send the reduced border of a spincolor vector
void start_communicate_lx_redspincolor_borders(redspincolor *in,redspincolor *out,MPI_Request *request)
{
  int nrequest=0;

  for(int i=0;i<4;i++) //sending the border containing backward derivate to the forward node
    if(paral_dir[i])
      {
	MPI_Irecv((void*)(in +bord_offset[i]           ),1,MPI_LXREDSPINCOLOR_BORD[i],rank_neighdw[i],83+i,cart_comm,&request[nrequest++]);
	MPI_Isend((void*)(out+bord_offset[i]+loc_bord/2),1,MPI_LXREDSPINCOLOR_BORD[i],rank_neighup[i],83+i,cart_comm,&request[nrequest++]);
      }
  
  for(int i=0;i<4;i++) //sending the border containing forward derivate to the backward node
    if(paral_dir[i])
      {
	MPI_Irecv((void*)(in +bord_offset[i]+loc_bord/2),1,MPI_LXREDSPINCOLOR_BORD[i],rank_neighup[i],87+i,cart_comm,&request[nrequest++]);
	MPI_Isend((void*)(out+bord_offset[i]           ),1,MPI_LXREDSPINCOLOR_BORD[i],rank_neighdw[i],87+i,cart_comm,&request[nrequest++]);
      }
}
