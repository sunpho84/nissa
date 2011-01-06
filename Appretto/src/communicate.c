#pragma once

#include "su3.c"

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
	MPI_Isend((void*)(data+start_lx_bord_send_up[i]*nbytes_per_site),1,MPI_BORD_SEND[i],rank_neighdw[i],83+i,
		  cart_comm,&request[nrequest++]);
	MPI_Irecv((void*)(data+start_lx_bord_rece_up[i]*nbytes_per_site),1,MPI_BORD_RECE[i],rank_neighup[i],83+i, 
		  cart_comm,&request[nrequest++]);
	
	//sending the lower border to the upper node
	MPI_Isend((void*)(data+start_lx_bord_send_dw[i]*nbytes_per_site),1,MPI_BORD_SEND[i],rank_neighup[i],87+i,
		  cart_comm,&request[nrequest++]);
	MPI_Irecv((void*)(data+start_lx_bord_rece_dw[i]*nbytes_per_site),1,MPI_BORD_RECE[i],rank_neighdw[i],87+i, 
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
	  MPI_Isend((void*)(data+send),1,MPI_EDGE_SEND[iedge],rank_neighdw[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  MPI_Irecv((void*)(data+rece),1,MPI_EDGE_RECE[iedge],rank_neighup[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  imessage++;
	  
	  //Send the i-j+ internal edge to the j+ site as i-j- external edge
	  send=(loc_vol+bord_offset[idir]+pos_edge_offset)*nbytes_per_site;
	  rece=(loc_vol+loc_bord+edge_offset[iedge])*nbytes_per_site;
	  MPI_Isend((void*)(data+send),1,MPI_EDGE_SEND[iedge],rank_neighup[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  MPI_Irecv((void*)(data+rece),1,MPI_EDGE_RECE[iedge],rank_neighdw[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  imessage++;
	  
	  //Send the i+j- internal edge to the j- site as i+j+ external edge
	  send=(loc_vol+bord_offset[idir]+loc_bord/2)*nbytes_per_site;
	  rece=(loc_vol+loc_bord+edge_offset[iedge]+3*loc_edge/4)*nbytes_per_site;
	  MPI_Isend((void*)(data+send),1,MPI_EDGE_SEND[iedge],rank_neighdw[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  MPI_Irecv((void*)(data+rece),1,MPI_EDGE_RECE[iedge],rank_neighup[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  imessage++;
	  
	  //Send the i+j+ internal edge to the j+ site as i+j- external edge
	  send=(loc_vol+bord_offset[idir]+loc_bord/2+pos_edge_offset)*nbytes_per_site;
	  rece=(loc_vol+loc_bord+edge_offset[iedge]+loc_edge/2)*nbytes_per_site;
	  MPI_Isend((void*)(data+send),1,MPI_EDGE_SEND[iedge],rank_neighup[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  MPI_Irecv((void*)(data+rece),1,MPI_EDGE_RECE[iedge],rank_neighdw[jdir],83+imessage,
		    cart_comm,&request[nrequest++]);
	  imessage++;
	}
  
  if(nrequest>0) MPI_Waitall(nrequest,request,status);
}
	 
//Send the borders of the gauge configuration
void communicate_gauge_borders(quad_su3 *conf)
{
  communicate_lx_borders((char*)conf,MPI_GAUGE_BORD_SEND,MPI_GAUGE_BORD_RECE,sizeof(quad_su3));
}

//Send the edges of the gauge configuration
void communicate_gauge_edges(quad_su3 *conf)
{
  communicate_lx_edges((char*)conf,MPI_GAUGE_EDGE_SEND,MPI_GAUGE_EDGE_RECE,sizeof(quad_su3));
}

//Send the borders of a spincolor vector
void communicate_lx_spincolor_borders(spincolor *s)
{
  communicate_lx_borders((char*)s,MPI_LXSPINCOLOR_BORD_SEND,MPI_LXSPINCOLOR_BORD_RECE,sizeof(spincolor));
}

//apply the projector over the internal border
void apply_lxspincolor_border_projector(redspincolor *p,spincolor *s)
{
#ifdef BGP
#pragma disjoint(*s,*p)
  static double _Complex S00,S01,S02;
  static double _Complex S10,S11,S12;
  static double _Complex S20,S21,S22;
  static double _Complex S30,S31,S32;

  static double _Complex P00,P01,P02;
  static double _Complex P10,P11,P12;

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
	  bgp_color_summ(P00,P01,P02, S00,S01,S02, S20,S21,S22);
	  bgp_color_summ(P10,P11,P12, S10,S11,S12, S30,S31,S32);
	  break;
	case 1:
	  bgp_color_subt(P00,P01,P02, S00,S01,S02, S20,S21,S22);
	  bgp_color_subt(P10,P11,P12, S10,S11,S12, S30,S31,S32);
	  break;
	case 2:
	  bgp_color_isumm(P00,P01,P02, S00,S01,S02, S30,S31,S32);
	  bgp_color_isumm(P10,P11,P12, S10,S11,S12, S20,S21,S22);
	  break;
	case 3:
	  bgp_color_isubt(P00,P01,P02, S00,S01,S02, S30,S31,S32);
	  bgp_color_isubt(P10,P11,P12, S10,S11,S12, S20,S21,S22);
	  break;
	case 4:
	  bgp_color_summ(P00,P01,P02, S00,S01,S02, S30,S31,S32);
	  bgp_color_subt(P10,P11,P12, S10,S11,S12, S20,S21,S22);
	  break;
	case 5:
	  bgp_color_subt(P00,P01,P02, S00,S01,S02, S30,S31,S32);
	  bgp_color_summ(P10,P11,P12, S10,S11,S12, S20,S21,S22);
	  break;
	case 6:
	  bgp_color_isumm(P00,P01,P02, S00,S01,S02, S20,S21,S22);
	  bgp_color_isubt(P10,P11,P12, S10,S11,S12, S30,S31,S32);
	  break;
	case 7:
	  bgp_color_isubt(P00,P01,P02, S00,S01,S02, S20,S21,S22);
	  bgp_color_isumm(P10,P11,P12, S10,S11,S12, S30,S31,S32);
	  break;
	}
      
      bgp_save_color(p[ibord][0],P00,P01,P02);
      bgp_save_color(p[ibord][1],P10,P11,P12);
    }
#else
  for(int ibord=0;ibord<loc_bord;ibord++)
    {
      int idir=dir_of_bordlx[ibord];
      int ivol=loclx_of_bordlx[ibord];
      switch(idir)
	{
	case 0:
	  color_summ(p[ibord][0],s[ivol][0],s[ivol][2]);
	  color_summ(p[ibord][1],s[ivol][1],s[ivol][3]);
	  break;
	case 1:
	  color_subt(p[ibord][0],s[ivol][0],s[ivol][2]);
	  color_subt(p[ibord][1],s[ivol][1],s[ivol][3]);
	  break;
	case 2:
	  color_isumm(p[ibord][0],s[ivol][0],s[ivol][3]);
	  color_isumm(p[ibord][1],s[ivol][1],s[ivol][2]);
	  break;
	case 3:
	  color_isubt(p[ibord][0],s[ivol][0],s[ivol][3]);
	  color_isubt(p[ibord][1],s[ivol][1],s[ivol][2]);
	  break;
	case 4:
	  color_summ(p[ibord][0],s[ivol][0],s[ivol][3]);
	  color_subt(p[ibord][1],s[ivol][1],s[ivol][2]);
	  break;
	case 5:
	  color_subt(p[ibord][0],s[ivol][0],s[ivol][3]);
	  color_summ(p[ibord][1],s[ivol][1],s[ivol][2]);
	  break;
	case 6:
	  color_isumm(p[ibord][0],s[ivol][0],s[ivol][2]);
	  color_isubt(p[ibord][1],s[ivol][1],s[ivol][3]);
	  break;
	case 7:
	  color_isubt(p[ibord][0],s[ivol][0],s[ivol][2]);
	  color_isumm(p[ibord][1],s[ivol][1],s[ivol][3]);
	  break;
	}
    }
#endif
}

//apply the projector over the internal border
void expand_projected_lxspincolor_border(spincolor *s,redspincolor *p)
{
  for(int ibord=0;ibord<loc_bord;ibord++)
    {
      memcpy(&(s[loc_vol+ibord][0]),p[ibord],sizeof(color)*2);
      memset(&(s[loc_vol+ibord][2]),0,sizeof(color)*2);
    }      
}

//Send the reduced border of a spincolor vector
void communicate_lx_redspincolor_borders(spincolor *s,redspincolor *temp_out,redspincolor *temp_in)
{
  int oall=0,iall=0;

  if(temp_out==NULL){oall=1; temp_out=(redspincolor*)malloc(loc_bord*sizeof(redspincolor));}
  if(temp_in==NULL) {iall=1; temp_in=(redspincolor*)malloc(loc_bord*sizeof(redspincolor));}

  apply_lxspincolor_border_projector(temp_out,s);

  int nrequest=0;
  MPI_Request request[16];
  MPI_Status status[16];

  for(int i=0;i<4;i++)
    if(paral_dir[i]!=0)
      {
	int pos_up=(start_lx_bord_rece_up[i]-loc_vol);
	int pos_dw=(start_lx_bord_rece_dw[i]-loc_vol);
	
       //sending the upper border to the lower node
       MPI_Isend((void*)(temp_out+pos_up),1,MPI_LXREDSPINCOLOR_BORD[i],rank_neighdw[i],83+i,cart_comm,&request[nrequest++]);
       MPI_Irecv((void*)(temp_in+pos_up),1,MPI_LXREDSPINCOLOR_BORD[i],rank_neighup[i],83+i,cart_comm,&request[nrequest++]);
	
       //sending the lower border to the upper node
       MPI_Isend((void*)(temp_out+pos_dw),1,MPI_LXREDSPINCOLOR_BORD[i],rank_neighup[i],87+i,cart_comm,&request[nrequest++]);
       MPI_Irecv((void*)(temp_in+pos_dw),1,MPI_LXREDSPINCOLOR_BORD[i],rank_neighdw[i],87+i,cart_comm,&request[nrequest++]);
      }
  
  if(nrequest>0) MPI_Waitall(nrequest,request,status);

  expand_projected_lxspincolor_border(s,temp_in);

  if(iall) free(temp_in);
  if(oall) free(temp_out);
}
