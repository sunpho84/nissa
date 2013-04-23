#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <mpi.h>

#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../geometry/geometry_lx.h"
#include "../new_types/new_types_definitions.h"
#include "../routines/ios.h"
#include "../routines/thread.h"

#include "buffered_borders.h"

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
void start_communicating_lx_borders(int *nrequest,MPI_Request *request,char *data,MPI_Datatype *MPI_BORDS_SEND,MPI_Datatype *MPI_BORDS_RECE,int nbytes_per_site)
{
  if(!check_borders_valid(data))
    {
      GET_THREAD_ID();
      
      if(!IS_MASTER_THREAD) (*nrequest)=4*nparal_dir;
      else
	{
	  (*nrequest)=0;
	  crash_if_borders_not_allocated(data);  
  
	  int imessage=654325;
	  
	  verbosity_lv3_master_printf("Start communicating borders of %s\n",get_vec_name((void*)data));
	  for(int mu=0;mu<4;mu++)
	    if(paral_dir[mu]!=0)
	      {
		//sending the upper border to the lower node
		MPI_Irecv(data+start_lx_bord_rece_up[mu]*nbytes_per_site,1,MPI_BORDS_RECE[mu],rank_neighup[mu],imessage,cart_comm,request+((*nrequest)++));
		MPI_Isend(data+start_lx_bord_send_up[mu]*nbytes_per_site,1,MPI_BORDS_SEND[mu],rank_neighdw[mu],imessage++,cart_comm,request+((*nrequest)++));
		
		//sending the lower border to the upper node
		MPI_Irecv(data+start_lx_bord_rece_dw[mu]*nbytes_per_site,1,MPI_BORDS_RECE[mu],rank_neighdw[mu],imessage,cart_comm,request+((*nrequest)++));
		MPI_Isend(data+start_lx_bord_send_dw[mu]*nbytes_per_site,1,MPI_BORDS_SEND[mu],rank_neighup[mu],imessage++,cart_comm,request+((*nrequest)++));
	      }
	}
    }
  else (*nrequest)=0;
}

//wait to finish communications
void finish_communicating_lx_borders(int *nrequest,MPI_Request *request,char *data)
{
  if((*nrequest)>0)
    {
      GET_THREAD_ID();
      
      if(IS_MASTER_THREAD)
	{
	  tot_nissa_comm_time-=take_time();
	  verbosity_lv3_master_printf("Waiting to finish %d communication of borders of vector %s\n",*nrequest,get_vec_name(data));
	  MPI_Status status[*nrequest];
	  MPI_Waitall(*nrequest,request,status);
	  tot_nissa_comm_time+=take_time();
	}
      
      set_borders_valid(data);
      (*nrequest)=0;
      set_edges_invalid(data);
    }
}

//merge start and finish
void communicate_lx_borders(char *data,MPI_Datatype *MPI_BORDS_SEND,MPI_Datatype *MPI_BORDS_RECE,int nbytes_per_site)
{
  MPI_Request request[16];
  int nrequest;
  start_communicating_lx_borders(&nrequest,request,data,MPI_BORDS_SEND,MPI_BORDS_RECE,nbytes_per_site);
  finish_communicating_lx_borders(&nrequest,request,data);
}

//Useful for gauge fixing and hyp
void communicate_lx_su3_borders(su3 *u)
{communicate_lx_borders((char*)u,MPI_LX_SU3_BORDS_SEND,MPI_LX_SU3_BORDS_RECE,sizeof(su3));}

//Send the borders of the gauge configuration
void communicate_lx_quad_su3_borders(quad_su3 *conf)
{
#ifdef BGQ
  buffered_communicate_lx_borders(conf,&buffered_lx_quad_su3_comm);
#else
  communicate_lx_borders((char*)conf,MPI_LX_QUAD_SU3_BORDS_SEND,MPI_LX_QUAD_SU3_BORDS_RECE,sizeof(quad_su3));
#endif
}

//Send the borders of a spincolor vector
void communicate_lx_spincolor_borders(spincolor *s)
{
#ifdef BGQ
  buffered_communicate_lx_borders(s,&buffered_lx_spincolor_comm);
#else
  communicate_lx_borders((char*)s,MPI_LX_SPINCOLOR_BORDS_SEND,MPI_LX_SPINCOLOR_BORDS_RECE,sizeof(spincolor));
#endif
}

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
void start_communicating_lx_spincolor_borders(int *nrequest,MPI_Request *request,spincolor *s)
{
#ifdef BGQ
  buffered_start_communicating_lx_borders(&buffered_lx_spincolor_comm,s);
#else
 start_communicating_lx_borders(nrequest,request,(char*)s,MPI_LX_SPINCOLOR_BORDS_SEND,MPI_LX_SPINCOLOR_BORDS_RECE,sizeof(spincolor));
#endif
}
void finish_communicating_lx_spincolor_borders(int *nrequest,MPI_Request *request,spincolor *s)
{
#ifdef BGQ
  buffered_finish_communicating_lx_borders(s,&buffered_lx_spincolor_comm);
#else
  finish_communicating_lx_borders(nrequest,request,(char*)s);
#endif
}

void finish_communicating_ev_borders(int *nrequest,MPI_Request *request,char *ev_data)
{
  if((*nrequest)>0)
    {
      GET_THREAD_ID();
      
      if(IS_MASTER_THREAD)
	{
	  tot_nissa_comm_time-=take_time();
	  verbosity_lv3_master_printf("Waiting to finish %d communication of ev borders of vector %s\n",*nrequest,get_vec_name(ev_data));
	  MPI_Status status[*nrequest];
	  MPI_Waitall(*nrequest,request,status);
	  tot_nissa_comm_time+=take_time();
	}
      
      set_borders_valid(ev_data);
      (*nrequest)=0;
      set_edges_invalid(ev_data);
    }
}

//Send the borders of the data
void start_communicating_ev_borders(int *nrequest,MPI_Request *request,char *ev_data,MPI_Datatype *MPI_EV_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_EV_BORDS_RECE,int nbytes_per_site)
{
  if(!check_borders_valid(ev_data))
    {
      GET_THREAD_ID();
      
      if(!IS_MASTER_THREAD) (*nrequest)=4*nparal_dir;
      else
	{
	  (*nrequest)=0;
	  crash_if_borders_not_allocated(ev_data);
	  
	  int imessage=534245;
	  
	  for(int mu=0;mu<3;mu++)
	    if(paral_dir[mu]!=0)
	      {
		//sending the upper border to the lower node
		MPI_Irecv(ev_data+start_eo_bord_rece_up[mu]*nbytes_per_site,1,MPI_EV_BORDS_RECE[mu],rank_neighup[mu],imessage,cart_comm,request+((*nrequest)++));
		MPI_Isend(ev_data+start_eo_bord_send_up[mu]*nbytes_per_site,1,MPI_EV_BORDS_SEND_TXY[mu],rank_neighdw[mu],imessage++,cart_comm,request+((*nrequest)++));
		
		//sending the lower border to the upper node
		MPI_Irecv(ev_data+start_eo_bord_rece_dw[mu]*nbytes_per_site,1,MPI_EV_BORDS_RECE[mu],rank_neighdw[mu],imessage,cart_comm,request+((*nrequest)++));
		MPI_Isend(ev_data+start_eo_bord_send_dw[mu]*nbytes_per_site,1,MPI_EV_BORDS_SEND_TXY[mu],rank_neighup[mu],imessage++,cart_comm,request+((*nrequest)++));
	      }
	  
	  if(paral_dir[3]!=0)
	    {
	      //sending the upper border to the lower node
	      MPI_Irecv(ev_data+start_eo_bord_rece_up[3]*nbytes_per_site,1,MPI_EV_BORDS_RECE[3],rank_neighup[3],imessage,cart_comm,request+((*nrequest)++));
	      MPI_Isend(ev_data,1,MPI_EV_BORDS_SEND_Z[0],rank_neighdw[3],imessage++,cart_comm,request+((*nrequest)++));
	      
	      //sending the lower border to the upper node
	      MPI_Irecv(ev_data+start_eo_bord_rece_dw[3]*nbytes_per_site,1,MPI_EV_BORDS_RECE[3],rank_neighdw[3],imessage,cart_comm,request+((*nrequest)++));
	      MPI_Isend(ev_data,1,MPI_EV_BORDS_SEND_Z[1],rank_neighup[3],imessage++,cart_comm,request+((*nrequest)++));
	    }
	  
	  if((*nrequest)>0) verbosity_lv3_master_printf("Starting communication of ev borders of vector %s\n",get_vec_name(ev_data));
	}
    }
  else (*nrequest)=0;
}

//Send the borders of the data
void communicate_ev_borders(char *ev_data,MPI_Datatype *MPI_EV_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_EV_BORDS_RECE,int nbytes_per_site)
{
  MPI_Request request[16];
  int nrequest;
  start_communicating_ev_borders(&nrequest,request,ev_data,MPI_EV_BORDS_SEND_TXY,MPI_EV_BORDS_SEND_Z,MPI_EV_BORDS_RECE,nbytes_per_site);
  finish_communicating_ev_borders(&nrequest,request,ev_data);
}

//Send the borders of the data
void communicate_od_borders(char *od_data,MPI_Datatype *MPI_EV_BORDS_SEND_TXY,MPI_Datatype *MPI_OD_BORDS_SEND_Z,MPI_Datatype *MPI_EV_BORDS_RECE,int nbytes_per_site)
{
  if(!check_borders_valid(od_data))
    {
      crash_if_borders_not_allocated(od_data);
      
      GET_THREAD_ID();
      
      if(IS_MASTER_THREAD)
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
	  tot_nissa_comm_time+=take_time();
	}
      
      set_borders_valid(od_data);
      set_edges_invalid(od_data);
    }
}

//Send the borders of the data
void communicate_eo_borders(char **data,MPI_Datatype *MPI_EO_BORDS_SEND_TXY,MPI_Datatype *MPI_EV_BORDS_SEND_Z,MPI_Datatype *MPI_OD_BORDS_SEND_Z,MPI_Datatype *MPI_EO_BORDS_RECE,int nbytes_per_site,int OPT=SEND_BACKWARD_BORD|SEND_FORWARD_BORD)
{
  if(!check_borders_valid(data[EVN])||!check_borders_valid(data[ODD]))
    {
      GET_THREAD_ID();
      
      if(IS_MASTER_THREAD)
	{
	  int nrequest=0;
	  MPI_Request request[32];
	  int imessage=63458729;
	  for(int par=0;par<2;par++)
	    {
	      crash_if_borders_not_allocated(data[par]);
	      
	      for(int mu=0;mu<3;mu++) //0,1,2
		if(paral_dir[mu]!=0)
		  {
		    //sending the upper border to the lower node
		    if(OPT & SEND_BACKWARD_BORD)
		      {
			imessage++;
			MPI_Irecv(data[par]+start_eo_bord_rece_up[mu]*nbytes_per_site,1,MPI_EO_BORDS_RECE[mu],rank_neighup[mu],imessage,
				  cart_comm,request+(nrequest++));
			MPI_Isend(data[par]+start_eo_bord_send_up[mu]*nbytes_per_site,1,MPI_EO_BORDS_SEND_TXY[mu],rank_neighdw[mu],imessage,
				  cart_comm,request+(nrequest++));
		      }
		    
		    //sending the lower border to the upper node
		    if(OPT & SEND_FORWARD_BORD)
		      {
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
	      //sending the upper border to the lower node
	      if(OPT & SEND_BACKWARD_BORD)
		{
		  imessage++;
		  MPI_Irecv(data[EVN]+start_eo_bord_rece_up[3]*nbytes_per_site,1,MPI_EO_BORDS_RECE[3],rank_neighup[3],imessage,cart_comm,request+(nrequest++));
		  MPI_Isend(data[EVN],1,MPI_EV_BORDS_SEND_Z[0],rank_neighdw[3],imessage,cart_comm,request+(nrequest++));
		}
	      
	      //sending the lower border to the upper node
	      if(OPT & SEND_FORWARD_BORD)
		{
		  imessage++;
		  MPI_Irecv(data[EVN]+start_eo_bord_rece_dw[3]*nbytes_per_site,1,MPI_EO_BORDS_RECE[3],rank_neighdw[3],imessage,cart_comm,request+(nrequest++));
		  MPI_Isend(data[EVN],1,MPI_EV_BORDS_SEND_Z[1],rank_neighup[3],imessage,cart_comm,request+(nrequest++));
		}
	      
	      //sending the upper border to the lower node
	      if(OPT & SEND_BACKWARD_BORD)
		{
		  imessage++;
		  MPI_Irecv(data[ODD]+start_eo_bord_rece_up[3]*nbytes_per_site,1,MPI_EO_BORDS_RECE[3],rank_neighup[3],imessage,cart_comm,request+(nrequest++));
		  MPI_Isend(data[ODD],1,MPI_OD_BORDS_SEND_Z[0],rank_neighdw[3],imessage,cart_comm,request+(nrequest++));
		}
	      
	      //sending the lower border to the upper node
	      if(OPT & SEND_FORWARD_BORD)
		{
		  imessage++;
		  MPI_Irecv(data[ODD]+start_eo_bord_rece_dw[3]*nbytes_per_site,1,MPI_EO_BORDS_RECE[3],rank_neighdw[3],imessage,cart_comm,request+(nrequest++));
		  MPI_Isend(data[ODD],1,MPI_OD_BORDS_SEND_Z[1],rank_neighup[3],imessage,cart_comm,request+(nrequest++));
		}
	    }
	
	  tot_nissa_comm_time-=take_time();
	  if(nrequest>0)
	    {
	      if(nrequest>0)
		{
		  verbosity_lv3_master_printf("Communication of evn & odd borders of vector %s\n",get_vec_name(data[0]));
		  MPI_Waitall(nrequest,request,MPI_STATUS_IGNORE);
		  verbosity_lv3_master_printf("Communication finished\n",get_vec_name(data[0]));
		}
	      nrequest=0;
	    }
	  tot_nissa_comm_time+=take_time();
	}
      
      verbosity_lv3_master_printf("Setting EVN borders valid (thread_pool_locked: %d)\n",thread_pool_locked);
      set_borders_valid(data[EVN]);
      verbosity_lv3_master_printf("Setting ODD borders valid\n");
      set_borders_valid(data[ODD]);
    }
  else
    verbosity_lv3_master_printf("Borders valid\n");
}

//Send the borders of the gauge configuration
void communicate_eo_quad_su3_borders(quad_su3 **eo_conf,int OPT=SEND_BACKWARD_BORD|SEND_FORWARD_BORD)
{
#ifdef BGQ
  buffered_communicate_ev_and_od_borders((void**)eo_conf,&buffered_lx_quad_su3_comm);
#else
  communicate_eo_borders((char**)eo_conf,MPI_EO_QUAD_SU3_BORDS_SEND_TXY,MPI_EV_QUAD_SU3_BORDS_SEND_Z,MPI_OD_QUAD_SU3_BORDS_SEND_Z,MPI_EO_QUAD_SU3_BORDS_RECE,sizeof(quad_su3),OPT);
#endif
}

//Send the borders of an even color
void communicate_ev_color_borders(color *ev)
{
#ifdef BGQ
  buffered_communicate_ev_or_od_borders(ev,&buffered_eo_color_comm,EVN);
#else
  communicate_ev_borders((char*)ev,MPI_EO_COLOR_BORDS_SEND_TXY,MPI_EV_COLOR_BORDS_SEND_Z,MPI_EO_COLOR_BORDS_RECE,sizeof(color));
#endif
}
void communicate_od_color_borders(color *od)
{
#ifdef BGQ
  buffered_communicate_ev_or_od_borders(od,&buffered_eo_color_comm,ODD);
#else
  communicate_od_borders((char*)od,MPI_EO_COLOR_BORDS_SEND_TXY,MPI_OD_COLOR_BORDS_SEND_Z,MPI_EO_COLOR_BORDS_RECE,sizeof(color));
#endif
}
void communicate_eo_color_borders(color **eos)
{communicate_eo_borders((char**)eos,MPI_EO_COLOR_BORDS_SEND_TXY,MPI_EV_COLOR_BORDS_SEND_Z,MPI_OD_COLOR_BORDS_SEND_Z,MPI_EO_COLOR_BORDS_RECE,sizeof(color));}
void start_communicating_ev_color_borders(int *nrequest,MPI_Request *request,color *ev)
{
#ifdef BGQ
  buffered_start_communicating_ev_or_od_borders(&buffered_eo_color_comm,ev,EVN);
#else
  start_communicating_ev_borders(nrequest,request,(char*)ev,MPI_EO_COLOR_BORDS_SEND_TXY,MPI_EV_COLOR_BORDS_SEND_Z,MPI_EO_COLOR_BORDS_RECE,sizeof(color));
#endif
}
void finish_communicating_ev_color_borders(int *nrequest,MPI_Request *request,color *ev)
{
#ifdef BGQ
  buffered_finish_communicating_ev_or_od_borders(ev,&buffered_eo_color_comm);
#else
  finish_communicating_ev_borders(nrequest,request,(char*)ev);
#endif
}

//Send the borders of an even spincolor
void communicate_ev_spincolor_borders(spincolor *ev)
{communicate_ev_borders((char*)ev,MPI_EO_SPINCOLOR_BORDS_SEND_TXY,MPI_EV_SPINCOLOR_BORDS_SEND_Z,MPI_EO_SPINCOLOR_BORDS_RECE,sizeof(spincolor));}
void communicate_od_spincolor_borders(spincolor *od)
{communicate_od_borders((char*)od,MPI_EO_SPINCOLOR_BORDS_SEND_TXY,MPI_OD_SPINCOLOR_BORDS_SEND_Z,MPI_EO_SPINCOLOR_BORDS_RECE,sizeof(spincolor));}
void start_communicating_ev_spincolor_borders(int *nrequest,MPI_Request *request,spincolor *ev)
{start_communicating_ev_borders(nrequest,request,(char*)ev,MPI_EO_SPINCOLOR_BORDS_SEND_TXY,MPI_EV_SPINCOLOR_BORDS_SEND_Z,MPI_EO_SPINCOLOR_BORDS_RECE,sizeof(spincolor));}
void finish_communicating_ev_spincolor_borders(int *nrequest,MPI_Request *request,spincolor *ev)
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

