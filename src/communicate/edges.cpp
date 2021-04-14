#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>

#include "communicate.hpp"

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
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
    const LocLxSite min_size=nbytes_per_site*(edge_vol()+bordVol()+locVol());
    crash_if_edges_not_allocated(data,min_size.nastyConvert());
    communicate_lx_borders(data,bord_comm);
    
    if(not check_edges_valid(data))
      {
	if(IS_MASTER_THREAD)
	  {
	    int nrequest=0;
	    MPI_Request request[NDIM*(NDIM-1)*4];
	    MPI_Status status[NDIM*(NDIM-1)*4];
	    LocLxSite send,rece;
	    int imessage=0;
	    LocCoords x;
	    FOR_ALL_DIRECTIONS(mu)
	      x(mu)=0;
	    
	    FOR_ALL_DIRECTIONS(idir)
	      for(Direction jdir=idir+1;jdir<NDIM;jdir++)
		if(paral_dir(idir) and paral_dir(jdir))
		  {
		    int iedge=edge_numb[idir.nastyConvert()][jdir.nastyConvert()];
		    LocLxSite pos_edge_offset;
		    
		    //take the starting point of the border
		    x(jdir)=locSize(jdir)-1;
		    pos_edge_offset=spatLxOfProjectedCoords<LocLxSite>(x,idir);
		    x(jdir)=0;
		    
		    //Send the i-j- internal edge to the j- rank as i-j+ external edge
		    send=(locVol.nastyConvert()+bord_offset[idir.nastyConvert()].nastyConvert())*nbytes_per_site;
		    rece=(locVol.nastyConvert()+bordVol.nastyConvert()+edge_offset[iedge].nastyConvert()+edge_vol.nastyConvert()/4)*nbytes_per_site;
		    MPI_Irecv(data+rece.nastyConvert(),1,MPI_EDGES_RECE[iedge],rank_neighup(jdir)(),imessage,cart_comm,request+(nrequest++));
		    MPI_Isend(data+send.nastyConvert(),1,MPI_EDGES_SEND[iedge],rank_neighdw(jdir)(),imessage++,cart_comm,request+(nrequest++));
		    
		    //Send the i-j+ internal edge to the j+ rank as i-j- external edge
		    send=(locVol.nastyConvert()+bord_offset[idir.nastyConvert()].nastyConvert()+pos_edge_offset)*nbytes_per_site;
		    rece=(locVol.nastyConvert()+bordVol.nastyConvert()+edge_offset[iedge])*nbytes_per_site;
		    MPI_Irecv(data+rece.nastyConvert(),1,MPI_EDGES_RECE[iedge],rank_neighdw(jdir)(),imessage,cart_comm,request+(nrequest++));
		    MPI_Isend(data+send.nastyConvert(),1,MPI_EDGES_SEND[iedge],rank_neighup(jdir)(),imessage++,cart_comm,request+(nrequest++));
		    
		    //Send the i+j- internal edge to the j- rank as i+j+ external edge
		    send=(locVol.nastyConvert()+bord_offset[idir.nastyConvert()].nastyConvert()+bordVol.nastyConvert()/2)*nbytes_per_site;
		    rece=(locVol+bordVol.nastyConvert()+edge_offset[iedge].nastyConvert()+3*edge_vol.nastyConvert()/4)*nbytes_per_site;
		    MPI_Irecv(data+rece.nastyConvert(),1,MPI_EDGES_RECE[iedge],rank_neighup(jdir)(),imessage,cart_comm,request+(nrequest++));
		    MPI_Isend(data+send.nastyConvert(),1,MPI_EDGES_SEND[iedge],rank_neighdw(jdir)(),imessage++,cart_comm,request+(nrequest++));
		    
		    //Send the i+j+ internal edge to the j+ rank as i+j- external edge
		    send=(locVol.nastyConvert()+bord_offset[idir.nastyConvert()].nastyConvert()+bordVol.nastyConvert()/2+pos_edge_offset)*nbytes_per_site;
		    rece=(locVol+bordVol.nastyConvert()+edge_offset[iedge].nastyConvert()+edge_vol.nastyConvert()/2)*nbytes_per_site;
		    MPI_Irecv(data+rece.nastyConvert(),1,MPI_EDGES_RECE[iedge],rank_neighdw(jdir)(),imessage,cart_comm,request+(nrequest++));
		    MPI_Isend(data+send.nastyConvert(),1,MPI_EDGES_SEND[iedge],rank_neighup(jdir)(),imessage++,cart_comm,request+(nrequest++));
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
  void communicate_eo_edges(eo_ptr<void> data,comm_t &bord_comm,MPI_Datatype *MPI_EDGES_SEND,MPI_Datatype *MPI_EDGES_RECE,int nbytes_per_site)
  {
    if(!check_edges_valid(data[EVN])||!check_edges_valid(data[ODD]))
      {
	communicate_ev_and_od_borders(data,bord_comm);
	
	
	if(IS_MASTER_THREAD)
	  {
	    int nrequest=0;
	    int nrequest_tot=0;
	    FOR_ALL_DIRECTIONS(mu)
	      for(Direction nu=mu+1;nu<NDIM;nu++)
		if(paral_dir(mu) and paral_dir(nu))
		  nrequest_tot+=16;
	    
	    MPI_Request request[nrequest_tot];
	    MPI_Status status[nrequest_tot];
	    int imessage=0;
	    
	    verbosity_lv3_master_printf("Communicating edges of %s\n",get_vect_name(data[EVN]));
	    
	    for(int par=0;par<2;par++)
	      {
		//check_edges_allocated(data);
		int64_t min_size=nbytes_per_site*(edge_volh+bord_volh+locVolh).nastyConvert();
		crash_if_edges_not_allocated(data[par],min_size);
		
		//"v" refer to the verse of the dir
		for(int vmu=0;vmu<2;vmu++)
		  FOR_ALL_DIRECTIONS(mu)
		    for(int vnu=0;vnu<2;vnu++)
		      for(Direction nu=mu+1;nu<NDIM;nu++)
			if(paral_dir(mu) and paral_dir(nu))
			  {
			    const int iedge=edge_numb[mu.nastyConvert()][nu.nastyConvert()];
			    
			    //communicators verse refers to the internal edge
			    int icomm_send=((par*2+vmu)*2+vnu)*6+iedge;
			    int icomm_recv=iedge;
			    
			    //Send the (mu,nu) internal edge to the nu rank as (mu,-nu) external edge
			    //Receive the (mu,-nu) external edge from the -nu rank (mu,nu) internal edge
			    LocLxSite ext_edge_recv_start=(locVolh+bord_volh+edge_volh/4*(vmu*2+!vnu)+edge_offset[iedge].nastyConvert()/2).nastyConvert()*nbytes_per_site;
			    int int_edge_send_start=locVolh.nastyConvert()*nbytes_per_site;
			    MPI_Irecv((char*)(data[par])+ext_edge_recv_start.nastyConvert(),1,MPI_EDGES_RECE[icomm_recv],rank_neigh[!vnu](nu)(),imessage,cart_comm,request+(nrequest++));
			    MPI_Isend((char*)(data[par])+int_edge_send_start,1,MPI_EDGES_SEND[icomm_send],rank_neigh[vnu](nu)(),imessage,cart_comm,request+(nrequest++));
			    
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
  void communicate_eo_quad_su3_edges(eo_ptr<quad_su3> conf)
  {communicate_eo_edges({conf[EVN],conf[ODD]},lx_quad_su3_comm,MPI_EO_QUAD_SU3_EDGES_SEND,MPI_EO_QUAD_SU3_EDGES_RECE,sizeof(quad_su3));}
  
  void communicate_eo_as2t_su3_edges(eo_ptr<as2t_su3> a)
  {communicate_eo_edges({a[EVN],a[ODD]},lx_as2t_su3_comm,MPI_EO_AS2T_SU3_EDGES_SEND,MPI_EO_AS2T_SU3_EDGES_RECE,sizeof(as2t_su3));}
}
