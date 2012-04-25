#include "../base/global_variables.h"
#include "../base/debug.h"
#include "../base/vectors.h"
#include "../base/routines.h"
#include "../new_types/new_types_definitions.h"
#include "../new_types/su3.h"

//set the eo geometry
void set_eo_geometry()
{
  if(nissa_eo_geom_inited) crash("E/O Geometry already initialized!");
  
  //check that all local sizes are multiples of 2
  int ok=1;
  for(int mu=0;mu<4;mu++) ok&=(loc_size[mu]%2==0);
  if(!ok) crash("local lattice size odd!");
  
  //set half the vol, bord and edge size
  glb_volh=glb_vol/2;
  loc_volh=loc_vol/2;
  bord_volh=bord_vol/2;
  edge_volh=edge_vol/2;
  
  //set the various time-slice types
  loclx_parity=nissa_malloc("loclx_parity",loc_vol+bord_vol+edge_vol,int);
  ignore_borders_communications_warning(loclx_parity);
  
  loceo_of_loclx=nissa_malloc("loceo_of_loclx",loc_vol+bord_vol+edge_vol,int);
  ignore_borders_communications_warning(loceo_of_loclx);
  
  for(int par=0;par<2;par++) loclx_of_loceo[par]=nissa_malloc("loclx_of_loceo",loc_volh+bord_volh+edge_volh,int);
  for(int par=0;par<2;par++) loceo_neighup[par]=nissa_malloc("loceo_neighup",loc_volh+bord_volh+edge_volh,coords);
  for(int par=0;par<2;par++) loceo_neighdw[par]=nissa_malloc("loceo_neighdw",loc_volh+bord_volh+edge_volh,coords);
  for(int par=0;par<2;par++) ignore_borders_communications_warning(loclx_of_loceo[par]);
  for(int par=0;par<2;par++) ignore_borders_communications_warning(loceo_neighup[par]);
  for(int par=0;par<2;par++) ignore_borders_communications_warning(loceo_neighdw[par]);

  //Label the bulk sites
  int iloc_eo[2]={0,0};
  for(int loclx=0;loclx<loc_vol+bord_vol+edge_vol;loclx++)
    {
      //calculate global coord and parity
      int par=0;
      for(int mu=0;mu<4;mu++)
	par+=glb_coord_of_loclx[loclx][mu];
      
      par%=2;
      
      //fix parity of local index
      loclx_parity[loclx]=par;
      
      //associate the e/o index to lx sites and vice-versa
      loceo_of_loclx[loclx]=iloc_eo[par];
      loclx_of_loceo[par][iloc_eo[par]]=loclx;
      iloc_eo[par]++;
    }
  
  //Fix the movements among e/o ordered sites
  for(int loclx=0;loclx<loc_vol+bord_vol;loclx++)
    for(int mu=0;mu<4;mu++)
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
	else
	  loceo_neighdw[par][loceo][mu]=loclx_dw;
      }
  
  //init sender and receiver points for borders
  for(int mu=0;mu<4;mu++)
    if(paral_dir[mu]!=0)
      {
	start_eo_bord_send_up[mu]=loceo_of_loclx[start_lx_bord_send_up[mu]];
	start_eo_bord_rece_up[mu]=loceo_of_loclx[start_lx_bord_rece_up[mu]];
	start_eo_bord_send_dw[mu]=loceo_of_loclx[start_lx_bord_send_dw[mu]];
	start_eo_bord_rece_dw[mu]=loceo_of_loclx[start_lx_bord_rece_dw[mu]];
      }
  
  master_printf("E/O Geometry intialized\n");
  
  nissa_eo_geom_inited=1;
}

//definitions of e/o split sender for borders
void initialize_eo_bord_senders_of_kind(MPI_Datatype *MPI_EO_BORD_SEND_TXY,MPI_Datatype *MPI_EV_BORD_SEND_Z,MPI_Datatype *MPI_OD_BORD_SEND_Z,MPI_Datatype *base)
{
  //Various type useful for edges and sub-borders
  MPI_Datatype MPI_EO_3_SLICE;
  MPI_Datatype MPI_EO_23_SLICE;
  MPI_Type_contiguous(loc_size[3]/2,*base,&MPI_EO_3_SLICE);
  MPI_Type_contiguous(loc_size[2]*loc_size[3]/2,*base,&MPI_EO_23_SLICE);
  
  ///////////define the sender for the 4 kinds of borders////////////
  MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3]/2,*base,&(MPI_EO_BORD_SEND_TXY[0]));
  MPI_Type_vector(loc_size[0],1,loc_size[1],MPI_EO_23_SLICE,&(MPI_EO_BORD_SEND_TXY[1]));
  MPI_Type_vector(loc_size[0]*loc_size[1],1,loc_size[2],MPI_EO_3_SLICE,&(MPI_EO_BORD_SEND_TXY[2]));
  //Commit
  for(int ibord=0;ibord<3;ibord++) MPI_Type_commit(&(MPI_EO_BORD_SEND_TXY[ibord]));
  
  //the z sending border is a mess
  int eo_bord_z_size=loc_volh/loc_size[3];
  int *ev_bord_z_pos_disp_dw=nissa_malloc("ev_bord_z_disp_dw",eo_bord_z_size,int);
  int *ev_bord_z_pos_disp_up=nissa_malloc("ev_bord_z_disp_up",eo_bord_z_size,int);
  int *od_bord_z_pos_disp_dw=nissa_malloc("od_bord_z_disp_dw",eo_bord_z_size,int);
  int *od_bord_z_pos_disp_up=nissa_malloc("od_bord_z_disp_up",eo_bord_z_size,int);
  int *single=nissa_malloc("single",eo_bord_z_size,int);
  int ev_izdw=0,ev_izup=0;
  int od_izdw=0,od_izup=0;
  nissa_loc_volh_loop(ieo)
    {
      int ev_ilx=loclx_of_loceo[0][ieo];
      int od_ilx=loclx_of_loceo[1][ieo];
      int ev_x3=loc_coord_of_loclx[ev_ilx][3];
      int od_x3=loc_coord_of_loclx[od_ilx][3];
      if(ev_x3==0) ev_bord_z_pos_disp_dw[ev_izdw++]=ieo;
      if(ev_x3==loc_size[3]-1) ev_bord_z_pos_disp_up[ev_izup++]=ieo;
      if(od_x3==0) od_bord_z_pos_disp_dw[od_izdw++]=ieo;
      if(od_x3==loc_size[3]-1) od_bord_z_pos_disp_up[od_izup++]=ieo;
    }
  for(int ibord_z=0;ibord_z<eo_bord_z_size;ibord_z++)
    single[ibord_z]=1;
  
  MPI_Type_indexed(eo_bord_z_size,single,ev_bord_z_pos_disp_dw,*base,&(MPI_EV_BORD_SEND_Z[0]));
  MPI_Type_indexed(eo_bord_z_size,single,ev_bord_z_pos_disp_up,*base,&(MPI_EV_BORD_SEND_Z[1]));
  
  MPI_Type_indexed(eo_bord_z_size,single,od_bord_z_pos_disp_dw,*base,&(MPI_OD_BORD_SEND_Z[0]));
  MPI_Type_indexed(eo_bord_z_size,single,od_bord_z_pos_disp_up,*base,&(MPI_OD_BORD_SEND_Z[1]));
  
  //commit the mess
  MPI_Type_commit(&(MPI_EV_BORD_SEND_Z[0]));
  MPI_Type_commit(&(MPI_EV_BORD_SEND_Z[1]));
  
  MPI_Type_commit(&(MPI_OD_BORD_SEND_Z[0]));
  MPI_Type_commit(&(MPI_OD_BORD_SEND_Z[1]));
  
  nissa_free(single);
  nissa_free(ev_bord_z_pos_disp_dw);
  nissa_free(ev_bord_z_pos_disp_up);
  nissa_free(od_bord_z_pos_disp_dw);
  nissa_free(od_bord_z_pos_disp_up);
}

//definitions of e/o split receivers for borders
void initialize_eo_bord_receivers_of_kind(MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base)
{
  //define the 4 dir borders receivers, which are contiguous in memory
  MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3]/2,*base,&(MPI_EO_BORD_RECE[0]));
  MPI_Type_contiguous(loc_size[0]*loc_size[2]*loc_size[3]/2,*base,&(MPI_EO_BORD_RECE[1]));
  MPI_Type_contiguous(loc_size[0]*loc_size[1]*loc_size[3]/2,*base,&(MPI_EO_BORD_RECE[2]));
  MPI_Type_contiguous(loc_size[0]*loc_size[1]*loc_size[2]/2,*base,&(MPI_EO_BORD_RECE[3]));
  for(int ibord=0;ibord<4;ibord++) MPI_Type_commit(&(MPI_EO_BORD_RECE[ibord]));
}

//initalize senders and receivers for borders of e/o split ordered vectors
void set_eo_bord_senders_and_receivers(MPI_Datatype *MPI_EO_BORD_SEND_TXY,MPI_Datatype *MPI_EV_BORD_SEND_Z,MPI_Datatype *MPI_OD_BORD_SEND_Z,MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base)
{
  initialize_eo_bord_senders_of_kind(MPI_EO_BORD_SEND_TXY,MPI_EV_BORD_SEND_Z,MPI_OD_BORD_SEND_Z,base);
  initialize_eo_bord_receivers_of_kind(MPI_EO_BORD_RECE,base);
}

//definitions of e/o split sender edges
void initialize_eo_edge_senders_of_kind(MPI_Datatype *MPI_EO_EDGES_SEND,MPI_Datatype *base)
{
  for(int par=0;par<2;par++)
    for(int vmu=0;vmu<2;vmu++)
      for(int mu=0;mu<4;mu++)
	for(int vnu=0;vnu<2;vnu++)
	  for(int nu=mu+1;nu<4;nu++)
	    if(paral_dir[mu] && paral_dir[nu])
	      {
		int iedge=edge_numb[mu][nu];
		int icomm=((par*2+vmu)*2+vnu)*6+iedge;
		
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
void initialize_eo_edge_receivers_of_kind(MPI_Datatype *MPI_EDGES_RECE,MPI_Datatype *base)
{
  //define the 6 edges receivers, which are contiguous in memory
  MPI_Type_contiguous(loc_size[2]*loc_size[3]/2,*base,&(MPI_EDGES_RECE[0]));
  MPI_Type_contiguous(loc_size[1]*loc_size[3]/2,*base,&(MPI_EDGES_RECE[1]));
  MPI_Type_contiguous(loc_size[1]*loc_size[2]/2,*base,&(MPI_EDGES_RECE[2]));
  MPI_Type_contiguous(loc_size[0]*loc_size[3]/2,*base,&(MPI_EDGES_RECE[3]));
  MPI_Type_contiguous(loc_size[0]*loc_size[2]/2,*base,&(MPI_EDGES_RECE[4]));
  MPI_Type_contiguous(loc_size[0]*loc_size[1]/2,*base,&(MPI_EDGES_RECE[5]));
  for(int iedge=0;iedge<6;iedge++) MPI_Type_commit(&(MPI_EDGES_RECE[iedge]));
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
  if(!nissa_eo_geom_inited)
    crash("asking to unset never initialized E/O Geometry!");

  master_printf("Unsetting E/O Geometry\n");
  
  for(int par=0;par<2;par++)
    {
      nissa_free(loclx_of_loceo[par]);
      nissa_free(loceo_neighup[par]);
      nissa_free(loceo_neighdw[par]);
    }
  nissa_free(loclx_parity);
  nissa_free(loceo_of_loclx);
  
  nissa_eo_geom_inited=0;
}

//add or remove staggered phases to/from a conf
void addrem_stagphases_to_eo_conf(quad_su3 **eo_conf)
{
  if(!nissa_eo_geom_inited) set_eo_geometry();
  
  //work also on borders and edges if allocated and valid
  int ending=loc_volh;
  if(check_borders_allocated(eo_conf[0]) && check_borders_allocated(eo_conf[1]) && check_borders_valid(eo_conf[0]) && check_borders_valid(eo_conf[1])) ending+=bord_volh;
  if(check_edges_allocated(eo_conf[0])   && check_edges_allocated(eo_conf[1])   && check_edges_valid(eo_conf[0])   && check_edges_valid(eo_conf[1]))   ending+=edge_volh;
  
  for(int par=0;par<2;par++)
    {
      for(int ivol_eo=0;ivol_eo<ending;ivol_eo++)
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
      if(ending<loc_volh+bord_volh) set_borders_invalid(eo_conf[par]);
      else if(ending<loc_volh+bord_volh+edge_volh) set_edges_invalid(eo_conf[par]);
    }
}
