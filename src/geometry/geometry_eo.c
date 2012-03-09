#pragma once

//set the eo geometry
void set_eo_geometry()
{
  if(nissa_eo_geom_inited!=0) crash("E/O Geometry already initialized!");
  
  //check that all local sizes are multiples of 2
  int ok=1;
  for(int mu=0;mu<4;mu++) ok&=(loc_size[mu]%2==0);
  if(!ok) crash("local lattice size odd!");
  
  //set half the vol, bord and edge size
  glb_volh=glb_vol/2;
  loc_volh=loc_vol/2;
  loc_bordh=loc_bord/2;
  loc_edgeh=loc_edge/2;
  
  //set the various time-slice types
  loclx_parity=nissa_malloc("loclx_parity",loc_vol+loc_bord+loc_edge,int);
  
  loceo_of_loclx=nissa_malloc("loceo_of_loclx",loc_vol+loc_bord+loc_edge,int);
  loclx_of_loceo[0]=nissa_malloc("loclx_of_loceo",loc_vol+loc_bord+loc_edge,int);
  loclx_of_loceo[1]=loclx_of_loceo[0]+loc_volh+loc_bordh+loc_edgeh;
  loceo_neighup[0]=nissa_malloc("loceo_neighup",loc_vol+loc_bord+loc_edge,coords);
  loceo_neighdw[0]=nissa_malloc("loceo_neighdw",loc_vol+loc_bord+loc_edge,coords);
  
  loceo_neighup[1]=loceo_neighup[0]+loc_volh+loc_bordh;
  loceo_neighdw[1]=loceo_neighdw[0]+loc_volh+loc_bordh;
  
  //Label the bulk sites
  int iloc_eo[2]={0,0};
  for(int loclx=0;loclx<loc_vol+loc_bord+loc_edge;loclx++)
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
  for(int loclx=0;loclx<loc_vol+loc_bord;loclx++)
    for(int mu=0;mu<4;mu++)
      {
	//take parity and e/o corresponding site
	int par=loclx_parity[loclx];
	int loceo=loceo_of_loclx[loclx];
	
	//up movements
	int loclx_up=loclx_neighup[loclx][mu];
	
	if(loclx_up>=0 && loclx_up<loc_vol+loc_bord+loc_edge)
	  loceo_neighup[par][loceo][mu]=loceo_of_loclx[loclx_up];
	
	//dw movements
	int loclx_dw=loclx_neighdw[loclx][mu];
	if(loclx_dw>=0 && loclx_dw<loc_vol+loc_bord+loc_edge)
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
  for(int ieo=0;ieo<loc_volh;ieo++)
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

//unset the eo geometry
void unset_eo_geometry()
{
  if(nissa_eo_geom_inited==0)
    crash("asking to unset never initialized E/O Geometry!");

  master_printf("Unsetting E/O Geometry\n");
  
  nissa_free(loclx_of_loceo[0]);
  nissa_free(loceo_neighup[0]);
  nissa_free(loceo_neighdw[0]);
  nissa_free(loclx_parity);
  nissa_free(loceo_of_loclx);
  
  nissa_eo_geom_inited=0;
}

//add or remove staggered phases to/from a conf
void addrem_stagphases_to_eo_conf(quad_su3 **eo_conf)
{
  if(nissa_eo_geom_inited==0) set_eo_geometry();
  
  for(int ieo=0;ieo<2;ieo++)
    for(int ivol_eo=0;ivol_eo<loc_volh;ivol_eo++)
      {
	int ivol_lx=loclx_of_loceo[ieo][ivol_eo];
	
	int d=0;
	
	//phase in direction 1 is always 0 so nothing has to be done in that dir
	//if(d%2==1) su3_prod_double(eo_conf[ieo][ivol_eo][1],eo_conf[ieo][ivol_eo][1],-1);
	
	//direction 2
	d+=glb_coord_of_loclx[ivol_lx][1];
	if(d%2==1) su3_prod_double(eo_conf[ieo][ivol_eo][2],eo_conf[ieo][ivol_eo][2],-1);
	
	//direction 3
	d+=glb_coord_of_loclx[ivol_lx][2];
	if(d%2==1) su3_prod_double(eo_conf[ieo][ivol_eo][3],eo_conf[ieo][ivol_eo][3],-1);
	
	//direction 0
	d+=glb_coord_of_loclx[ivol_lx][3];
	//debug: putting the anti-periodic condition on the temporal border
	//in future remove it!!!
	if(glb_coord_of_loclx[ivol_lx][0]==glb_size[0]-1) d+=1;
	if(d%2==1) su3_prod_double(eo_conf[ieo][ivol_eo][0],eo_conf[ieo][ivol_eo][0],-1);
      }
}

//separate the even or odd part of a vector
void take_e_or_o_part_of_lx_vector(char *out_e_or_o,char *in_lx,int nsites,int bps,int par)
{
  if(nsites<0||nsites>loc_vol+loc_bord+loc_edge) crash("too large or small vector: %d",nsites);
  
  //extract
  for(int loclx=0;loclx<nsites;loclx++)
    if(par==loclx_parity[loclx])
      memcpy(out_e_or_o+bps*loceo_of_loclx[loclx],in_lx+bps*loclx,bps);
}
//wrappers
void take_e_part_of_lx_vector(char *out_e,char *in_lx,int nsites,int bps)
{take_e_or_o_part_of_lx_vector(out_e,in_lx,nsites,bps,EVN);}
void take_o_part_of_lx_vector(char *out_o,char *in_lx,int nsites,int bps)
{take_e_or_o_part_of_lx_vector(out_o,in_lx,nsites,bps,ODD);}
void take_e_part_of_lx_color(color *out_e,color *in_lx,int nsites)
{take_e_part_of_lx_vector((char*)out_e,(char*)in_lx,nsites,sizeof(color));}
void take_o_part_of_lx_color(color *out_o,color *in_lx,int nsites)
{take_o_part_of_lx_vector((char*)out_o,(char*)in_lx,nsites,sizeof(color));}

//separate the even and odd part of a vector
void split_lx_vector_into_eo_parts(char *out_ev,char *out_od,char *in_lx,int nsites,int bps)
{
  if(nsites<0||nsites>loc_vol+loc_bord+loc_edge) crash("too large or small vector: %d",nsites);
  
  //collect even and odd sites pointer into a more convenient structure
  char *out[2]={out_ev,out_od};
  
  //split
  for(int loclx=0;loclx<nsites;loclx++)
    memcpy(out[loclx_parity[loclx]]+bps*loceo_of_loclx[loclx],in_lx+bps*loclx,bps);
}

//paste the even and odd parts of a vector into a full lx vector
void paste_eo_parts_into_lx_vector(char *out_lx,char *in_ev,char *in_od,int nsites,int bps)
{
  if(nsites<0||nsites>loc_vol+loc_bord+loc_edge) crash("too large or small vector: %d",nsites);
  
  //collect even and odd sites pointer into a more convenient structure
  char *in[2]={in_ev,in_od};
  
  //paste
  for(int loclx=0;loclx<nsites;loclx++)
    {
      int eo=loceo_of_loclx[loclx];
      int par=loclx_parity[loclx];
      memcpy(out_lx+bps*loclx,in[par]+bps*eo,bps);
    }
}

//wrappers
void split_lx_conf_into_eo_parts(quad_su3 **eo_out,quad_su3 *lx_in,int nsites)
{split_lx_vector_into_eo_parts((char*)(eo_out[EVN]),(char*)(eo_out[ODD]),(char*)lx_in,nsites,sizeof(quad_su3));}
void split_lx_color_into_eo_parts(color **eo_out,color *lx_in,int nsites)
{split_lx_vector_into_eo_parts((char*)(eo_out[EVN]),(char*)(eo_out[ODD]),(char*)lx_in,nsites,sizeof(color));}
void split_lx_spincolor_into_eo_parts(spincolor **eo_out,spincolor *lx_in,int nsites)
{split_lx_vector_into_eo_parts((char*)(eo_out[EVN]),(char*)(eo_out[ODD]),(char*)lx_in,nsites,sizeof(spincolor));}

void paste_eo_parts_into_lx_spincolor(spincolor *out_lx,spincolor *in_ev,spincolor *in_od,int nsites)
{paste_eo_parts_into_lx_vector((char*)out_lx,(char*)in_ev,(char*)in_od,nsites,sizeof(spincolor));}
