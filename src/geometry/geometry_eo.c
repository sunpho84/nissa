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
void initialize_eo_bord_senders_of_kind(MPI_Datatype *MPI_EO_BORD_SEND,MPI_Datatype *base)
{
  //Various type useful for edges and sub-borders
  MPI_Datatype MPI_EO_3_SLICE;
  MPI_Datatype MPI_EO_23_SLICE;
  MPI_Type_contiguous(loc_size[3]/2,*base,&MPI_EO_3_SLICE);
  MPI_Type_contiguous(loc_size[2]*loc_size[3]/2,*base,&MPI_EO_23_SLICE);

  ///////////define the sender for the 4 kinds of borders////////////
  MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3]/2,*base,&(MPI_EO_BORD_SEND[0]));
  MPI_Type_vector(loc_size[0],1,loc_size[1],MPI_EO_23_SLICE,&(MPI_EO_BORD_SEND[1]));
  MPI_Type_vector(loc_size[0]*loc_size[1],1,loc_size[2],MPI_EO_3_SLICE,&(MPI_EO_BORD_SEND[2]));
  MPI_Type_commit(&(MPI_EO_BORD_SEND[2]));
  MPI_Type_vector(loc_size[0]*loc_size[1]*loc_size[2],1,loc_size[3]/2,*base,&(MPI_EO_BORD_SEND[3]));
  //Commit
  for(int ibord=0;ibord<4;ibord++) MPI_Type_commit(&(MPI_EO_BORD_SEND[ibord]));
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
void set_eo_bord_senders_and_receivers(MPI_Datatype *MPI_EO_BORD_SEND,MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base)
{
  initialize_eo_bord_senders_of_kind(MPI_EO_BORD_SEND,base);
  initialize_eo_bord_receivers_of_kind(MPI_EO_BORD_RECE,base);
}

/*
//definitions of lexical ordered receivers for edges
void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base)
{
  //define the 6 edges receivers, which are contiguous in memory
  MPI_Type_contiguous(loc_size[2]*loc_size[3],*base,&(MPI_EDGE_RECE[0]));
  MPI_Type_contiguous(loc_size[1]*loc_size[3],*base,&(MPI_EDGE_RECE[1]));
  MPI_Type_contiguous(loc_size[1]*loc_size[2],*base,&(MPI_EDGE_RECE[2]));
  MPI_Type_contiguous(loc_size[0]*loc_size[3],*base,&(MPI_EDGE_RECE[3]));
  MPI_Type_contiguous(loc_size[0]*loc_size[2],*base,&(MPI_EDGE_RECE[4]));
  MPI_Type_contiguous(loc_size[0]*loc_size[1],*base,&(MPI_EDGE_RECE[5]));
  for(int iedge=0;iedge<6;iedge++) MPI_Type_commit(&(MPI_EDGE_RECE[iedge]));
}
*/

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
