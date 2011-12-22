#pragma once

//set the eo geometry
void set_eo_geometry()
{
  if(nissa_eo_geom_inited==0)
    {
      //check that all size are multiples of 2
      int ok=1;
      for(int mu=0;mu<4;mu++) ok&=(loc_size[mu]%2==0);
      if(!ok) crash("local lattice size odd!");
      
      //set the various time-slice types
      loclx_parity=nissa_malloc("loclx_parity",loc_vol+loc_bord,int);
      
      loceo_of_loclx=nissa_malloc("loceo_of_loclx",loc_vol+loc_bord,int);
      loclx_of_loceo[0]=nissa_malloc("loclx_of_loceo",loc_vol+loc_bord,int);
      loclx_of_loceo[1]=loclx_of_loceo[0]+(loc_vol+loc_bord)/2;
      loceo_neighup[0]=nissa_malloc("loceo_neighup",loc_vol+loc_bord,int*);
      loceo_neighdw[0]=nissa_malloc("loceo_neighdw",loc_vol+loc_bord,int*);
      loceo_neighup[0][0]=nissa_malloc("loceo_neighup[0]",4*(loc_vol+loc_bord),int);
      loceo_neighdw[0][0]=nissa_malloc("loceo_neighdw[0]",4*(loc_vol+loc_bord),int);

      loceo_neighup[1]=loceo_neighup[0]+(loc_vol+loc_bord)/2;
      loceo_neighdw[1]=loceo_neighdw[0]+(loc_vol+loc_bord)/2;      
      for(int eo=0;eo<2;eo++)
	for(int loceo=1;loceo<(loc_vol+loc_bord)/2;loceo++)
	  {
	    loceo_neighup[eo][loceo]=loceo_neighup[eo][loceo-1]+4;
	    loceo_neighdw[eo][loceo]=loceo_neighdw[eo][loceo-1]+4;
	  }
      
      //int iloce=0,iloco=0;
      //Label the bulk sites
      for(int loclx=0;loclx<loc_vol+loc_bord;loclx++)
	{
	  //calculate global coord and parity
	  int par=0;
	  for(int idir=0;idir<4;idir++) par+=glb_coord_of_loclx[loclx][idir];
	  par%=2;
	  
	  //fix parity of local index
	  loclx_parity[loclx]=par;
	  /*
	  //assign the lx site of the e/o site, and
	  //increment the number of even or odd sites
	  if(par==0)
	    {
	      loceo_of_loclx[loclx]=iloce;
	      loclx_of_loce[iloce]=loclx;
	      iloce++;
	    }
	  else
	    {
	      loceo_of_loclx[loclx]=iloco;
	      loclx_of_loco[iloco]=loclx;
	      iloco++;
	    }
	  */
	}
      
      //////////////////neighbours search//////////////////////
      /*
      memset(bord_offset_eo,0,sizeof(int)*2*8);
      
      //now fill the neighbours of sites
      for(int loclx=0;loclx<loc_vol+loc_bord;loclx++)
	{
	  int loceo=loceo_of_loclx[loclx];
	  for(int idir=0;idir<4;idir++)
	    if(loclx_parity[loclx]==0)
	      {
		loce_neighdw[loceo][idir]=loceo_of_loclx[loclx_neighdw[loclx][idir]];
		loce_neighup[loceo][idir]=loceo_of_loclx[loclx_neighup[loclx][idir]];
	      }
	    else
	      {
		loco_neighdw[loceo][idir]=loceo_of_loclx[loclx_neighdw[loclx][idir]];
		loco_neighup[loceo][idir]=loceo_of_loclx[loclx_neighup[loclx][idir]];
	      }
	  
	  //count the number of points in the even and odd border
	  if(loclx>=loc_vol)
	    {
	      int ibord=loclx-loc_vol;
	      int ibord_dir=0;
	      if(loclx>=(loc_vol+loc_bord/2))
		{
		  ibord-=loc_vol/2;
		  ibord_dir+=4;
		}
	      ibord_dir+=dir_of_bord[ibord];
	      
	      int eo=loclx_parity[loclx];
	      bord_offset_eo[eo][ibord_dir]++;
	    }	      
	}
      */
      
      master_printf("E/O Geometry intialized\n");
      
      nissa_eo_geom_inited=1;
    }
  else crash("E/O Geometry already initialized!");
  
  //if(rank==0) for(int idir=0;idir<8;idir++) printf("%d %d %d\n",idir,bord_offset_eo[0][idir],bord_offset_eo[1][idir]);
}

/*

//definitions of lexical ordered sender for borders
void initialize_e_bord_senders_of_kind(MPI_Datatype *MPI_E_BORD_SEND,MPI_Datatype *base)
{
  //Various type useful for edges and sub-borders
  MPI_Datatype MPI_EO_3_SLICE;
  MPI_Datatype MPI_EO_23_SLICE;
  MPI_Type_contiguous(loc_size[3]/2,*base,&MPI_EO_3_SLICE);
  MPI_Type_contiguous(loc_size[2]*loc_size[3]/2,*base,&MPI_EO_23_SLICE);

  ///////////define the sender for the 4 kinds of borders////////////
  MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3]/2,*base,&(MPI_EO_BORD_SEND[0]));
  MPI_Type_vector(loc_size[0],1,loc_size[1],MPI_EO_23_SLICE,&(MPI_EO_BORD_SEND[1]));
  MPI_Type_vector(loc_size[0]*loc_size[1],1,loc_size[2],MPI_EO_3_SLICE,&(MPI_BORD_EO_SEND[2]));
  MPI_Type_commit(&(MPI_EO_BORD_SEND[2]));
  MPI_Type_vector(loc_size[0]*loc_size[1]*loc_size[2]/2,1,loc_size[3],*base,&(MPI_BORD_EO_SEND[3]));
  //Commit
  for(int ibord=0;ibord<4;ibord++) MPI_Type_commit(&(MPI_EO_BORD_SEND[ibord]));
}

//definitions of lexical ordered receivers for borders
void initialize_lx_bord_receivers_of_kind(MPI_Datatype *MPI_BORD_RECE,MPI_Datatype *base)
{
  //define the 4 dir borders receivers, which are contiguous in memory
  MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3],*base,&(MPI_BORD_RECE[0]));
  MPI_Type_contiguous(loc_size[0]*loc_size[2]*loc_size[3],*base,&(MPI_BORD_RECE[1]));
  MPI_Type_contiguous(loc_size[0]*loc_size[1]*loc_size[3],*base,&(MPI_BORD_RECE[2]));
  MPI_Type_contiguous(loc_size[0]*loc_size[1]*loc_size[2],*base,&(MPI_BORD_RECE[3]));
  for(int ibord=0;ibord<4;ibord++) MPI_Type_commit(&(MPI_BORD_RECE[ibord]));
}

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
  if(nissa_eo_geom_inited==0) crash("asking to unset never initialized E/O Geometry!");

  master_printf("Unsetting E/O Geometry\n");
  
  nissa_free(loceo_neighup[0][0]);
  nissa_free(loceo_neighdw[0][0]);
  nissa_free(loclx_of_loceo[0]);
  nissa_free(loceo_neighup[0]);
  nissa_free(loceo_neighdw[0]);
  nissa_free(loclx_parity);
  nissa_free(loceo_of_loclx);
  
  nissa_eo_geom_inited=0;
}
