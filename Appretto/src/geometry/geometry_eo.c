#pragma once

//indexes run as t,x,y,z (faster:z)
void set_eo_geometry()
{
  if(appretto_eo_geom_inited==0)
    {
      appretto_eo_geom_inited=1;
      
      //set the various time-slice types
      loclx_parity=(int*)malloc(sizeof(int)*(loc_vol+loc_bord));
      /*
      loceo_of_loclx=(int*)malloc(sizeof(int)*(loc_vol+loc_bord));
      loclx_of_loce=(int*)malloc(sizeof(int)*(loc_vol+loc_bord)/2);
      loclx_of_loco=(int*)malloc(sizeof(int)*(loc_vol+loc_bord)/2);
      loce_neighup=(int**)malloc(sizeof(int*)*(loc_vol+loc_bord)/2);
      loce_neighdw=(int**)malloc(sizeof(int*)*(loc_vol+loc_bord)/2);
      loco_neighup=(int**)malloc(sizeof(int*)*(loc_vol+loc_bord)/2);
      loco_neighdw=(int**)malloc(sizeof(int*)*(loc_vol+loc_bord)/2);
      
      for(int loceo=0;loceo<(loc_vol+loc_bord)/2;loceo++)
	{
	  loce_neighup[loceo]=(int*)malloc(sizeof(int)*4);
	  loce_neighdw[loceo]=(int*)malloc(sizeof(int)*4);
	  loco_neighup[loceo]=(int*)malloc(sizeof(int)*4);
	  loco_neighdw[loceo]=(int*)malloc(sizeof(int)*4);
	}
      */

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
      
      if(rank==0 && debug_lvl) printf("E/O Geometry intialized\n");
      */
    }
  
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
