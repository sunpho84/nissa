#pragma once

#include <mpi.h>
#include <stdio.h>

//indexes run as t,z,y,x (faster:x)
void set_eo_geometry()
{
  //set the various time-slice types
  loclx_parity=(int*)malloc(sizeof(int)*loc_vol);
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

  int iloce=0,iloco=0;
  //Label the bulk sites
  for(int loclx=0;loclx<loc_vol+loc_bord;loclx++)
    {
      //calculate global coord and parity
      int par=0;
      for(int idir=0;idir<4;idir++) par+=glb_coord_of_loclx[loclx][idir];
      par%=2;
      
      //fix parity of local index
      loclx_parity[loclx]=par;
      
      //assign the lx site of the e/o site, and
      //increment the number of even or odd sites
      if(par==0) loclx_of_loce[iloce++]=loclx;
      else       loclx_of_loce[iloco++]=loclx;
    }
  
  //////////////////neighbours search//////////////////////
  
  //now fill the neighbours of sites
  for(int loclx=0;loclx<loc_vol+loc_bord;loclx++)
    {
      int loceo=loceo_of_loclx[loclx];
      int **loceo_neighdw,**loceo_neighup;
      if(loclx_parity[loclx]==0)
	{
	  loceo_neighdw=loce_neighdw;
	  loceo_neighup=loce_neighup;
	}
      else
	{
	  loceo_neighdw=loco_neighdw;
	  loceo_neighup=loco_neighup;
	}
      
      //Direction on the whole iper-cube
      for(int idir=0;idir<4;idir++)
	{
	  loceo_neighdw[loceo][idir]=loceo_of_loclx[loclx_neighdw[loclx][idir]];
	  loceo_neighup[loceo][idir]=loceo_of_loclx[loclx_neighup[loclx][idir]];
	}
    }

  if(rank==0 && debug) printf("E/O Geometry intialized\n");
}
