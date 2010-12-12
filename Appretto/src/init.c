#pragma once

#include <mpi.h>
#include <stdio.h>

#include "endianess.c"
#include "random.c"
#include "dirac.c"

//Return the index of site of coord x in the border idir,jdir
int edgelx_of_coord(int *x,int idir,int jdir)
{
  int ilx=0;
  
  for(int i=0;i<4;i++) if(i!=idir && i!=jdir) ilx=ilx*loc_size[i]+x[i];
  
  return ilx;
}

//Return the index of site of coord x in the border idir
int bordlx_of_coord(int *x,int idir)
{
  int ilx=0;
  
  for(int i=0;i<4;i++) if(i!=idir) ilx=ilx*loc_size[i]+x[i];
  
  return ilx;
}

int bordlx_of_coord_list(int x0,int x1,int x2,int x3,int idir)
{
  int x[4]={x0,x1,x2,x3};
  return bordlx_of_coord(x,idir);
}

//Return the index of site of coord x in a box of sides s
int lx_of_coord(int *x,int *s)
{
  int ilx=0;
  
  for(int i=0;i<4;i++) ilx=ilx*s[i]+x[i];
  
  return ilx;
}

//wrappers
int loclx_of_coord(int *x)
{
  return lx_of_coord(x,loc_size);
}
  
//wrappers
int loclx_of_coord_list(int x0,int x1,int x2,int x3)
{
  int x[4]={x0,x1,x2,x3};
  return lx_of_coord(x,loc_size);
}
  
//wrappers
int glblx_of_coord(int *x)
{
  return lx_of_coord(x,glb_size);
}
  
void init_appretto()
{
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&rank_tot);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //define the gauge link
  MPI_Type_contiguous(18,MPI_DOUBLE,&MPI_SU3);
  MPI_Type_commit(&MPI_SU3);

  //four links starting from a single point
  MPI_Type_contiguous(4,MPI_SU3,&MPI_QUAD_SU3);
  MPI_Type_commit(&MPI_QUAD_SU3);

  check_endianess();
  init_base_gamma();
}

//indexes run as t,z,y,x (faster:x)
void set_geometry()
{
  //set the various time-slice types

  //find the rank of the neighbour in the various dir
  for(int idir=0;idir<4;idir++)
    MPI_Cart_shift(cart_comm,idir,1,&(rank_neighdw[idir]),&(rank_neighup[idir]));
  
  loc_coord_of_loclx=(int**)malloc(sizeof(int*)*loc_vol);
  glb_coord_of_loclx=(int**)malloc(sizeof(int*)*loc_vol);
  loclx_neighup=(int**)malloc(sizeof(int*)*(loc_vol+loc_bord));
  loclx_neighdw=(int**)malloc(sizeof(int*)*(loc_vol+loc_bord));
  for(int loc_ind=0;loc_ind<loc_vol;loc_ind++)
    {
      loc_coord_of_loclx[loc_ind]=(int*)malloc(sizeof(int)*4);
      glb_coord_of_loclx[loc_ind]=(int*)malloc(sizeof(int)*4);
    }

  for(int loc_ind=0;loc_ind<(loc_vol+loc_bord);loc_ind++)
    {
      loclx_neighup[loc_ind]=(int*)malloc(sizeof(int)*4);
      loclx_neighdw[loc_ind]=(int*)malloc(sizeof(int)*4);
    }

  glblx_of_loclx=(int*)malloc(sizeof(int)*loc_vol);

  //Label the sites
  int x[4],gx[4];
  for(x[0]=0;x[0]<loc_size[0];x[0]++) 
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    for(int i=0;i<4;i++) gx[i]=x[i]+proc_coord[i]*loc_size[i];

	    int loc_ind,glb_ind;
	    loc_ind=loclx_of_coord(x);
	    glb_ind=glblx_of_coord(gx);

	    for(int i=0;i<4;i++)
	      {
		loc_coord_of_loclx[loc_ind][i]=x[i];
		glb_coord_of_loclx[loc_ind][i]=gx[i];
	      }

	    glblx_of_loclx[loc_ind]=glb_ind;
          }
  

  //////////////////neighbours search//////////////////////

  //now fill the neighbours of sites of the bulk, and the defined
  //neighbours of the sites of the external borders
  for(int iloc=0;iloc<loc_vol;iloc++)
    {
      for(int idir=0;idir<4;idir++) x[idir]=loc_coord_of_loclx[iloc][idir];
      
      //Direction on the whole iper-cube
      for(int idir=0;idir<4;idir++)
	{
	  //Down direction
	  if(x[idir]!=0)
	    {
	      x[idir]--;
	      loclx_neighdw[iloc][idir]=loclx_of_coord(x);
	      x[idir]++;
	    }
	  else //border
	    {
	      int ibord=loc_vol+bord_offset[idir]+bordlx_of_coord(x,idir);
	      loclx_neighdw[iloc][idir]=ibord;
	      loclx_neighup[ibord][idir]=iloc;
	      //the movement along the down direction from the border are not defined 
	      loclx_neighdw[ibord][idir]=-1;

	      //This is the bad moment: the movent inside the cube
	      for(int jdir=0;jdir<4;jdir++)
		if(idir!=jdir) 
		  {
		    //Down direction
		    if(x[jdir]!=0)
		      {
			x[jdir]--;
			int ibord2=loc_vol+bord_offset[idir]+bordlx_of_coord(x,idir);
			x[jdir]++;
			loclx_neighdw[ibord][jdir]=ibord2;
		      }
		    else //edge
		      {
			int iedge=loc_vol+loc_bord+edge_offset[edge_numb[idir][jdir]]+edgelx_of_coord(x,idir,jdir);
			loclx_neighdw[ibord][jdir]=iedge;
		      }
		    //Upper direction
		    if(x[jdir]!=loc_size[jdir]-1)
		      {
			x[jdir]++;
			int ibord2=loc_vol+bord_offset[idir]+bordlx_of_coord(x,idir);
			x[jdir]--;
			loclx_neighup[ibord][jdir]=ibord2;
		      }
		    else
		      {
			int iedge=loc_vol+loc_bord+edge_offset[edge_numb[idir][jdir]]+edgelx_of_coord(x,idir,jdir);
			if(idir<jdir) iedge+=loc_edge/4; //edge i-j+
			else iedge+=loc_edge/2; //edge j+i-
			loclx_neighup[ibord][jdir]=iedge;
		      }
		  }
	    }
	  
	  //Upper direction
	  if(x[idir]!=loc_size[idir]-1)
	    {
	      x[idir]++;
	      loclx_neighup[iloc][idir]=loclx_of_coord(x);
	      x[idir]--;
	    }
	  else //border
	    {
	      int ibord=loc_vol+bord_offset[idir]+loc_bord/2+bordlx_of_coord(x,idir);
	      loclx_neighup[iloc][idir]=ibord;
	      loclx_neighdw[ibord][idir]=iloc;
	      //the movement along the up direction from the border are not defined 
	      loclx_neighup[ibord][idir]=-1;

	      //Another very bad moment: the movent inside the cube
	      for(int jdir=0;jdir<4;jdir++)
		if(idir!=jdir) 
		  {
		    //Down direction
		    if(x[jdir]!=0)
		      {
			x[jdir]--;
			int ibord2=loc_vol+bord_offset[idir]+loc_bord/2+bordlx_of_coord(x,idir);
			x[jdir]++;
			loclx_neighdw[ibord][jdir]=ibord2;
		      }
		    else //edge
		      {
			int iedge=loc_vol+loc_bord+edge_offset[edge_numb[idir][jdir]]+edgelx_of_coord(x,idir,jdir);
			if(idir<jdir) iedge+=loc_edge/2; //edge i-j+
			else iedge+=loc_edge/4; //edge j+i-

			loclx_neighdw[ibord][jdir]=iedge;
		      }
		    //Upper direction
		    if(x[jdir]!=loc_size[jdir]-1)
		      {
			x[jdir]++;
			int ibord2=loc_vol+bord_offset[idir]+loc_bord/2+bordlx_of_coord(x,idir);
			x[jdir]--;
			loclx_neighup[ibord][jdir]=ibord2;
		      }
		    else //edge
		      {
			int iedge=loc_vol+loc_bord+edge_offset[edge_numb[idir][jdir]]+3*loc_edge/4+edgelx_of_coord(x,idir,jdir);
			loclx_neighup[ibord][jdir]=iedge;
		      }
		  }

	    }
	}
    }
    
  if(rank==0 && debug) printf("Geometry intialized\n");
}

void init_grid()
{
  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      tic=MPI_Wtime();
    }


  int periods[4]={1,1,1,1};
  char proc_name[1024];
  int proc_name_length;

  glb_size[2]=glb_size[3]=glb_size[1];
  MPI_Bcast(glb_size,4,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  //calculate global volume, initialize local one
  glb_vol=1;
  for(int idir=0;idir<4;idir++)
    {
      loc_size[idir]=glb_size[idir];
      glb_vol*=glb_size[idir];
    }

  if(rank==0)
    {
      printf("\nNumber of running processes: %d\n",rank_tot);
      printf("Global lattice:\t%dx%dx%dx%d = %d\n",glb_size[0],glb_size[1],glb_size[2],glb_size[3],glb_vol);
    }

  MPI_Get_processor_name(proc_name,&proc_name_length);
  MPI_Dims_create(rank_tot,4,nproc_dir);
  if(rank==0 && debug==1)
    {
      printf("\nCreating grid:\t%dx%dx%dx%d\n",nproc_dir[0],nproc_dir[1],nproc_dir[2],nproc_dir[3]);
      fflush(stdout);
    }

  int ok=1;
  for(int idir=0;idir<4;idir++)
    {
      ok=ok && (nproc_dir[idir]>0);
      ok=ok && (glb_size[idir]%nproc_dir[idir]==0);
    }

  if(!ok && rank==0)
    {
      fprintf(stderr,"The lattice is incommensurable with the total processor amount\n");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  //Calculate local volume
  for(int idir=0;idir<4;idir++) loc_size[idir]=glb_size[idir]/nproc_dir[idir];
  loc_vol=glb_vol/rank_tot;
  
  //Calculate the border size
  loc_bord=0;
  bord_offset[0]=0;
  for(int idir=0;idir<4;idir++)
    {
      //bord size along the idir dir
      bord_dir_vol[idir]=loc_vol/loc_size[idir];

      //total bord
      loc_bord+=bord_dir_vol[idir];

      //summ of the border extent up to dir idir
      if(idir>0) bord_offset[idir]=bord_offset[idir-1]+bord_dir_vol[idir-1];
    }
  loc_bord*=2;
  
  //Calculate the egdes size
  loc_edge=0;
  edge_offset[0]=0;
  int iedge=0;
  for(int idir=0;idir<4;idir++)
    for(int jdir=idir+1;jdir<4;jdir++)
      {
	//edge among the i and j dir
	edge_dir_vol[iedge]=bord_dir_vol[idir]/loc_size[jdir];
	
	//total edge
	loc_edge+=edge_dir_vol[iedge];
	
	//summ of the border extent up to dir i
	if(iedge>0)
	  edge_offset[iedge]=edge_offset[iedge-1]+edge_dir_vol[iedge-1];
	iedge++;
    }
  loc_edge*=4;
  
  if(rank==0)
    {
      printf("Local volume\t%dx%dx%dx%d = %d\n",loc_size[0],loc_size[1],loc_size[2],loc_size[3],loc_vol);
      if(debug>1) printf("Border size: %d\n",loc_bord);
      if(debug>1) printf("Edge size: %d\n",loc_edge);
      if(debug>2) 
	for(int idir=0;idir<4;idir++)
	  printf("Border offset for dir %d: %d\n",idir,bord_offset[idir]);
      if(debug>2)
	for(int iedge=0;iedge<6;iedge++)
	  printf("Border offset for edge %d: %d\n",iedge,edge_offset[iedge]);
    }
  MPI_Cart_create(MPI_COMM_WORLD,4,nproc_dir,periods,1,&cart_comm);
  MPI_Comm_rank(cart_comm,&cart_rank);
  MPI_Cart_coords(cart_comm,cart_rank,4,proc_coord);

  if(debug>2)
    {
      printf("Process %d of %d on %s: %d (%d %d %d %d)\n",rank,rank_tot,
	     proc_name,cart_rank,proc_coord[0],proc_coord[1],proc_coord[2],proc_coord[3]);
      fflush(stdout);
  
      MPI_Barrier(MPI_COMM_WORLD);
    }

  //Various type useful for edges and sub-borders
  MPI_Datatype MPI_GAUGE_3_SLICE;
  MPI_Datatype MPI_GAUGE_23_SLICE;
  MPI_Type_contiguous(loc_size[3],MPI_QUAD_SU3,&MPI_GAUGE_3_SLICE);
  MPI_Type_contiguous(loc_size[2]*loc_size[3],MPI_QUAD_SU3,&MPI_GAUGE_23_SLICE);
  MPI_Type_commit(&MPI_GAUGE_3_SLICE);
  MPI_Type_commit(&MPI_GAUGE_23_SLICE);

  //define the 0-dir slice
  MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3],MPI_QUAD_SU3,&(MPI_GAUGE_SLICE_RECE[0]));
  MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3],MPI_QUAD_SU3,&(MPI_GAUGE_SLICE_SEND[0]));
  MPI_Type_commit(&(MPI_GAUGE_SLICE_RECE[0]));
  MPI_Type_commit(&(MPI_GAUGE_SLICE_SEND[0]));
  
  //define the 1-dir slice
  MPI_Type_contiguous(loc_size[0]*loc_size[2]*loc_size[3],MPI_QUAD_SU3,&(MPI_GAUGE_SLICE_RECE[1]));
  MPI_Type_vector(loc_size[0],1,loc_size[1],MPI_GAUGE_23_SLICE,&(MPI_GAUGE_SLICE_SEND[1]));
  MPI_Type_commit(&(MPI_GAUGE_SLICE_RECE[1]));
  MPI_Type_commit(&(MPI_GAUGE_SLICE_SEND[1]));

  //define the 2-dir slice
  MPI_Type_contiguous(loc_size[0]*loc_size[1]*loc_size[3],MPI_QUAD_SU3,&(MPI_GAUGE_SLICE_RECE[2]));
  MPI_Type_vector(loc_size[0]*loc_size[1],1,loc_size[2],MPI_GAUGE_3_SLICE,&(MPI_GAUGE_SLICE_SEND[2]));
  MPI_Type_commit(&(MPI_GAUGE_SLICE_RECE[2]));
  MPI_Type_commit(&(MPI_GAUGE_SLICE_SEND[2]));
  
  //define the 3-dir slice
  MPI_Type_contiguous(loc_size[0]*loc_size[1]*loc_size[2],MPI_QUAD_SU3,&(MPI_GAUGE_SLICE_RECE[3]));
  MPI_Type_vector(loc_size[0]*loc_size[1]*loc_size[2],1,loc_size[3],MPI_QUAD_SU3,&(MPI_GAUGE_SLICE_SEND[3]));
  MPI_Type_commit(&(MPI_GAUGE_SLICE_RECE[3]));
  MPI_Type_commit(&(MPI_GAUGE_SLICE_SEND[3]));
  
  ///////////define the sender and receiver for the 6 kinds of edges////////////
  
  //this is the 01 sender, that is simply a vector of L[2] vector of L[3] 
  MPI_Type_contiguous(loc_size[2]*loc_size[3],MPI_QUAD_SU3,&(MPI_GAUGE_EDGE_SEND[0]));
  MPI_Type_contiguous(loc_size[2]*loc_size[3],MPI_QUAD_SU3,&(MPI_GAUGE_EDGE_RECE[0]));
  
  //this is the 02 sender, a vector of L[1] segment of
  //length L[3] (already defined) separated by L[2] of them
  MPI_Type_vector(loc_size[1],1,loc_size[2],MPI_GAUGE_3_SLICE,&(MPI_GAUGE_EDGE_SEND[1]));
  MPI_Type_contiguous(loc_size[1]*loc_size[3],MPI_QUAD_SU3,&(MPI_GAUGE_EDGE_RECE[1]));
  
  //this is the 03 sender, a segment of length L[1] of a
  //vector of L[2] single elements, separated by L[3] of them
  MPI_Datatype MPI_EDGE_23_SUB_SLICE;
  MPI_Type_vector(loc_size[2],1,loc_size[3],MPI_QUAD_SU3,&(MPI_EDGE_23_SUB_SLICE)); 
  MPI_Type_commit(&MPI_EDGE_23_SUB_SLICE);
  MPI_Type_contiguous(loc_size[1],MPI_EDGE_23_SUB_SLICE,&(MPI_GAUGE_EDGE_SEND[2]));
  MPI_Type_contiguous(loc_size[1]*loc_size[2],MPI_QUAD_SU3,&(MPI_GAUGE_EDGE_RECE[2]));
  
  //this is the 12 sender, should be equal to the 02 sender, with 1->0
  MPI_Type_vector(loc_size[0],1,loc_size[2],MPI_GAUGE_3_SLICE,&(MPI_GAUGE_EDGE_SEND[3]));
  MPI_Type_contiguous(loc_size[0]*loc_size[3],MPI_QUAD_SU3,&(MPI_GAUGE_EDGE_RECE[3]));
  
  //this is the 13 sender, should be equal to the 03 sender, with 1->0
  MPI_Type_contiguous(loc_size[0],MPI_EDGE_23_SUB_SLICE,&(MPI_GAUGE_EDGE_SEND[4]));
  MPI_Type_contiguous(loc_size[0]*loc_size[2],MPI_QUAD_SU3,&(MPI_GAUGE_EDGE_RECE[4]));
  
  //this is the 23 sender, should be equal to the 03 sender with 1<->2
  MPI_Datatype MPI_EDGE_13_SUB_SLICE;
  MPI_Type_vector(loc_size[1],1,loc_size[3],MPI_QUAD_SU3,&(MPI_EDGE_13_SUB_SLICE)); 
  MPI_Type_commit(&MPI_EDGE_13_SUB_SLICE);
  MPI_Type_contiguous(loc_size[0],MPI_EDGE_13_SUB_SLICE,&(MPI_GAUGE_EDGE_SEND[5]));
  MPI_Type_contiguous(loc_size[0]*loc_size[1],MPI_QUAD_SU3,&(MPI_GAUGE_EDGE_RECE[5]));
  
  for(int iedge=0;iedge<6;iedge++)
    {
      MPI_Type_commit(&(MPI_GAUGE_EDGE_SEND[iedge]));
      MPI_Type_commit(&(MPI_GAUGE_EDGE_RECE[iedge]));
    }

  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      tac=MPI_Wtime();

      if(rank==0) printf("Time elapsed for MPI inizialization: %f s\n",tac-tic);
    }

  set_geometry();
}
