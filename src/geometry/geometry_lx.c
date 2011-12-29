#pragma once

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
{return lx_of_coord(x,loc_size);}
  
//wrappers
int loclx_of_coord_list(int x0,int x1,int x2,int x3)
{
  int x[4]={x0,x1,x2,x3};
  return lx_of_coord(x,loc_size);
}
  
//wrappers
int glblx_of_coord(int *x)
{return lx_of_coord(x,glb_size);}

//Return the coordinate of the proc containing the global coord
void proc_coord_of_site_of_coord(int *proc_coord,int *glb_coord)
{for(int mu=0;mu<4;mu++) proc_coord[mu]=glb_coord[mu]/loc_size[mu];}

//Return the rank of passed coord
int proc_of_coord(int *coord)
{return lx_of_coord(coord,nproc_dir);}

//Return the rank containing the global coordinates
int rank_hosting_site_of_coord(int *x)
{
  int p[4];
  proc_coord_of_site_of_coord(p,x);
  
  return proc_of_coord(p);
}

//indexes run as t,x,y,z (faster:z)
void set_lx_geometry()
{
  if(nissa_lx_geom_inited==1) crash("cartesian geometry already intialized!");
  
  //find the rank of the neighbour in the various dir
  for(int idir=0;idir<4;idir++)
    MPI_Cart_shift(cart_comm,idir,1,&(rank_neighdw[idir]),&(rank_neighup[idir]));
  
  loc_coord_of_loclx=nissa_malloc("loc_coord_of_loclx",loc_vol,int*);
  glb_coord_of_loclx=nissa_malloc("glb_coord_of_loclx",loc_vol+loc_bord,int*);
  loclx_neighup=nissa_malloc("loclx_neighup",loc_vol+loc_bord,int*);
  loclx_neighdw=nissa_malloc("loclx_neighdw",loc_vol+loc_bord,int*);
  
  loc_coord_of_loclx[0]=nissa_malloc("loc_coord_of_loclx",4*loc_vol,int);
  glb_coord_of_loclx[0]=nissa_malloc("glb_coord_of_loclx",4*(loc_vol+loc_bord),int);
  loclx_neighup[0]=nissa_malloc("loclx_neighup",4*(loc_vol+loc_bord),int);
  loclx_neighdw[0]=nissa_malloc("loclx_neighdw",4*(loc_vol+loc_bord),int);
  
  for(int loc_ind=1;loc_ind<loc_vol;loc_ind++)
    loc_coord_of_loclx[loc_ind]=loc_coord_of_loclx[loc_ind-1]+4;
  
  for(int loc_ind=1;loc_ind<loc_vol+loc_bord;loc_ind++)
    {
      glb_coord_of_loclx[loc_ind]=glb_coord_of_loclx[loc_ind-1]+4;
      loclx_neighup[loc_ind]=loclx_neighup[loc_ind-1]+4;
      loclx_neighdw[loc_ind]=loclx_neighdw[loc_ind-1]+4;
    }
  
  //local to global
  glblx_of_loclx=nissa_malloc("glblx_of_loclx",loc_vol,int);
      
  //borders
  glblx_of_bordlx=nissa_malloc("glblx_of_bordlx",loc_bord,int);
  loclx_of_bordlx=nissa_malloc("loclx_of_bordlx",loc_bord,int);
  dir_of_bord=nissa_malloc("dir_of_bord",loc_bord,int);
      
  //edges
  glblx_of_edgelx=nissa_malloc("glblx_of_edgelx",loc_edge,int);
  
  //Label the sites
  int x[4],gx[4]={0,0,0,0};
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
      for(int idir=0;idir<4;idir++)
	{
	  x[idir]=loc_coord_of_loclx[iloc][idir];
	  gx[idir]=glb_coord_of_loclx[iloc][idir];
	}
      
      //Direction on the whole iper-cube
      for(int idir=0;idir<4;idir++)
	{
	  //Down direction
	  if(x[idir]!=0 || paral_dir[idir]==0)
	    {
	      int temp_xdir=x[idir];
	      x[idir]=(x[idir]+loc_size[idir]-1)%loc_size[idir];
	      loclx_neighdw[iloc][idir]=loclx_of_coord(x);
	      x[idir]=temp_xdir;
	    }
	  else //border
	    {
	      int raw_ibord=bordlx_of_coord(x,idir)+bord_offset[idir];
	      int ibord=loc_vol+raw_ibord;
	      loclx_neighdw[iloc][idir]=ibord;
	      loclx_neighup[ibord][idir]=iloc;
	      //the movement along the down direction from the border are not defined 
	      loclx_neighdw[ibord][idir]=-1;
	      //number the elements of the border
	      gx[idir]=(gx[idir]-1+glb_size[idir])%glb_size[idir];
	      x[idir]=(x[idir]-1+loc_size[idir])%loc_size[idir];
	      for(int idir=0;idir<4;idir++) glb_coord_of_loclx[ibord][idir]=gx[idir];
	      glblx_of_bordlx[raw_ibord]=glblx_of_coord(gx);
	      loclx_of_bordlx[raw_ibord]=loclx_of_coord(x);
	      dir_of_bord[raw_ibord]=2*idir+1;
	      gx[idir]=glb_coord_of_loclx[iloc][idir];
	      x[idir]=loc_coord_of_loclx[iloc][idir];
	      
	      //This is the bad moment: the movents inside the cube
	      for(int jdir=0;jdir<4;jdir++)
		if(idir!=jdir) 
		  {
		    //Down direction
		    if(x[jdir]!=0 || paral_dir[jdir]==0)
		      {
			int temp_xdir=x[jdir];
			x[jdir]=(x[jdir]+loc_size[jdir]-1)%loc_size[jdir];
			int ibord2=loc_vol+bord_offset[idir]+bordlx_of_coord(x,idir);
			x[jdir]=temp_xdir;
			loclx_neighdw[ibord][jdir]=ibord2;
		      }
		    else //i-j- edge
		      {
			int raw_iedge=edge_offset[edge_numb[idir][jdir]]+edgelx_of_coord(x,idir,jdir);
			//put the value of the neighbour
			int iedge=loc_vol+loc_bord+raw_iedge;
			loclx_neighdw[ibord][jdir]=iedge;
			//number the elements of the edge
			gx[idir]=(gx[idir]-1+glb_size[idir])%glb_size[idir];
			gx[jdir]=(gx[jdir]-1+glb_size[jdir])%glb_size[jdir];
			glblx_of_edgelx[raw_iedge]=glblx_of_coord(gx);
			gx[idir]=glb_coord_of_loclx[iloc][idir];
			gx[jdir]=glb_coord_of_loclx[iloc][jdir];
		      }
		    //Upper direction
		    if(x[jdir]!=loc_size[jdir]-1 || paral_dir[jdir]==0)
		      {
			int temp_xdir=x[jdir];
			x[jdir]=(x[jdir]+1)%loc_size[jdir];
			int ibord2=loc_vol+bord_offset[idir]+bordlx_of_coord(x,idir);
			x[jdir]=temp_xdir;
			loclx_neighup[ibord][jdir]=ibord2;
		      }
		    else
		      {
			int raw_iedge=edge_offset[edge_numb[idir][jdir]]+edgelx_of_coord(x,idir,jdir);
			if(idir<jdir) raw_iedge+=loc_edge/4; //edge i-j+
			else raw_iedge+=loc_edge/2; //edge j+i-
			//put the value of the neighbour
			int iedge=loc_vol+loc_bord+raw_iedge;
			loclx_neighup[ibord][jdir]=iedge;
			//number the elements of the edge
			gx[idir]=(gx[idir]-1+glb_size[idir])%glb_size[idir];
			gx[jdir]=(gx[jdir]+1+glb_size[jdir])%glb_size[jdir];
			glblx_of_edgelx[raw_iedge]=glblx_of_coord(gx);
			gx[idir]=glb_coord_of_loclx[iloc][idir];
			gx[jdir]=glb_coord_of_loclx[iloc][jdir];
		      }
		  }
	    }
	  
	  //Upper direction
	  if(x[idir]!=loc_size[idir]-1 || paral_dir[idir]==0)
	    {
	      int temp_xdir=x[idir];
	      x[idir]=(x[idir]+1)%loc_size[idir];
	      loclx_neighup[iloc][idir]=loclx_of_coord(x);
	      x[idir]=temp_xdir;
	    }
	  else //border
	    {
	      int raw_ibord=loc_bord/2+bordlx_of_coord(x,idir)+bord_offset[idir];
	      int ibord=loc_vol+raw_ibord;
	      loclx_neighup[iloc][idir]=ibord;
	      loclx_neighdw[ibord][idir]=iloc;
	      //the movement along the up direction from the border are not defined 
	      loclx_neighup[ibord][idir]=-1;
	      //number the elements of the border
	      gx[idir]=(gx[idir]+1+glb_size[idir])%glb_size[idir];
	      x[idir]=(x[idir]+1+loc_size[idir])%loc_size[idir];
	      for(int idir=0;idir<4;idir++) glb_coord_of_loclx[ibord][idir]=gx[idir];
	      glblx_of_bordlx[raw_ibord]=glblx_of_coord(gx);
	      loclx_of_bordlx[raw_ibord]=loclx_of_coord(x);
	      dir_of_bord[raw_ibord]=2*idir;
	      gx[idir]=glb_coord_of_loclx[iloc][idir];
	      x[idir]=loc_coord_of_loclx[iloc][idir];
	      
	      //Another very bad moment: the movents inside the cube
	      for(int jdir=0;jdir<4;jdir++)
		if(idir!=jdir) 
		  {
		    //Down direction
		    if(x[jdir]!=0 || paral_dir[jdir]==0)
		      {
			int temp_xdir=x[jdir];
			x[jdir]=(x[jdir]+loc_size[jdir]-1)%loc_size[jdir];
			int ibord2=loc_vol+bord_offset[idir]+loc_bord/2+bordlx_of_coord(x,idir);
			x[jdir]=temp_xdir;
			loclx_neighdw[ibord][jdir]=ibord2;
		      }
		    else //edge
		      {
			int raw_iedge=edge_offset[edge_numb[idir][jdir]]+edgelx_of_coord(x,idir,jdir);
			if(idir<jdir) raw_iedge+=loc_edge/2; //edge i-j+
			else raw_iedge+=loc_edge/4; //edge j+i-
			//put the value of the neighbour
			int iedge=loc_vol+loc_bord+raw_iedge;
			loclx_neighdw[ibord][jdir]=iedge;
			//number the elements of the edge
			gx[idir]=(gx[idir]+1+glb_size[idir])%glb_size[idir];
			gx[jdir]=(gx[jdir]-1+glb_size[jdir])%glb_size[jdir];
			glblx_of_edgelx[raw_iedge]=glblx_of_coord(gx);
			gx[idir]=glb_coord_of_loclx[iloc][idir];
			gx[jdir]=glb_coord_of_loclx[iloc][jdir];
		      }
		    //Upper direction
		    if(x[jdir]!=loc_size[jdir]-1 || paral_dir[jdir]==0)
		      {
			int temp_xdir=x[jdir];
			x[jdir]=(x[jdir]+1)%loc_size[jdir];
			int ibord2=loc_vol+bord_offset[idir]+loc_bord/2+bordlx_of_coord(x,idir);
			x[jdir]=temp_xdir;
			loclx_neighup[ibord][jdir]=ibord2;
		      }
		    else //edge
		      {
			int raw_iedge=edge_offset[edge_numb[idir][jdir]]+3*loc_edge/4+edgelx_of_coord(x,idir,jdir);
			//put the value of the neighbour
			int iedge=loc_vol+loc_bord+raw_iedge;
			loclx_neighup[ibord][jdir]=iedge;
			//number the elements of the edge
			gx[idir]=(gx[idir]+1+glb_size[idir])%glb_size[idir];
			gx[jdir]=(gx[jdir]+1+glb_size[jdir])%glb_size[jdir];
			glblx_of_edgelx[raw_iedge]=glblx_of_coord(gx);
			gx[idir]=glb_coord_of_loclx[iloc][idir];
			gx[jdir]=glb_coord_of_loclx[iloc][jdir];
		      }
		  }
	    }
	}
    }
  
  //init sender and receiver points for borders
  for(int i=0;i<4;i++)
    if(paral_dir[i]!=0)
      {
	start_lx_bord_send_up[i]=loclx_of_coord_list(0,0,0,0);
	start_lx_bord_rece_up[i]=(loc_vol+bord_offset[i]+loc_bord/2);
	int x[4];
	for(int jdir=0;jdir<4;jdir++)
	  if(jdir==i) x[jdir]=loc_size[i]-1;
	  else x[jdir]=0;
	start_lx_bord_send_dw[i]=loclx_of_coord(x);
	start_lx_bord_rece_dw[i]=loc_vol+bord_offset[i];
      }
  
  nissa_lx_geom_inited=1;
  master_printf("Cartesian geometry intialized\n");
}

//unset cartesian geometry
void unset_lx_geometry()
{
  if(nissa_lx_geom_inited!=1) crash("cartesian geometry not initialized!");

  master_printf("Unsetting cartesian geometry\n");
  nissa_lx_geom_inited=0;
  
  nissa_free(loc_coord_of_loclx[0]);
  nissa_free(glb_coord_of_loclx[0]);
  nissa_free(loclx_neighup[0]);
  nissa_free(loclx_neighdw[0]);
  
  nissa_free(loc_coord_of_loclx);
  nissa_free(glb_coord_of_loclx);
  nissa_free(loclx_neighup);
  nissa_free(loclx_neighdw);
  
  nissa_free(glblx_of_loclx);
  nissa_free(glblx_of_bordlx);
  nissa_free(loclx_of_bordlx);
  nissa_free(dir_of_bord);
  nissa_free(glblx_of_edgelx);
}

//definitions of lexical ordered sender for borders
void initialize_lx_bord_senders_of_kind(MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *base)
{
  //Various type useful for edges and sub-borders
  MPI_Datatype MPI_3_SLICE;
  MPI_Datatype MPI_23_SLICE;
  MPI_Type_contiguous(loc_size[3],*base,&MPI_3_SLICE);
  MPI_Type_contiguous(loc_size[2]*loc_size[3],*base,&MPI_23_SLICE);

  ///////////define the sender for the 4 kinds of borders////////////
  MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3],*base,&(MPI_BORD_SEND[0]));
  MPI_Type_vector(loc_size[0],1,loc_size[1],MPI_23_SLICE,&(MPI_BORD_SEND[1]));
  MPI_Type_vector(loc_size[0]*loc_size[1],1,loc_size[2],MPI_3_SLICE,&(MPI_BORD_SEND[2]));
  MPI_Type_commit(&(MPI_BORD_SEND[2]));
  MPI_Type_vector(loc_size[0]*loc_size[1]*loc_size[2],1,loc_size[3],*base,&(MPI_BORD_SEND[3]));
  //Commit
  for(int ibord=0;ibord<4;ibord++) MPI_Type_commit(&(MPI_BORD_SEND[ibord]));
}

//definitions of lexical ordered senders for edges
void initialize_lx_edge_senders_of_kind(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *base)
{
  //Various type useful for edges and sub-borders
  MPI_Datatype MPI_3_SLICE;
  MPI_Type_contiguous(loc_size[3],*base,&MPI_3_SLICE);

  ///////////define the sender for the 6 kinds of edges////////////
  //the 01 sender, that is simply a vector of L[2] vector of L[3] 
  MPI_Type_contiguous(loc_size[2]*loc_size[3],*base,&(MPI_EDGE_SEND[0]));
  //the 02 sender is a vector of L[1] segment of length L[3] (already defined) separated by L[2] of them
  MPI_Type_vector(loc_size[1],1,loc_size[2],MPI_3_SLICE,&(MPI_EDGE_SEND[1]));
  //the 03 sender is a vector of length L[1]xL[2] of single elements, separated by L[3] of them
  MPI_Type_vector(loc_size[1]*loc_size[2],1,loc_size[3],*base,&(MPI_EDGE_SEND[2]));
  //the 12 sender should be equal to the 02 sender, with 1->0
  MPI_Type_vector(loc_size[0],1,loc_size[2],MPI_3_SLICE,&(MPI_EDGE_SEND[3]));
  //the 13 sender should be equal to the 03 sender, with 1->0
  MPI_Type_vector(loc_size[0]*loc_size[2],1,loc_size[3],*base,&(MPI_EDGE_SEND[4]));
  //the 23 sender should be equal to the 03 sender with 1<->2
  MPI_Type_vector(loc_size[0]*loc_size[1],1,loc_size[3],*base,&(MPI_EDGE_SEND[5]));
  //Commit
  for(int iedge=0;iedge<6;iedge++) MPI_Type_commit(&(MPI_EDGE_SEND[iedge]));
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

//initalize senders and receivers for borders of lexically ordered vectors
void set_lx_bord_senders_and_receivers(MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *MPI_BORD_RECE,MPI_Datatype *base)
{
  initialize_lx_bord_senders_of_kind(MPI_BORD_SEND,base);
  initialize_lx_bord_receivers_of_kind(MPI_BORD_RECE,base);
}

//initalize senders and receivers for edges of lexically ordered vectors
void set_lx_edge_senders_and_receivers(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base)
{
  initialize_lx_edge_senders_of_kind(MPI_EDGE_SEND,base);
  initialize_lx_edge_receivers_of_kind(MPI_EDGE_RECE,base);
}
