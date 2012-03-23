#pragma once

//Return the index of site of coord x in the border mu,nu
int edgelx_of_coord(int *x,int mu,int nu)
{
  int ilx=0;
  
  for(int rho=0;rho<4;rho++)
    if(rho!=mu && rho!=nu)
      ilx=ilx*loc_size[rho]+x[rho];
  
  return ilx;
}

//Return the index of site of coord x in the border mu
int bordlx_of_coord(int *x,int mu)
{
  int ilx=0;
  
  for(int nu=0;nu<4;nu++) if(nu!=mu) ilx=ilx*loc_size[nu]+x[nu];
  
  return ilx;
}

//wrapper
int bordlx_of_coord_list(int x0,int x1,int x2,int x3,int mu)
{
  coords x={x0,x1,x2,x3};
  return bordlx_of_coord(x,mu);
}

//Return the index of site of coord x in a box of sides s
int lx_of_coord(coords x,coords s)
{
  int ilx=0;
  
  for(int mu=0;mu<4;mu++)
    ilx=ilx*s[mu]+x[mu];
  
  return ilx;
}

//wrappers
int loclx_of_coord(coords x)
{return lx_of_coord(x,loc_size);}
  
//wrappers
int loclx_of_coord_list(int x0,int x1,int x2,int x3)
{
  coords x={x0,x1,x2,x3};
  return lx_of_coord(x,loc_size);
}
  
//wrappers
int glblx_of_coord(coords x)
{return lx_of_coord(x,glb_size);}

//wrappers
int glblx_of_coord_list(int x0,int x1,int x2,int x3)
{
  coords x={x0,x1,x2,x3};
  return lx_of_coord(x,glb_size);
}
  
//Return the coordinate of the rank containing the global coord
void rank_coord_of_site_of_coord(coords rank_coord,coords glb_coord)
{for(int mu=0;mu<4;mu++) rank_coord[mu]=glb_coord[mu]/loc_size[mu];}

//Return the rank of passed coord
int rank_of_coord(coords x)
{return lx_of_coord(x,nrank_dir);}

//Return the rank containing the global coordinates
int rank_hosting_site_of_coord(coords x)
{
  coords p;
  rank_coord_of_site_of_coord(p,x);
  
  return rank_of_coord(p);
}

//Return the local site and rank containing the global coordinates
void get_loclx_and_rank_of_coord(int *ivol,int *rank,coords g)
{
  coords l,p;
  for(int mu=0;mu<4;mu++)
    {
      p[mu]=g[mu]/loc_size[mu];
      l[mu]=g[mu]-p[mu]*loc_size[mu];
    }
  
  (*rank)=rank_of_coord(p);
  (*ivol)=loclx_of_coord(l);
}

//indexes run as t,x,y,z (faster:z)
void set_lx_geometry()
{
  if(nissa_lx_geom_inited==1) crash("cartesian geometry already intialized!");
  
  //find the rank of the neighbour in the various dir
  for(int mu=0;mu<4;mu++)
    MPI_Cart_shift(cart_comm,mu,1,&(rank_neighdw[mu]),&(rank_neighup[mu]));
  
  loc_coord_of_loclx=nissa_malloc("loc_coord_of_loclx",loc_vol,coords);
  glb_coord_of_loclx=nissa_malloc("glb_coord_of_loclx",loc_vol+loc_bord+loc_edge,coords);
  loclx_neighup=nissa_malloc("loclx_neighup",loc_vol+loc_bord,coords);
  loclx_neighdw=nissa_malloc("loclx_neighdw",loc_vol+loc_bord,coords);
  
  //local to global
  glblx_of_loclx=nissa_malloc("glblx_of_loclx",loc_vol,int);
      
  //borders
  glblx_of_bordlx=nissa_malloc("glblx_of_bordlx",loc_bord,int);
  loclx_of_bordlx=nissa_malloc("loclx_of_bordlx",loc_bord,int);
  dir_of_bord=nissa_malloc("dir_of_bord",loc_bord,int);
      
  //edges
  glblx_of_edgelx=nissa_malloc("glblx_of_edgelx",loc_edge,int);
  
  //Label the sites
  coords x,gx={0,0,0,0};
  for(x[0]=0;x[0]<loc_size[0];x[0]++) 
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    for(int mu=0;mu<4;mu++)
	      gx[mu]=x[mu]+rank_coord[mu]*loc_size[mu];
	    
	    //find local and global sites
	    int loc_ind=loclx_of_coord(x);
	    int glb_ind=glblx_of_coord(gx);
	    
	    //save local and global coordinates
	    for(int mu=0;mu<4;mu++)
	      {
		loc_coord_of_loclx[loc_ind][mu]=x[mu];
		glb_coord_of_loclx[loc_ind][mu]=gx[mu];
	      }
	    
	    glblx_of_loclx[loc_ind]=glb_ind;
	  }
  
  //////////////////neighbours search//////////////////////
  
  //now fill the neighbours of sites of the bulk, and the defined
  //neighbours of the sites of the external borders
  nissa_loc_vol_loop(ivol)
    {
      //takes local and global coords of site ivol
      for(int mu=0;mu<4;mu++)
	{
	  x[mu]=loc_coord_of_loclx[ivol][mu];
	  gx[mu]=glb_coord_of_loclx[ivol][mu];
	}
      
      //Direction on the whole iper-cube
      for(int mu=0;mu<4;mu++)
	{
	  //Down direction
	  if(x[mu]!=0 || paral_dir[mu]==0)
	    {
	      int temp_xdir=x[mu];
	      x[mu]=(x[mu]+loc_size[mu]-1)%loc_size[mu];
	      loclx_neighdw[ivol][mu]=loclx_of_coord(x);
	      x[mu]=temp_xdir;
	    }
	  else //border
	    {
	      int raw_ibord=bordlx_of_coord(x,mu)+bord_offset[mu];
	      int ibord=loc_vol+raw_ibord;
	      loclx_neighdw[ivol][mu]=ibord;
	      loclx_neighup[ibord][mu]=ivol;
	      //the movement along the down direction from the border are not defined 
	      loclx_neighdw[ibord][mu]=-1;
	      //number the elements of the border
	      gx[mu]=(gx[mu]-1+glb_size[mu])%glb_size[mu];
	      x[mu]=(x[mu]-1+loc_size[mu])%loc_size[mu];
	      for(int mu=0;mu<4;mu++) glb_coord_of_loclx[ibord][mu]=gx[mu];
	      glblx_of_bordlx[raw_ibord]=glblx_of_coord(gx);
	      loclx_of_bordlx[raw_ibord]=loclx_of_coord(x);
	      dir_of_bord[raw_ibord]=2*mu+1;
	      gx[mu]=glb_coord_of_loclx[ivol][mu];
	      x[mu]=loc_coord_of_loclx[ivol][mu];
	      
	      //This is the bad moment: the movents inside the cube
	      for(int nu=0;nu<4;nu++)
		if(mu!=nu) 
		  {
		    //Down direction
		    if(x[nu]!=0 || paral_dir[nu]==0)
		      {
			int temp_xdir=x[nu];
			x[nu]=(x[nu]+loc_size[nu]-1)%loc_size[nu];
			int ibord2=loc_vol+bord_offset[mu]+bordlx_of_coord(x,mu);
			x[nu]=temp_xdir;
			loclx_neighdw[ibord][nu]=ibord2;
		      }
		    else //i-j- edge
		      {
			int raw_iedge=edge_offset[edge_numb[mu][nu]]+edgelx_of_coord(x,mu,nu);
			//put the value of the neighbour
			int iedge=loc_vol+loc_bord+raw_iedge;
			loclx_neighdw[ibord][nu]=iedge;
			//number the elements of the edge
			gx[mu]=(gx[mu]-1+glb_size[mu])%glb_size[mu];
			gx[nu]=(gx[nu]-1+glb_size[nu])%glb_size[nu];
			glblx_of_edgelx[raw_iedge]=glblx_of_coord(gx);
			for(int rho=0;rho<4;rho++) glb_coord_of_loclx[iedge][rho]=gx[rho];
			gx[mu]=glb_coord_of_loclx[ivol][mu];
			gx[nu]=glb_coord_of_loclx[ivol][nu];
		      }
		    //Upper direction
		    if(x[nu]!=loc_size[nu]-1 || paral_dir[nu]==0)
		      {
			int temp_xdir=x[nu];
			x[nu]=(x[nu]+1)%loc_size[nu];
			int ibord2=loc_vol+bord_offset[mu]+bordlx_of_coord(x,mu);
			x[nu]=temp_xdir;
			loclx_neighup[ibord][nu]=ibord2;
		      }
		    else
		      {
			int raw_iedge=edge_offset[edge_numb[mu][nu]]+edgelx_of_coord(x,mu,nu);
			if(mu<nu) raw_iedge+=loc_edge/4; //edge mu-nu+
			else raw_iedge+=loc_edge/2; //edge nu+mu-
			//put the value of the neighbour
			int iedge=loc_vol+loc_bord+raw_iedge;
			loclx_neighup[ibord][nu]=iedge;
			//number the elements of the edge
			gx[mu]=(gx[mu]-1+glb_size[mu])%glb_size[mu];
			gx[nu]=(gx[nu]+1+glb_size[nu])%glb_size[nu];
			glblx_of_edgelx[raw_iedge]=glblx_of_coord(gx);
			for(int rho=0;rho<4;rho++) glb_coord_of_loclx[iedge][rho]=gx[rho];
			gx[mu]=glb_coord_of_loclx[ivol][mu];
			gx[nu]=glb_coord_of_loclx[ivol][nu];
		      }
		  }
	    }
	  
	  //Upper direction
	  if(x[mu]!=loc_size[mu]-1 || paral_dir[mu]==0)
	    {
	      int temp_xdir=x[mu];
	      x[mu]=(x[mu]+1)%loc_size[mu];
	      loclx_neighup[ivol][mu]=loclx_of_coord(x);
	      x[mu]=temp_xdir;
	    }
	  else //border
	    {
	      int raw_ibord=loc_bord/2+bordlx_of_coord(x,mu)+bord_offset[mu];
	      int ibord=loc_vol+raw_ibord;
	      loclx_neighup[ivol][mu]=ibord;
	      loclx_neighdw[ibord][mu]=ivol;
	      //the movement along the up direction from the border are not defined 
	      loclx_neighup[ibord][mu]=-1;
	      //number the elements of the border
	      gx[mu]=(gx[mu]+1+glb_size[mu])%glb_size[mu];
	      x[mu]=(x[mu]+1+loc_size[mu])%loc_size[mu];
	      for(int mu=0;mu<4;mu++) glb_coord_of_loclx[ibord][mu]=gx[mu];
	      glblx_of_bordlx[raw_ibord]=glblx_of_coord(gx);
	      loclx_of_bordlx[raw_ibord]=loclx_of_coord(x);
	      dir_of_bord[raw_ibord]=2*mu;
	      gx[mu]=glb_coord_of_loclx[ivol][mu];
	      x[mu]=loc_coord_of_loclx[ivol][mu];
	      
	      //Another very bad moment: the movents inside the cube
	      for(int nu=0;nu<4;nu++)
		if(mu!=nu) 
		  {
		    //Down direction
		    if(x[nu]!=0 || paral_dir[nu]==0)
		      {
			int temp_xdir=x[nu];
			x[nu]=(x[nu]+loc_size[nu]-1)%loc_size[nu];
			int ibord2=loc_vol+bord_offset[mu]+loc_bord/2+bordlx_of_coord(x,mu);
			x[nu]=temp_xdir;
			loclx_neighdw[ibord][nu]=ibord2;
		      }
		    else //edge
		      {
			int raw_iedge=edge_offset[edge_numb[mu][nu]]+edgelx_of_coord(x,mu,nu);
			if(mu<nu) raw_iedge+=loc_edge/2; //edge mu-nu+
			else raw_iedge+=loc_edge/4; //edge nu+mu-
			//put the value of the neighbour
			int iedge=loc_vol+loc_bord+raw_iedge;
			loclx_neighdw[ibord][nu]=iedge;
			//number the elements of the edge
			gx[mu]=(gx[mu]+1+glb_size[mu])%glb_size[mu];
			gx[nu]=(gx[nu]-1+glb_size[nu])%glb_size[nu];
			glblx_of_edgelx[raw_iedge]=glblx_of_coord(gx);
			for(int rho=0;rho<4;rho++) glb_coord_of_loclx[iedge][rho]=gx[rho];
			gx[mu]=glb_coord_of_loclx[ivol][mu];
			gx[nu]=glb_coord_of_loclx[ivol][nu];
		      }
		    //Upper direction
		    if(x[nu]!=loc_size[nu]-1 || paral_dir[nu]==0)
		      {
			int temp_xdir=x[nu];
			x[nu]=(x[nu]+1)%loc_size[nu];
			int ibord2=loc_vol+bord_offset[mu]+loc_bord/2+bordlx_of_coord(x,mu);
			x[nu]=temp_xdir;
			loclx_neighup[ibord][nu]=ibord2;
		      }
		    else //edge
		      {
			int raw_iedge=edge_offset[edge_numb[mu][nu]]+3*loc_edge/4+edgelx_of_coord(x,mu,nu);
			//put the value of the neighbour
			int iedge=loc_vol+loc_bord+raw_iedge;
			loclx_neighup[ibord][nu]=iedge;
			//number the elements of the edge
			gx[mu]=(gx[mu]+1+glb_size[mu])%glb_size[mu];
			gx[nu]=(gx[nu]+1+glb_size[nu])%glb_size[nu];
			glblx_of_edgelx[raw_iedge]=glblx_of_coord(gx);
			for(int rho=0;rho<4;rho++) glb_coord_of_loclx[iedge][rho]=gx[rho];
			gx[mu]=glb_coord_of_loclx[ivol][mu];
			gx[nu]=glb_coord_of_loclx[ivol][nu];
		      }
		  }
	    }
	}
    }
  
  //init sender and receiver points for borders
  for(int mu=0;mu<4;mu++)
    if(paral_dir[mu]!=0)
      {
	start_lx_bord_send_up[mu]=loclx_of_coord_list(0,0,0,0);
	start_lx_bord_rece_up[mu]=(loc_vol+bord_offset[mu]+loc_bord/2);
	coords x;
	for(int nu=0;nu<4;nu++)
	  if(nu==mu) x[nu]=loc_size[mu]-1;
	  else x[nu]=0;
	start_lx_bord_send_dw[mu]=loclx_of_coord(x);
	start_lx_bord_rece_dw[mu]=loc_vol+bord_offset[mu];
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
