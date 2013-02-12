#include <math.h>
#include <string.h>

#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../base/debug.h"
#include "../base/routines.h"

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
extern "C" int loclx_of_coord_list(int x0,int x1,int x2,int x3)
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

//combine two points
int glblx_of_comb(int b,int wb,int c,int wc)
{
  coords co;
  for(int mu=0;mu<4;mu++)
    {
      co[mu]=glb_coord_of_loclx[b][mu]*wb+glb_coord_of_loclx[c][mu]*wc;
      while(co[mu]<0) co[mu]+=glb_size[mu];
      co[mu]%=glb_size[mu];
    }

  return glblx_of_coord(co);
}

void glb_coord_of_glblx(coords x,int gx)
{
  for(int mu=3;mu>=0;mu--)
    {
      int next=gx/glb_size[mu];
      x[mu]=gx-next*glb_size[mu];
      gx=next;
    }
}

int glblx_of_diff(int b,int c)
{return glblx_of_comb(b,+1,c,-1);}
  
int glblx_of_summ(int b,int c)
{return glblx_of_comb(b,+1,c,+1);}

int glblx_opp(int b)
{return glblx_of_diff(0,b);}
  
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
//Retrun the rank containing the glblx passed
int rank_hosting_glblx(int gx)
{
  coords c;
  glb_coord_of_glblx(c,gx);
  return rank_hosting_site_of_coord(c);
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

//return the index of the site of passed "pseudolocal" coordinate
//if the coordinates are local, return the index according to the function loclx_of_coord
//if exactly one of the coordinate is just out return its index according to bordlx_of_coord, incremented of previous border and loc_vol
//if exactly two coordinates are outside, return its index according to edgelx_of_coord, incremented as before stated
int full_lx_of_coords(coords ext_x)
{
  int ort_dir_bord[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
  int ort_dir_edge[6][2]={{2,3},{1,3},{1,2},{0,3},{0,2},{0,1}};
  
  //pseudo-localize it
  coords x;
  for(int mu=0;mu<4;mu++)
    {
      x[mu]=ext_x[mu];
      while(x[mu]<0) x[mu]+=glb_size[mu];
      while(x[mu]>=glb_size[mu]) x[mu]-=glb_size[mu];
    }
  
  //check locality
  int isloc=1;
  for(int mu=0;mu<4;mu++)
    {
      isloc&=(x[mu]>=0);
      isloc&=(x[mu]<loc_size[mu]);
    }
  
  if(isloc) return loclx_of_coord(x);
  
  //check borderity
  int isbord[4];
  for(int mu=0;mu<4;mu++)
    {
      isbord[mu]=0;
      if(paral_dir[mu])
	{
	  if(x[mu]==glb_size[mu]-1) isbord[mu]=-1;
	  if(x[mu]==loc_size[mu]) isbord[mu]=+1;
	}
    }
  
  //check if it is in one of the 4 forward or backward borders
  for(int mu=0;mu<4;mu++)
    if((isbord[ort_dir_bord[mu][0]]==0)&&(isbord[ort_dir_bord[mu][1]]==0)&&(isbord[ort_dir_bord[mu][2]]==0))
      {
	if(isbord[mu]==-1) return loc_vol+bord_offset[mu]+bordlx_of_coord(x,mu);             //backward border comes first
	if(isbord[mu]==+1) return loc_vol+bord_vol/2+bord_offset[mu]+bordlx_of_coord(x,mu);  //forward border comes after
      }
  
  //check if it is in one of the 6 --,-+,+-,++ edges
  for(int mu=0;mu<4;mu++)
    for(int inu=0;inu<3;inu++)
      {
	int nu=ort_dir_bord[mu][inu];
	int iedge=edge_numb[mu][nu];
	
	//order mu,nu
	int al=(mu<nu)?mu:nu;
	int be=(mu>nu)?mu:nu;
	
	if((isbord[ort_dir_edge[iedge][0]]==0)&&(isbord[ort_dir_edge[iedge][1]]==0))
	  {
	    if((isbord[al]==-1)&&(isbord[be]==-1)) return loc_vol+bord_vol+edge_offset[iedge]+0*edge_vol/4+edgelx_of_coord(x,mu,nu);
	    if((isbord[al]==-1)&&(isbord[be]==+1)) return loc_vol+bord_vol+edge_offset[iedge]+1*edge_vol/4+edgelx_of_coord(x,mu,nu);
	    if((isbord[al]==+1)&&(isbord[be]==-1)) return loc_vol+bord_vol+edge_offset[iedge]+2*edge_vol/4+edgelx_of_coord(x,mu,nu);
	    if((isbord[al]==+1)&&(isbord[be]==+1)) return loc_vol+bord_vol+edge_offset[iedge]+3*edge_vol/4+edgelx_of_coord(x,mu,nu);
	  }
      }
  
  return -1;
}

//return the border site of opposite sites at surface
int bordlx_of_surflx(int loclx,int mu)
{
  if(!paral_dir[mu]) return -1;
  if(loc_size[mu]<2) crash("not working if one dir is smaller than 2");
  
  //copy the coords
  coords x={loc_coord_of_loclx[loclx][0],loc_coord_of_loclx[loclx][1],loc_coord_of_loclx[loclx][2],loc_coord_of_loclx[loclx][3]};
  
  //dw surf
  if(x[mu]==0) 
    {
      x[mu]=loc_size[mu];
      return full_lx_of_coords(x);
    }

  //up surf
  if(x[mu]==loc_size[mu]-1)
    {
      x[mu]=-1;
      return full_lx_of_coords(x);
    }
  
  return -1;
}

//wrapper
extern "C" int full_lx_of_coords_list(const int t,const int x,const int y,const int z)
{
  coords c={t,x,y,z};
  return full_lx_of_coords(c);
}

//label all the sites: bulk, border and edge
void label_all_sites()
{
  coords x;
  
  for(x[0]=-paral_dir[0];x[0]<loc_size[0]+paral_dir[0];x[0]++) 
    for(x[1]=-paral_dir[1];x[1]<loc_size[1]+paral_dir[1];x[1]++) 
      for(x[2]=-paral_dir[2];x[2]<loc_size[2]+paral_dir[2];x[2]++) 
	for(x[3]=-paral_dir[3];x[3]<loc_size[3]+paral_dir[3];x[3]++) 
	  {
	    //check if it is defined
	    int iloc=full_lx_of_coords(x);
	    if(iloc!=-1)
	      {
		//compute global coordinates, assigning
		for(int nu=0;nu<4;nu++)
		  glb_coord_of_loclx[iloc][nu]=(x[nu]+rank_coord[nu]*loc_size[nu]+glb_size[nu])%glb_size[nu];
		
		//find the global index
		int iglb=glblx_of_coord(glb_coord_of_loclx[iloc]);
		
		//if it is on the bulk store it
		if(iloc<loc_vol)
		  {
		    for(int nu=0;nu<4;nu++) loc_coord_of_loclx[iloc][nu]=x[nu];
		    glblx_of_loclx[iloc]=iglb;
		  }
		
		//if it is on the border store it
		if(iloc>=loc_vol&&iloc<loc_vol+bord_vol)
		  {
		    int ibord=iloc-loc_vol;
		    glblx_of_bordlx[ibord]=iglb;
		    loclx_of_bordlx[ibord]=iloc;
		  }
		
		//if it is on the edge store it
		if(iloc>=loc_vol+bord_vol)
		  {
		    int iedge=iloc-loc_vol-bord_vol;
		    glblx_of_edgelx[iedge]=iglb;
		    //loclx_of_edgelx[iedge]=iloc;
		  }
	      }
	  }
}

//find the neighbours
void find_neighbouring_sites()
{
  //loop over the four directions
  for(int ivol=0;ivol<loc_vol+bord_vol+edge_vol;ivol++)
    for(int mu=0;mu<4;mu++)
      {
	//copy the coords
	coords n;
	for(int nu=0;nu<4;nu++) n[nu]=glb_coord_of_loclx[ivol][nu]-loc_size[nu]*rank_coord[nu];
	
	//move forward
	n[mu]++;
	int nup=full_lx_of_coords(n);
	//move backward
	n[mu]-=2;
	int ndw=full_lx_of_coords(n);
	
	//if "local" assign it (automatically -1 otherwise)
	loclx_neighup[ivol][mu]=nup;
	loclx_neighdw[ivol][mu]=ndw;
      }
}

//finds how to fill the borders with opposite surface (up b->dw s)
void find_surf_of_bord()
{
  nissa_loc_vol_loop(loclx)
    for(int mu=0;mu<4;mu++)
      {
        int bordlx=bordlx_of_surflx(loclx,mu);
        if(bordlx!=-1) surflx_of_bordlx[bordlx]=loclx;
      }
}

//indexes run as t,x,y,z (faster:z)
void set_lx_geometry()
{
  if(nissa_lx_geom_inited==1) crash("cartesian geometry already intialized!");
  nissa_lx_geom_inited=1;
  
  if(nissa_grid_inited!=1) crash("grid not initialized!");
  
  //find the rank of the neighbour in the various dir
  for(int mu=0;mu<4;mu++)
    MPI_Cart_shift(cart_comm,mu,1,&(rank_neighdw[mu]),&(rank_neighup[mu]));
  memcpy(rank_neigh[0],rank_neighdw,sizeof(coords));
  memcpy(rank_neigh[1],rank_neighup,sizeof(coords));

  loc_coord_of_loclx=nissa_malloc("loc_coord_of_loclx",loc_vol,coords);
  glb_coord_of_loclx=nissa_malloc("glb_coord_of_loclx",loc_vol+bord_vol+edge_vol,coords);
  loclx_neigh[0]=loclx_neighdw=nissa_malloc("loclx_neighdw",loc_vol+bord_vol+edge_vol,coords);
  loclx_neigh[1]=loclx_neighup=nissa_malloc("loclx_neighup",loc_vol+bord_vol+edge_vol,coords);  
  ignore_borders_communications_warning(loc_coord_of_loclx);
  ignore_borders_communications_warning(glb_coord_of_loclx);
  ignore_borders_communications_warning(loclx_neighup);
  ignore_borders_communications_warning(loclx_neighdw);

  //local to global
  glblx_of_loclx=nissa_malloc("glblx_of_loclx",loc_vol,int);
      
  //borders
  glblx_of_bordlx=nissa_malloc("glblx_of_bordlx",bord_vol,int);
  loclx_of_bordlx=nissa_malloc("loclx_of_bordlx",bord_vol,int);
  surflx_of_bordlx=nissa_malloc("surflx_of_bordlx",bord_vol,int);
  
  //edges
  glblx_of_edgelx=nissa_malloc("glblx_of_edgelx",edge_vol,int);
  
  //label the sites and neighbours
  label_all_sites();
  find_neighbouring_sites();
  
  //matches surface and opposite border
  find_surf_of_bord();
  
  //init sender and receiver points for borders
  for(int mu=0;mu<4;mu++)
    if(paral_dir[mu]!=0)
      {
	start_lx_bord_send_up[mu]=loclx_of_coord_list(0,0,0,0);
	start_lx_bord_rece_up[mu]=(loc_vol+bord_offset[mu]+bord_vol/2);
	coords x;
	for(int nu=0;nu<4;nu++)
	  if(nu==mu) x[nu]=loc_size[mu]-1;
	  else x[nu]=0;
	start_lx_bord_send_dw[mu]=loclx_of_coord(x);
	start_lx_bord_rece_dw[mu]=loc_vol+bord_offset[mu];
      }
  
  master_printf("Cartesian geometry intialized\n");
}

//global movements
int glblx_neighup(int gx,int mu)
{
  coords c;
  glb_coord_of_glblx(c,gx);
  c[mu]=(c[mu]+1)%glb_size[mu];
  
  return glblx_of_coord(c);
}
int glblx_neighdw(int gx,int mu)
{
  coords c;
  glb_coord_of_glblx(c,gx);
  c[mu]=(c[mu]+glb_size[mu]-1)%glb_size[mu];
  
  return glblx_of_coord(c);
}

//wrapper for a previous defined function
void get_loclx_and_rank_of_glblx(int *lx,int *rx,int gx)
{
  coords c;
  glb_coord_of_glblx(c,gx);
  get_loclx_and_rank_of_coord(lx,rx,c);
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
  nissa_free(surflx_of_bordlx);
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

//define all the local lattice momenta
void define_local_momenta(momentum_t *k,double *k2,momentum_t *ktilde,double *ktilde2,momentum_t bc)
{
  if(!nissa_lx_geom_inited) set_lx_geometry();
  
  //first of all, defines the local momenta for the various directions
  nissa_loc_vol_loop(imom)
    {
      k2[imom]=ktilde2[imom]=0;
      for(int mu=0;mu<4;mu++)
	{
	  k[imom][mu]=M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu];
	  ktilde[imom][mu]=sin(k[imom][mu]);
	  
	  k2[imom]+=k[imom][mu]*k[imom][mu];
	  ktilde2[imom]+=ktilde[imom][mu]*ktilde[imom][mu];
	}
    }
}
