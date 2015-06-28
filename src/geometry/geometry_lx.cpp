#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#ifdef SPI
 #include <malloc.h>
 #include <stdlib.h>
#endif

namespace nissa
{
  //Return the index of site of coord x in the border mu,nu
  int edgelx_of_coord(int *x,int mu,int nu)
  {
    int ilx=0;
    
    for(int rho=0;rho<NDIM;rho++)
      if(rho!=mu && rho!=nu)
	ilx=ilx*loc_size[rho]+x[rho];
    
    return ilx;
  }
  
  //Return the index of site of coord x in the border mu
  int bordlx_of_coord(int *x,int mu)
  {
    int ilx=0;  
    for(int nu=0;nu<NDIM;nu++)
      if(nu!=mu)
	ilx=ilx*loc_size[nu]+x[nu];
    
    return ilx;
  }
  
  //Return the index of site of coord x in a box of sides s
  int lx_of_coord(coords x,coords s)
  {
    int ilx=0;
    
    for(int mu=0;mu<NDIM;mu++)
      ilx=ilx*s[mu]+x[mu];
    
    return ilx;
  }
  void coord_of_lx(coords x,int ilx,coords s)
  {
    for(int mu=NDIM-1;mu>=0;mu--)
      {
	x[mu]=ilx%s[mu];
	ilx/=s[mu];
      }
  }
  
  //wrappers
  int loclx_of_coord(coords x)
  {return lx_of_coord(x,loc_size);}
  
  //wrappers
  int glblx_of_coord(coords x)
  {return lx_of_coord(x,glb_size);}
  
  //combine two points
  int glblx_of_comb(int b,int wb,int c,int wc)
  {
    coords co;
    for(int mu=0;mu<NDIM;mu++)
      {
	co[mu]=glb_coord_of_loclx[b][mu]*wb+glb_coord_of_loclx[c][mu]*wc;
	while(co[mu]<0) co[mu]+=glb_size[mu];
	co[mu]%=glb_size[mu];
      }
    
    return glblx_of_coord(co);
  }
  
  void glb_coord_of_glblx(coords x,int gx)
  {
    for(int mu=NDIM-1;mu>=0;mu--)
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
  {for(int mu=0;mu<NDIM;mu++) rank_coord[mu]=glb_coord[mu]/loc_size[mu];}
  
  //Return the rank of passed coord
  int rank_of_coord(coords x)
  {return lx_of_coord(x,nrank_dir);}
  void coord_of_rank(coords c,int x)
  {coord_of_lx(c,x,nrank_dir);}
  
  //Return the rank containing the global coordinates
  int rank_hosting_site_of_coord(coords x)
  {
    coords p;
    rank_coord_of_site_of_coord(p,x);
    
    return rank_of_coord(p);
  }
  //Return the rank containing the glblx passed
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
    for(int mu=0;mu<NDIM;mu++)
      {
	p[mu]=g[mu]/loc_size[mu];
	l[mu]=g[mu]-p[mu]*loc_size[mu];
      }
    
    (*rank)=rank_of_coord(p);
    (*ivol)=loclx_of_coord(l);
  }
  
  //Return the global index of site addressed by rank and loclx
  int get_glblx_of_rank_and_loclx(int irank,int loclx)
  {
    coords p;
    coord_of_rank(p,irank);
    
    int iglblx=0;
    for(int mu=0;mu<NDIM;mu++)
      iglblx=iglblx*glb_size[mu]+loc_coord_of_loclx[loclx][mu];
    
    return iglblx;
  }
  
  //return the index of the site of passed "pseudolocal" coordinate
  //if the coordinates are local, return the index according to the function loclx_of_coord
  //if exactly one of the coordinate is just out return its index according to bordlx_of_coord, incremented of previous border and loc_vol
  //if exactly two coordinates are outside, return its index according to edgelx_of_coord, incremented as before stated
  int full_lx_of_coords(coords ext_x)
  {
    //pseudo-localize it
    coords x;
    for(int mu=0;mu<NDIM;mu++)
      {
	x[mu]=ext_x[mu];
	while(x[mu]<0) x[mu]+=glb_size[mu];
	while(x[mu]>=glb_size[mu]) x[mu]-=glb_size[mu];
      }
    
    //check locality
    int isloc=1;
    for(int mu=0;mu<NDIM;mu++)
      {
	isloc&=(x[mu]>=0);
	isloc&=(x[mu]<loc_size[mu]);
      }
    
    if(isloc) return loclx_of_coord(x);
    
    //check borderity
    coords is_bord;
    for(int mu=0;mu<NDIM;mu++)
      {
	is_bord[mu]=0;
	if(paral_dir[mu])
	  {
	    if(x[mu]==glb_size[mu]-1) is_bord[mu]=-1;
	    if(x[mu]==loc_size[mu]) is_bord[mu]=+1;
	  }
      }
    
    //check if it is in one of the NDIM forward or backward borders
    for(int mu=0;mu<NDIM;mu++)
      {
	bool is=is_bord[mu];
	for(int inu=0;inu<NDIM-1;inu++) is&=(is_bord[perp_dir[mu][inu]]==0);
	
	if(is)
	  {
	    if(is_bord[mu]==-1) return loc_vol+bord_offset[mu]+bordlx_of_coord(x,mu);             //backward border comes first
	    if(is_bord[mu]==+1) return loc_vol+bord_vol/2+bord_offset[mu]+bordlx_of_coord(x,mu);  //forward border comes after
	    crash("if is bord should not arrive here %d %d %d %d",ext_x[0],ext_x[1],ext_x[2],ext_x[3]);
	  }
      }
    
    //check if it is in one of the NDIM*(NDIM-1)/2 --,-+,+-,++ edges
    for(int mu=0;mu<NDIM;mu++)
      for(int inu=0;inu<NDIM-1;inu++)
	{
	  int nu=perp_dir[mu][inu];
	  
	  //order mu,nu
	  int al=(mu<nu)?mu:nu;
	  int be=(mu>nu)?mu:nu;
	  
	  bool is=is_bord[mu]&&is_bord[nu];
#if NDIM>=3
	  for(int irho=0;irho<NDIM-2;irho++) is&=(is_bord[perp2_dir[mu][inu][irho]]==0);
#endif
	  
	  if(is)
	    {
	      int iedge=edge_numb[mu][nu];
	      if((is_bord[al]==-1)&&(is_bord[be]==-1)) return loc_vol+bord_vol+edge_offset[iedge]+0*edge_vol/4+edgelx_of_coord(x,mu,nu);
	      if((is_bord[al]==-1)&&(is_bord[be]==+1)) return loc_vol+bord_vol+edge_offset[iedge]+1*edge_vol/4+edgelx_of_coord(x,mu,nu);
	      if((is_bord[al]==+1)&&(is_bord[be]==-1)) return loc_vol+bord_vol+edge_offset[iedge]+2*edge_vol/4+edgelx_of_coord(x,mu,nu);
	      if((is_bord[al]==+1)&&(is_bord[be]==+1)) return loc_vol+bord_vol+edge_offset[iedge]+3*edge_vol/4+edgelx_of_coord(x,mu,nu);
	      crash("Edge: %d, mu=%d, nu=%d %d %d %d %d",iedge,mu,nu,ext_x[0],ext_x[1],ext_x[2],ext_x[3]);
	    }
	}
    
    return -1;
  }
  
  //return the border site adiacent at surface
  int bordlx_of_surflx(int loclx,int mu)
  {
    if(!paral_dir[mu]) return -1;
    if(loc_size[mu]<2) crash("not working if one dir is smaller than 2");
    
    if(loc_coord_of_loclx[loclx][mu]==0) return loclx_neighdw[loclx][mu]-loc_vol;
    if(loc_coord_of_loclx[loclx][mu]==loc_size[mu]-1) return loclx_neighup[loclx][mu]-loc_vol;
    
    return -1;
  }
  
  //label all the sites: bulk, border and edge
  void label_all_sites()
  {
    //defined a box extending over borders/edges
    int extended_box_vol=1;
    coords extended_box_size;
    for(int mu=0;mu<NDIM;mu++)
      {
	extended_box_size[mu]=paral_dir[mu]*2+loc_size[mu];
	extended_box_vol*=extended_box_size[mu];
      }
    
    for(int ivol=0;ivol<extended_box_vol;ivol++)
      {
	//subtract by one if dir is parallelized
	coords x;
	coord_of_lx(x,ivol,extended_box_size);
	for(int mu=0;mu<NDIM;mu++) if(paral_dir[mu]) x[mu]--;
	
	//check if it is defined
	int iloc=full_lx_of_coords(x);
	if(iloc!=-1)
	  {
	    //compute global coordinates, assigning
	    for(int nu=0;nu<NDIM;nu++)
	      glb_coord_of_loclx[iloc][nu]=(x[nu]+rank_coord[nu]*loc_size[nu]+glb_size[nu])%glb_size[nu];
	    
	    //find the global index
	    int iglb=glblx_of_coord(glb_coord_of_loclx[iloc]);
	    
	    //if it is on the bulk store it
	    if(iloc<loc_vol)
	      {
		for(int nu=0;nu<NDIM;nu++) loc_coord_of_loclx[iloc][nu]=x[nu];
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
      for(int mu=0;mu<NDIM;mu++)
	{
	  //copy the coords
	  coords n;
	  for(int nu=0;nu<NDIM;nu++) n[nu]=glb_coord_of_loclx[ivol][nu]-loc_size[nu]*rank_coord[nu];
	  
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
    NISSA_LOC_VOL_LOOP(loclx)
      for(int mu=0;mu<NDIM;mu++)
	{
	  int bordlx=bordlx_of_surflx(loclx,mu);
	  if(bordlx!=-1) surflx_of_bordlx[bordlx]=loclx;
	}
  }
  
  //index all the sites on bulk
  void find_bulk_sites()
  {
    //check surfacity
    int ibulk=0,inon_fw_surf=0,inon_bw_surf=0;
    int isurf=0,ifw_surf=0,ibw_surf=0;
    NISSA_LOC_VOL_LOOP(ivol)
    {
      //find if it is on bulk or non_fw or non_bw surf
      int is_bulk=true,is_non_fw_surf=true,is_non_bw_surf=true;
      for(int mu=0;mu<NDIM;mu++)
	if(paral_dir[mu])
	  {
	    if(loc_coord_of_loclx[ivol][mu]==loc_size[mu]-1) is_bulk=is_non_fw_surf=false;
	    if(loc_coord_of_loclx[ivol][mu]==0)              is_bulk=is_non_bw_surf=false;
	  }
      
      //mark it
      if(is_bulk) loclx_of_bulklx[ibulk++]=ivol;
      else        loclx_of_surflx[isurf++]=ivol;
      if(is_non_fw_surf) loclx_of_non_fw_surflx[inon_fw_surf++]=ivol;
      else               loclx_of_fw_surflx[ifw_surf++]=ivol;
      if(is_non_bw_surf) loclx_of_non_bw_surflx[inon_bw_surf++]=ivol;
      else               loclx_of_bw_surflx[ibw_surf++]=ivol;
    }
    
    if(ibulk!=bulk_vol) crash("mismatch in bulk id");
    if(isurf!=surf_vol) crash("mismatch in surf id");
    if(inon_fw_surf!=non_fw_surf_vol) crash("mismatch in non_fw_surf id");
    if(inon_bw_surf!=non_bw_surf_vol) crash("mismatch in non_bw_surf id");
    if(ifw_surf!=fw_surf_vol) crash("mismatch in fw_surf id");
    if(ibw_surf!=bw_surf_vol) crash("mismatch in bw_surf id");
  }  
  
  //indexes run as t,x,y,z (faster:z)
  void set_lx_geometry()
  {
    if(lx_geom_inited==1) crash("cartesian geometry already intialized!");
    lx_geom_inited=1;
    
    if(grid_inited!=1) crash("grid not initialized!");
    
    //find the rank of the neighbour in the various dir
    for(int mu=0;mu<NDIM;mu++)
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
    
    //bulk and surfs
    loclx_of_bulklx=nissa_malloc("loclx_of_bulklx",bulk_vol,int);
    loclx_of_surflx=nissa_malloc("loclx_of_surflx",surf_vol,int);
    loclx_of_non_bw_surflx=nissa_malloc("loclx_of_non_bw_surflx",non_bw_surf_vol,int);
    loclx_of_non_fw_surflx=nissa_malloc("loclx_of_non_fw_surflx",non_fw_surf_vol,int);
    loclx_of_bw_surflx=nissa_malloc("loclx_of_bw_surflx",bw_surf_vol,int);
    loclx_of_fw_surflx=nissa_malloc("loclx_of_fw_surflx",fw_surf_vol,int);
    
    //edges
    glblx_of_edgelx=nissa_malloc("glblx_of_edgelx",edge_vol,int);
    
    //label the sites and neighbours
    label_all_sites();
    find_neighbouring_sites();
    
    //matches surface and opposite border
    find_surf_of_bord();
    
    //find bulk sites
    find_bulk_sites();
    
    //allocate a buffer large enough to allow communications of su3spinspin lx border
    recv_buf_size=std::max(recv_buf_size,bord_vol*sizeof(su3spinspin));
    send_buf_size=std::max(send_buf_size,bord_vol*sizeof(su3spinspin));
    
    //create the sweepers but do not fully initialize
    Wilson_sweeper=new gauge_sweeper_t;
    tlSym_sweeper=new gauge_sweeper_t;
    
    //set locd geom (one of the dimension local and fastest running, the other as usual)
    max_locd_size=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	remap_lx_to_locd[mu]=remap_locd_to_lx[mu]=NULL;
	max_locd_perp_size_per_dir[mu]=(glb_vol/glb_size[mu]+nranks-1)/nranks;
	max_locd_size=std::max(max_locd_size,max_locd_perp_size_per_dir[mu]*glb_size[mu]);
	locd_perp_size_per_dir[mu]=(int)std::min((int64_t)max_locd_perp_size_per_dir[mu],glb_vol-max_locd_perp_size_per_dir[mu]*rank);
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
    if(lx_geom_inited!=1) crash("cartesian geometry not initialized!");
    
    master_printf("Unsetting cartesian geometry\n");
    lx_geom_inited=0;
    
#if defined BGQ && defined SPI
    free(recv_buf);
    free(send_buf);
#else
    nissa_free(recv_buf);
    nissa_free(send_buf);
#endif
    
    nissa_free(loc_coord_of_loclx);
    nissa_free(glb_coord_of_loclx);
    nissa_free(loclx_neighup);
    nissa_free(loclx_neighdw);
    
    nissa_free(glblx_of_loclx);
    nissa_free(glblx_of_bordlx);
    nissa_free(loclx_of_bordlx);
    nissa_free(surflx_of_bordlx);
    nissa_free(glblx_of_edgelx);
    
    nissa_free(loclx_of_bulklx);
    nissa_free(loclx_of_surflx);
    nissa_free(loclx_of_non_fw_surflx);
    nissa_free(loclx_of_fw_surflx);
    nissa_free(loclx_of_non_bw_surflx);
    nissa_free(loclx_of_bw_surflx);
    
    delete Wilson_sweeper;
    delete tlSym_sweeper;
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
  
  //definitions of lexical ordered receivers for edges
  void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base)
  {
    //define the NDIM*(NDIM-1)/2 edges receivers, which are contiguous in memory
    int iedge=0;
    for(int mu=0;mu<NDIM;mu++)
      for(int nu=mu+1;nu<NDIM;nu++)
	{
	  MPI_Type_contiguous(loc_vol/loc_size[mu]/loc_size[nu],*base,&(MPI_EDGE_RECE[iedge]));
	  MPI_Type_commit(&(MPI_EDGE_RECE[iedge]));
	  iedge++;
	}
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
    if(!lx_geom_inited) set_lx_geometry();
    
    //first of all, defines the local momenta for the various directions
    NISSA_LOC_VOL_LOOP(imom)
    {
      k2[imom]=ktilde2[imom]=0;
      for(int mu=0;mu<NDIM;mu++)
	{
	  k[imom][mu]=M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu];
	  ktilde[imom][mu]=sin(k[imom][mu]);
	  
	  k2[imom]+=k[imom][mu]*k[imom][mu];
	  ktilde2[imom]+=ktilde[imom][mu]*ktilde[imom][mu];
	}
    }
  }
  
  //multiply the whole conf for stag phases
  THREADABLE_FUNCTION_1ARG(addrem_stagphases_to_lx_conf, quad_su3*,lx_conf)
  {
    //we must ensure that nobody is using the conf
    THREAD_BARRIER();
    
    //work also on borders and edges if allocated and valid
    int ending=loc_vol;
    if(check_borders_allocated(lx_conf) && check_borders_valid(lx_conf)) ending+=bord_vol;
    if(check_edges_allocated(lx_conf) && check_edges_valid(lx_conf)) ending+=edge_vol;
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,ending)
      {
	int d=0;
	
	//phase in direction 1 is always 0 so nothing has to be done in that dir
	//if(d%2==1) su3_prod_double(lx_conf[ivol][1],lx_conf[ivol][1],-1);
	
	//direction 2
	d+=glb_coord_of_loclx[ivol][1];
	if(d%2==1) su3_prod_double(lx_conf[ivol][2],lx_conf[ivol][2],-1);
	
	//direction 3
	d+=glb_coord_of_loclx[ivol][2];
	if(d%2==1) su3_prod_double(lx_conf[ivol][3],lx_conf[ivol][3],-1);
	
	//direction 0
	d+=glb_coord_of_loclx[ivol][3];
	//debug: putting the anti-periodic condition on the temporal border
	//in future remove it!!!
	if(glb_coord_of_loclx[ivol][0]==glb_size[0]-1) d+=1;
	if(d%2==1) su3_prod_double(lx_conf[ivol][0],lx_conf[ivol][0],-1);
      }
    
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //check that passed argument is between 0 and 15
  inline void crash_if_not_hypercubic_red(int hyp_red)
  {if(hyp_red<0||hyp_red>=16) crash("%d not a hyperucbic reduced point",hyp_red);}
  
  //return the coordinates inside the hypercube
  void red_coords_of_hypercubic_red_point(coords h,int hyp_red)
  {
    crash_if_not_hypercubic_red(hyp_red);
    
    for(int mu=NDIM-1;mu>=0;mu--)
      {
	h[mu]=hyp_red%2;
	hyp_red/=2;
      }
  }
  
  //takes the NDIM coordinates of the hypercube vertex one by one
  void lx_coords_of_hypercube_vertex(coords lx,int hyp_cube)
  {
    for(int mu=NDIM-1;mu>=0;mu--)
      {
	lx[mu]=2*(hyp_cube%(loc_size[mu]/2));
	hyp_cube/=loc_size[mu]/2;
      }
  }
  
  //return the point of passed coords in the hypercube
  int hypercubic_red_point_of_red_coords(coords h)
  {
    int hyp=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	if(h[mu]<0||h[mu]>=2) crash("coordinate %d not in the range [0,1]",h[mu]);
	hyp*=2;
	hyp+=h[mu];
      }
    
    return hyp;
  }
}
