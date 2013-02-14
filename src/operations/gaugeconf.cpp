#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>
#include <string.h>

#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../base/debug.h"
#include "../base/communicate.h"
#include "../base/random.h"
#include "../geometry/geometry_mix.h"
#include "../new_types/complex.h"
#include "../new_types/su3.h"
#include "../operations/remap_vector.h"
#include "../operations/su3_paths/plaquette.h"
#include "../routines/ios.h"

/*
  rotate a field anti-clockwise by 90 degrees

   0---1---2---0        0---6---3---0     
   |   |   |   |        |   |   |   |            	
   6---7---8---6        2---8---5---2
   |   |   |   |        |   |   |   |     	
   3---4---5---3        1---7---4---1     
   |   |   |   |        |   |   |   |     	
   O---1---2---0        O---6---3---0     

   d2
   O d1

   where d1=axis+1
   and   d2=d1+1
   
*/

void ac_rotate_vector(void *out,void *in,int axis,int bps)
{
  //find the two swapping direction 
  int d1=1+(axis-1+1)%3;
  int d2=1+(axis-1+2)%3;
  
  //check that the two directions have the same size and that we are not asking 0 as axis
  if(glb_size[d1]!=glb_size[d2]) crash("Rotation works only if dir %d and %d have the same size!",glb_size[d1],glb_size[d2]);
  if(axis==0) crash("Error, only spatial rotations implemented");
  int L=glb_size[d1];
  
  //allocate destinations and sources
  coords *xto=nissa_malloc("xto",loc_vol,coords);
  coords *xfr=nissa_malloc("xfr",loc_vol,coords);
  
  //scan all local sites to see where to send and from where to expect data
  nissa_loc_vol_loop(ivol)
    {
      //copy 0 and axis coord to "to" and "from" sites
      xto[ivol][0]=xfr[ivol][0]=glb_coord_of_loclx[ivol][0];
      xto[ivol][axis]=xfr[ivol][axis]=glb_coord_of_loclx[ivol][axis];
      
      //find reamining coord of "to" site
      xto[ivol][d1]=(L-glb_coord_of_loclx[ivol][d2])%L;
      xto[ivol][d2]=glb_coord_of_loclx[ivol][d1];
      
      //find remaining coord of "from" site
      xfr[ivol][d1]=glb_coord_of_loclx[ivol][d2];
      xfr[ivol][d2]=(L-glb_coord_of_loclx[ivol][d1])%L;
    }
  
  //call the remapping
  remap_vector((char*)out,(char*)in,xto,xfr,bps);
  
  //free vectors
  nissa_free(xfr);
  nissa_free(xto);
}


/*
  rotate the gauge configuration anti-clockwise by 90 degrees
  this is more complicated than a single vector because of link swaps
  therefore the rotation is accomplished through 2 separates steps
  
   .---.---.---.     .---.---.---.       .---.---.---.     
   |           |     |           |       |           |            	
   .   B 3 C   .     .   B 2'C   .       .   C 4'D   .     
   |   2   4   |     |   1   3   |       |   3   1   |     	
   .   A 1 D   .     .   A 4'D   .       .   B 2'A   .     
   |           |     |           |       |           |     	
   O---.---.---.     O---.---.---.       O---.---.---.     

   d2
   O d1

*/

void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis)
{
  int d0=0;
  int d1=1+(axis-1+1)%3;
  int d2=1+(axis-1+2)%3;
  int d3=axis;
  
  //allocate a temporary conf with borders
  quad_su3 *temp_conf=nissa_malloc("temp_conf",loc_vol+bord_vol,quad_su3);
  memcpy(temp_conf,in,loc_vol*sizeof(quad_su3));
  communicate_lx_quad_su3_borders(temp_conf);
  
  //now reorder links
  nissa_loc_vol_loop(ivol)
      {
	//copy temporal direction and axis
	memcpy(out[ivol][d0],temp_conf[ivol][d0],sizeof(su3));
	memcpy(out[ivol][d3],temp_conf[ivol][d3],sizeof(su3));
	//swap the other two
	unsafe_su3_hermitian(out[ivol][d1],temp_conf[loclx_neighdw[ivol][d2]][d2]);
	memcpy(out[ivol][d2],temp_conf[ivol][d1],sizeof(su3));
      }
  
  //rotate rigidly
  ac_rotate_vector(out,out,axis,sizeof(quad_su3));
}

//put boundary conditions on the gauge conf
void put_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges)
{
  complex theta[4];
  for(int idir=0;idir<4;idir++)
    {
      theta[idir][0]=cos(theta_in_pi[idir]*M_PI/glb_size[idir]);
      theta[idir][1]=sin(theta_in_pi[idir]*M_PI/glb_size[idir]);
    }
  
  int nsite=loc_vol;
  if(putonbords) nsite+=bord_vol;
  if(putonedges) nsite+=edge_vol;

  for(int ivol=0;ivol<nsite;ivol++)
    for(int idir=0;idir<4;idir++) safe_su3_prod_complex(conf[ivol][idir],conf[ivol][idir],theta[idir]);
  
  if(!putonbords) set_borders_invalid(conf);
  if(!putonedges) set_edges_invalid(conf);
}

void rem_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges)
{
  double meno_theta_in_pi[4]={-theta_in_pi[0],-theta_in_pi[1],-theta_in_pi[2],-theta_in_pi[3]};
  put_boundaries_conditions(conf,meno_theta_in_pi,putonbords,putonedges);
}

//Adapt the border condition
void adapt_theta(quad_su3 *conf,double *old_theta,double *put_theta,int putonbords,int putonedges)
{
  double diff_theta[4];
  int adapt=0;

  for(int idir=0;idir<4;idir++)
    {
      adapt=adapt || (old_theta[idir]!=put_theta[idir]);
      diff_theta[idir]=put_theta[idir]-old_theta[idir];
      old_theta[idir]=put_theta[idir];
    }
  
  if(adapt)
    {
      master_printf("Necessary to add boundary condition: %f %f %f %f\n",diff_theta[0],diff_theta[1],diff_theta[2],diff_theta[3]);
      put_boundaries_conditions(conf,diff_theta,putonbords,putonedges);
    }
}

//generate an identical conf
void generate_cold_eo_conf(quad_su3 **conf)
{
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
        for(int mu=0;mu<4;mu++)
          su3_put_to_id(conf[par][ivol][mu]);

      set_borders_invalid(conf[par]);
    }
}

//generate a random conf
void generate_hot_eo_conf(quad_su3 **conf)
{
  if(nissa_loc_rnd_gen_inited==0) crash("random number generator not inited");
  
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ieo)
        {
	  int ilx=loclx_of_loceo[par][ieo];
	  for(int mu=0;mu<4;mu++)
	    su3_put_to_rnd(conf[par][ieo][mu],loc_rnd_gen[ilx]);
	}
      
      set_borders_invalid(conf[par]);
    }
}

//heatbath or overrelax algorithm for the quenched simulation case, Wilson action
void heatbath_or_overrelax_conf_Wilson_action(quad_su3 **eo_conf,theory_pars_type *theory_pars,pure_gauge_evol_pars_type *evol_pars,int heat_over)
{
  //loop on directions and on parity
  for(int mu=0;mu<4;mu++)
    for(int par=0;par<2;par++)
      {
        nissa_loc_volh_loop(ieo)
	  {
	    //find the lex index of the point and catch the random gen
	    int ilx=loclx_of_loceo[par][ieo];
	    rnd_gen *gen=&(loc_rnd_gen[ilx]);
	    
	    //compute the staples
	    su3 staples;
	    compute_point_staples_eo_conf_single_dir(staples,eo_conf,ilx,mu);
	    
	    //compute heatbath or overrelax link
	    su3 new_link;
	    if(heat_over==0) su3_find_heatbath(new_link,eo_conf[par][ieo][mu],staples,theory_pars->beta,evol_pars->nhb_hits,gen);
	    else             su3_find_overrelaxed(new_link,eo_conf[par][ieo][mu],staples,evol_pars->nov_hits);
	    
	    //change it
	    su3_copy(eo_conf[par][ieo][mu],new_link);
	  }
	
        //set the borders invalid: since we split conf in e/o, only now needed                                                                                              
        set_borders_invalid(eo_conf[par]);
      }
}

//cool the configuration
void cool_conf(quad_su3 **eo_conf,int over_flag,double over_exp)
{
  //loop on parity and directions
  for(int mu=0;mu<4;mu++)
    for(int par=0;par<2;par++)
      {
        nissa_loc_volh_loop(ieo)
        {
          //find the transformation
          su3 u;
          su3_find_cooled(u,eo_conf,par,ieo,mu);
            
          //overrelax if needed
          if(over_flag)
            {
              //find the transformation
              su3 temp1;
              unsafe_su3_prod_su3_dag(temp1,u,eo_conf[par][ieo][mu]);
                
              //exponentiate it and re-unitarize
              su3 temp2;
              su3_overrelax(temp2,temp1,over_exp);
                
              //find the transformed link
              unsafe_su3_prod_su3(u,temp2,eo_conf[par][ieo][mu]);
            }
            
          //change the link
          su3_copy(eo_conf[par][ieo][mu],u);
        }
        
        //now set the borders invalid: since we split conf in e/o, only now needed
        set_borders_invalid(eo_conf[par]);
      }
}

//used below, to check if the square is to be asked and in the case (if at 2 iter) note it
inline int check_add_square_staple(int *isquare_staples_to_ask,int &nsquare_staple_to_ask,int ivol,int dir,int verse,int iter)
{
  if(ivol>=loc_vol)
    {
      int *note=isquare_staples_to_ask+(nsquare_staple_to_ask++)*3;
      if(iter==1)
	{
	  note[0]=ivol;
	  note[1]=dir;
	  note[2]=verse;
	}
      return 1;
    }
  else return 0;
}

//heatbath or overrelax algorithm for the quenched simulation case, tlSym action
void heatbath_or_overrelax_conf_tlSym_action(quad_su3 **eo_conf,theory_pars_type *theory_pars,pure_gauge_evol_pars_type *evol_pars,int heat_over)
{
  //check that local volume is a multpiple of 4
  for(int mu=0;mu<4;mu++)
    if(loc_size[mu]%4!=0)
      crash("Direction %d has not local size multiple of 4",mu);
  
  //allocate the list of link to update
  int nlinks_to_update_exp=loc_vol/8;
  int *links_to_update=nissa_malloc("links_to_update",nlinks_to_update_exp,int);
  
  //allocate the 24 squared staples to be computed:
  //-6 horizontal, 2 mu-perp verses, 3 mu-perp directions at distance 1mu-perp
  //-6 horizontal, 2 mu-perp verses, 3 mu-perp directions at distance 2mu-perp
  //-12 vertical,  2 mu-perm verses, 3 mu-perp directions, 2 mu-paral verses at distance 1mu-paral
  int nsquare_staples=nlinks_to_update_exp*24;
  su3 *square_staples=nissa_malloc("square_staples",nsquare_staples,su3);
  
  //space where to take not of staples to ask
  int *isquare_staples_to_ask=NULL;
  
  //loop on directions
  for(int mu=0;mu<4;mu++)
    {
      //take the other three directions
      int nu=(mu+1)%4;
      int rho=(nu+1)%4;
      int sig=(rho+1)%4;
      
      //loop on parity of the mu coord direction
      for(int eo_mu=0;eo_mu<2;eo_mu++)
	//loop on quadruple parity
	for(int qpar=0;qpar<4;qpar++)
	  {
	    //print debug info
	    master_printf("mu=%d eo_mu=%d qpar=%d\n",mu,eo_mu,qpar);	   
	    
	    //find the links to be updated: the point having 
	    //x_mu%2==eo_mu and (x'_nu+x'_rho+x'_sigma)%4==qpar
	    //where x'=x%4
	    vector_reset(links_to_update);
	    int nlinks_to_update=0;
	    nissa_loc_vol_loop(ivol)
	      if(glb_coord_of_loclx[ivol][mu]%2==eo_mu && 
		 (glb_coord_of_loclx[ivol][nu]%4+
		  glb_coord_of_loclx[ivol][rho]%4+
		  glb_coord_of_loclx[ivol][sig]%4)%4==qpar)
		links_to_update[nlinks_to_update++]=ivol;
	    
	    //compute buffer size: needs to ask                   .. ..   
	    //for a staple, among the dotted ones              ..:__:__:.. 
	    //at first iteration count and allocate bu        :..|__I__|..: 
	    //at second take iteration take note                 |__|__|
	    int nsquare_to_ask_to_rank[8]={0}; //positives, than negative directions
	    for(int de=0;de<4;de++)
	      if(paral_dir[de])
		{
		  if(de!=mu) nsquare_to_ask_to_rank[de]=nsquare_to_ask_to_rank[4+de]=loc_vol/loc_size[de]/4/2;
		  else
		    {
		      //only if points are on inner mu+ border needs external nu/rho/sigma 
		      if(eo_mu==1) nsquare_to_ask_to_rank[mu]=2*3*loc_vol/loc_size[mu]/4; 
		      nsquare_to_ask_to_rank[4+mu]=0;
		    }
		}
	    
	    for(int k=0;k<8;k++)
	      master_printf("%d ",nsquare_to_ask_to_rank[k]);
	    master_printf("\n");
	    
	    int nsquare_staple_to_ask=0;
	    int iter=0;
	    //recheck in another way
	    for(int ilink=0;ilink<nlinks_to_update;ilink++)
	      {
		int ivol=links_to_update[ilink];

		for(int nu=0;nu<4;nu++)
		  if(nu!=mu)
		    {
		      //find Left, TopLeft, Top and Right points
		      int L=loclx_neighdw[ivol][nu];
		      int T=loclx_neighup[ivol][mu];
		      int TL=loclx_neighdw[T][nu];
		      int R=loclx_neighup[ivol][nu];
		      
		      //check left staple: it has to be asked if L is on the border
		      if(check_add_square_staple(isquare_staples_to_ask,nsquare_staple_to_ask,L,mu,-1,iter)) nsquare_to_ask_to_rank[4+nu]--;
		      
		      //check top left staple: it has to be asked if TL is on border or edge
		      if(check_add_square_staple(isquare_staples_to_ask,nsquare_staple_to_ask,TL,nu,+1,iter)) nsquare_to_ask_to_rank[mu]--;
		      
		      //check top right staple: it has to be asked it T is on border
		      if(check_add_square_staple(isquare_staples_to_ask,nsquare_staple_to_ask,T,nu,+1,iter)) nsquare_to_ask_to_rank[mu]--;
		      
		      //check right staple: it has to be asked if R is on boder
		      if(check_add_square_staple(isquare_staples_to_ask,nsquare_staple_to_ask,R,mu,+1,iter)) nsquare_to_ask_to_rank[nu]--;
		    }
	      }
	    
	    for(int k=0;k<8;k++)
	      master_printf("%d ",nsquare_to_ask_to_rank[k]);
	    master_printf("\n");
	    
	    //allocate 3 int per square staple ot ask
	    //if(iter==0) isquare_staples_to_ask=nissa_malloc("istta",nsquare_staple_to_ask*3,int);
	  }

      
      //initialize space for square staples to zero
      //vector_reset(square_staples);
      
      //allocate the in and out buffers: they have the same size
      
      //compute the staples to be sent and store them in the buffer
      //location of staples to be sent is just the same of the staples to 
      //ask, since the lattice is translation invariant w.r.t nodes
      
      //start sending the staples
      
      //compute local staples
      
      //wait to receive staples
      
      //put incoming staples in place
      
      //compute rectangular staples
    }
  
  nissa_free(links_to_update);
}

//heatbath or overrelax algorithm for the quenched simulation case
void heatbath_or_overrelax_conf(quad_su3 **eo_conf,theory_pars_type *theory_pars,pure_gauge_evol_pars_type *evol_pars,int heat_over)
{
  switch(theory_pars->gac_type)
    {
    case Wilson_action:
      heatbath_or_overrelax_conf_Wilson_action(eo_conf,theory_pars,evol_pars,heat_over);
      break;
   case tlSym_action:
      heatbath_or_overrelax_conf_tlSym_action(eo_conf,theory_pars,evol_pars,heat_over);
      break;
    default:
      crash("Unknown action");
    }
}

//heatbath algorithm for the quenched simulation case
void heatbath_conf(quad_su3 **eo_conf,theory_pars_type *theory_pars,pure_gauge_evol_pars_type *evol_pars)
{heatbath_or_overrelax_conf(eo_conf,theory_pars,evol_pars,0);}

//overrelax algorithm for the quenched simulation case
void overrelax_conf(quad_su3 **eo_conf,theory_pars_type *theory_pars,pure_gauge_evol_pars_type *evol_pars)
{heatbath_or_overrelax_conf(eo_conf,theory_pars,evol_pars,1);}
