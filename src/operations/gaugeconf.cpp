#include <math.h>
#include <string.h>

#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../base/debug.h"
#include "../base/communicate.h"
#include "../base/random.h"
#include "../base/routines.h"
#include "../new_types/complex.h"
#include "../new_types/su3.h"
#include "../operations/remap_vector.h"
#include "../operations/su3_paths.h"
#include "../geometry/geometry_mix.h"

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
  
  if(putonbords) set_borders_invalid(conf);
  if(putonedges) set_edges_invalid(conf);
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
      master_printf("Necesary to add boundary condition: %f %f %f %f\n",diff_theta[0],diff_theta[1],diff_theta[2],diff_theta[3]);
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

//return a cooled copy of the passed link
void find_cooled_link(su3 u,quad_su3 **eo_conf,int par,int ieo,int mu)
{
  //compute the staple
  su3 staple;
  compute_point_staples_eo_conf_single_dir(staple,eo_conf,loclx_of_loceo[par][ieo],mu);
  
  //find the link that maximize the plaquette
  su3_unitarize_orthonormalizing(u,staple);
}

//cool the configuration
void cool_conf(quad_su3 **eo_conf,int over_flag,double over_exp)
{
  //loop first on parity and then on directions
  for(int mu=0;mu<4;mu++)
    for(int par=0;par<2;par++)
      {
	nissa_loc_volh_loop(ieo)
	  {
	    //find the transformation
	    su3 u;
	    find_cooled_link(u,eo_conf,par,ieo,mu);
	    
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
