#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "../../communicate/communicate.h"
#include "../../base/debug.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../geometry/geometry_lx.h"
#include "../../geometry/geometry_mix.h"
#include "../../new_types/complex.h"
#include "../../new_types/dirac.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/spin.h"
#include "../../new_types/su3.h"
#include "../shift.h"
#include "../../routines/mpi_routines.h"

//compute the polyakov loop
void average_polyakov_loop(complex tra,quad_su3 *conf,int mu)
{
  su3 *u=nissa_malloc("u",loc_vol+bord_vol,su3);
  
  communicate_lx_quad_su3_borders(conf);
  
  //reset the link product
  nissa_loc_vol_loop(ivol)
    su3_put_to_id(u[ivol]);

  //move along +mu
  for(int i=0;i<glb_size[mu];i++)
    {
      //take the product
      nissa_loc_vol_loop(ivol)
        safe_su3_prod_su3(u[ivol],u[ivol],conf[ivol][mu]);
      set_borders_invalid(u);
      
      su3_vec_single_shift(u,mu,+1);
    }
  
  //compute the trace; since we reduce over all the volume there are glb_size[mu] replica
  complex loc_tra={0,0};
  nissa_loc_vol_loop(ivol)
    su3_summ_the_trace(loc_tra,u[ivol]);
  glb_reduce_complex(tra,loc_tra);
  complex_prodassign_double(tra,1.0/glb_vol/3.0);
  
  nissa_free(u);
}

//definition in case of eos conf
void average_polyakov_loop_of_eos_conf(complex tra,quad_su3 **eo_conf,int mu)
{
  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol,quad_su3);
  paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
  
  average_polyakov_loop(tra,lx_conf,mu);
  
  nissa_free(lx_conf);
}

//Compute the Pline in a certain direction mu, starting from xmu_start.
//Between xmu_start and (xmu_start+glb_size[mu]/2) it contains forward line
//Between xmu_start and (xmu_start-glb_size[mu]/2) it contains backward line
//At (xmu_start-glb_size[mu]/2) is done twice and left with the forward
//Actually since what is needed for Wprop is the revert, it is dag
void compute_Pline_dag_internal(su3 *pline,quad_su3 *conf,int mu,int xmu_start)
{
  communicate_lx_quad_su3_borders(conf);
  
  //Loop simultaneously forward and backward
  for(int xmu_shift=1;xmu_shift<=glb_size[mu]/2;xmu_shift++)
    {
      communicate_lx_su3_borders(pline);
      
      nissa_loc_vol_loop(x)
      {
	int x_back=loclx_neighdw[x][mu];
	int x_forw=loclx_neighup[x][mu];
          
	//consider +xmu_shift point
	if(glb_coord_of_loclx[x][mu]==(xmu_start+xmu_shift)%glb_size[mu]) unsafe_su3_dag_prod_su3(pline[x],conf[x_back][mu],pline[x_back]);
	//consider -xmu_shift point
	if((glb_coord_of_loclx[x][mu]+xmu_shift)%glb_size[mu]==xmu_start) unsafe_su3_prod_su3(pline[x],conf[x][mu],pline[x_forw]);
      }
      
      set_borders_invalid(pline);
    }
}

///compute the pline daggered stemming from the all the points at xmu=xmu_start along dir mu
void compute_Pline_dag_wall(su3 *pline,quad_su3 *conf,int mu,int xmu_start)
{
  //reset the link product, putting id at x such that xmu==xmu_start
  vector_reset(pline);
  nissa_loc_vol_loop(x) if(glb_coord_of_loclx[x][mu]==xmu_start) su3_put_to_id(pline[x]);

  compute_Pline_dag_internal(pline,conf,mu,xmu_start);
}

//compute the Pline daggered_stemming from x_start along dir mu
void compute_Pline_dag_point(su3 *pline,quad_su3 *conf,int mu,coords glb_x_start)
{
  //get the rank and loc site x
  int loc_x_start,rank_hosting_x;
  get_loclx_and_rank_of_coord(&loc_x_start,&rank_hosting_x,glb_x_start);
  
  //reset the link product, putting id at x_start
  vector_reset(pline);
  if(rank==rank_hosting_x) su3_put_to_id(pline[loc_x_start]);
  
  compute_Pline_dag_internal(pline,conf,mu,glb_x_start[mu]);
}

//Compute the stochastic Pline, using a color as source
void compute_stoch_Pline_dag(color *pline,quad_su3 *conf,int mu,int xmu_start,color *source)
{
  communicate_lx_quad_su3_borders(conf);
  
  //Reset the link product, putting id at xmu_start
  vector_reset(pline);
  nissa_loc_vol_loop(ivol)
    if(glb_coord_of_loclx[ivol][mu]==xmu_start)
      color_copy(pline[ivol],source[ivol]);
  
  //Loop simultaneously forward and backward
  for(int xmu_shift=1;xmu_shift<=glb_size[mu]/2;xmu_shift++)
    {
      communicate_lx_color_borders(pline);
      
      nissa_loc_vol_loop(x)
      {
	int x_back=loclx_neighdw[x][mu];
	int x_forw=loclx_neighup[x][mu];
          
	//consider +xmu_shift point
	if(glb_coord_of_loclx[x][mu]==(xmu_start+xmu_shift)%glb_size[mu]) unsafe_su3_dag_prod_color(pline[x],conf[x_back][mu],pline[x_back]);
	//consider -xmu_shift point
	if((glb_coord_of_loclx[x][mu]+xmu_shift)%glb_size[mu]==xmu_start) unsafe_su3_prod_color(pline[x],conf[x][mu],pline[x_forw]);
      }
      
      set_borders_invalid(pline);
    }
}

//compute the static propagator, putting the 1+gamma_mu in place
void compute_Wstat_prop_finalize(su3spinspin *prop,quad_su3 *conf,int mu,int xmu_start,su3 *pline)
{
  //reset the output
  vector_reset(prop);
  
  //take the gamma
  dirac_matr *gamma_mu=base_gamma+nissa_map_mu[mu];
  
  nissa_loc_vol_loop(x)
  {
    int xmu=glb_coord_of_loclx[x][mu];
    int dist=fabs(xmu-xmu_start);
    int ord=(xmu>=xmu_start);
      
    for(int ic1=0;ic1<3;ic1++)
      for(int ic2=0;ic2<3;ic2++)
	{
	  spinspin_dirac_summ_the_prod_complex(prop[x][ic1][ic2],base_gamma+0,pline[x][ic1][ic2]);
            
	  //sign of 1+-gamma_mu
	  if((ord==1 && dist<=glb_size[mu]/2)||(ord==0 && dist>=glb_size[mu]/2)) spinspin_dirac_summ_the_prod_complex(prop[x][ic1][ic2],gamma_mu,pline[x][ic1][ic2]); //forward
	  else                                                                   spinspin_dirac_subt_the_prod_complex(prop[x][ic1][ic2],gamma_mu,pline[x][ic1][ic2]); //backward
            
	  spinspin_prodassign_double(prop[x][ic1][ic2],0.5);
	}
  }
  
  set_borders_invalid(prop);
}

//compute the static propagator
void compute_Wstat_prop_wall(su3spinspin *prop,quad_su3 *conf,int mu,int xmu_start)
{
  su3 *pline=nissa_malloc("pline",loc_vol+bord_vol,su3);
  
  //version with pline stemming from a wall
  compute_Pline_dag_wall(pline,conf,mu,xmu_start);
  
  compute_Wstat_prop_finalize(prop,conf,mu,xmu_start,pline);

  nissa_free(pline);
}

//version with pline by a point
void compute_Wstat_prop_point(su3spinspin *prop,quad_su3 *conf,int mu,coords x_start)
{
  su3 *pline=nissa_malloc("pline",loc_vol+bord_vol,su3);

  //compute pline stemming from a point
  compute_Pline_dag_point(pline,conf,mu,x_start);
  
  compute_Wstat_prop_finalize(prop,conf,mu,x_start[mu],pline);

  nissa_free(pline);
}

//compute the stochastic static propagator, putting the 1+gamma_mu in place
void compute_Wstat_stoch_prop(colorspinspin *prop,quad_su3 *conf,int mu,int xmu_start,color *source)
{
  color *pline=nissa_malloc("pline",loc_vol+bord_vol,color);

  //compute stocasthic pline
  compute_stoch_Pline_dag(pline,conf,mu,xmu_start,source);

  //reset the output
  vector_reset(prop);

  //take the gamma
  dirac_matr *gamma_mu=base_gamma+nissa_map_mu[mu];

  nissa_loc_vol_loop(x)
  {
    int xmu=glb_coord_of_loclx[x][mu];
    int dist=fabs(xmu-xmu_start);
    int ord=(xmu>=xmu_start);

    for(int ic=0;ic<3;ic++)
      {
        spinspin_dirac_summ_the_prod_complex(prop[x][ic],base_gamma+0,pline[x][ic]);

        if((ord==1 && dist<=glb_size[mu]/2)||(ord==0 && dist>=glb_size[mu]/2)) spinspin_dirac_summ_the_prod_complex(prop[x][ic],gamma_mu,pline[x][ic]); //forward
        else                                                                   spinspin_dirac_subt_the_prod_complex(prop[x][ic],gamma_mu,pline[x][ic]); //backward
        
        spinspin_prodassign_double(prop[x][ic],0.5);
      }
  }

  set_borders_invalid(prop);

  nissa_free(pline);
}
