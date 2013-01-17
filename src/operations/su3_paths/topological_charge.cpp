#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "../gaugeconf.h"

#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"
#include "../../base/communicate.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/complex.h"
#include "../../new_types/spin.h"
#include "../../new_types/su3.h"
#include "../../geometry/geometry_mix.h"

//This will calculate 2*a^2*ig*P_{mu,nu}
/*
  ^                   C--<-- B --<--Y 
  |                   |  2  | |  1  | 
  n                   |     | |     | 
  u                   D-->--\X/-->--A 
  |                   D--<--/X\--<--A 
  -----mu---->        |  3  | |  4  | 
  		      |     | |     | 
		      E-->-- F -->--G 
in order to have the anti-simmetric part, use
the routine "Pmunu_term"
*/
void four_leaves(as2t_su3 *Pmunu,quad_su3 *conf)
{
  communicate_lx_quad_su3_edges(conf);
  
  int A,B,C,D,E,F,G;
  int munu;

  su3 temp1,temp2,leaves_summ;

  nissa_loc_vol_loop(X)
    {
      as2t_su3_put_to_zero(Pmunu[X]);

      munu=0;
      for(int mu=0;mu<4;mu++)
	{
	  A=loclx_neighup[X][mu];
	  D=loclx_neighdw[X][mu];
	  
	  for(int nu=mu+1;nu<4;nu++)
	    {
	      B=loclx_neighup[X][nu];
	      F=loclx_neighdw[X][nu];
	      
	      C=loclx_neighup[D][nu];
	      E=loclx_neighdw[D][nu];
	      
	      G=loclx_neighdw[A][nu];
	      
	      //Put to 0 the summ of the leaves
	      su3_put_to_zero(leaves_summ);
	      
	      //Leaf 1
	      unsafe_su3_prod_su3(temp1,conf[X][mu],conf[A][nu]);         //    B--<--Y 
	      unsafe_su3_prod_su3_dag(temp2,temp1,conf[B][mu]);           //    |  1  | 
	      unsafe_su3_prod_su3_dag(temp1,temp2,conf[X][nu]);           //    |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);	                  //    X-->--A 
	      
	      //Leaf 2
	      unsafe_su3_prod_su3_dag(temp1,conf[X][nu],conf[C][mu]);     //    C--<--B
	      unsafe_su3_prod_su3_dag(temp2,temp1,conf[D][nu]);           //    |  2  | 
	      unsafe_su3_prod_su3(temp1,temp2,conf[D][mu]);		  //    |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);		          //    D-->--X
	      
	      //Leaf 3
	      unsafe_su3_dag_prod_su3_dag(temp1,conf[D][mu],conf[E][nu]);  //   D--<--X
	      unsafe_su3_prod_su3(temp2,temp1,conf[E][mu]);		   //   |  3  | 
	      unsafe_su3_prod_su3(temp1,temp2,conf[F][nu]);		   //   |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);		           //   E-->--F
	      
	      //Leaf 4
	      unsafe_su3_dag_prod_su3(temp1,conf[F][nu],conf[F][mu]);       //  X--<--A 
	      unsafe_su3_prod_su3(temp2,temp1,conf[G][nu]);                 //  |  4  | 
	      unsafe_su3_prod_su3_dag(temp1,temp2,conf[X][mu]);             //  |     |  
	      su3_summ(leaves_summ,leaves_summ,temp1);                      //  F-->--G 
	      
	      su3_copy(Pmunu[X][munu],leaves_summ);
	      
	      munu++;
	    }
	}
    }
  
  set_borders_invalid(Pmunu);
}

//takes the anti-simmetric part of the four-leaves
void Pmunu_term(as2t_su3 *Pmunu,quad_su3 *conf)
{
  four_leaves(Pmunu,conf);
  
  //calculate U-U^dagger
  nissa_loc_vol_loop(X)
    for(int munu=0;munu<6;munu++)
      {
	//bufferized antisimmetrization
	su3 leaves_summ;
	memcpy(leaves_summ,Pmunu[X][munu],sizeof(su3));
	
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    {
	      Pmunu[X][munu][ic1][ic2][0]=(leaves_summ[ic1][ic2][0]-leaves_summ[ic2][ic1][0])/4;
	      Pmunu[X][munu][ic1][ic2][1]=(leaves_summ[ic1][ic2][1]+leaves_summ[ic2][ic1][1])/4;
	    }
      }
}

//apply the chromo operator to the passed spinor site by site (not yet fully optimized)
void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,as2t_su3 Pmunu,spincolor in)
{
  color temp_d1;
  
  for(int d1=0;d1<4;d1++)
    {
      color_put_to_zero(out[d1]);
      for(int imunu=0;imunu<6;imunu++)
	{
	  unsafe_su3_prod_color(temp_d1,Pmunu[imunu],in[smunu_pos[d1][imunu]]);
	  for(int c=0;c<3;c++) complex_summ_the_prod(out[d1][c],smunu_entr[d1][imunu],temp_d1[c]);
	}
    }
}

//apply the chromo operator to the passed spinor to the whole volume
void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,as2t_su3 *Pmunu,spincolor *in)
{
  nissa_loc_vol_loop(ivol) unsafe_apply_point_chromo_operator_to_spincolor(out[ivol],Pmunu[ivol],in[ivol]);
}

//apply the chromo operator to the passed colorspinspin
//normalization as in ape next
void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,as2t_su3 *Pmunu,colorspinspin *in)
{
  spincolor temp1,temp2;
  
  nissa_loc_vol_loop(ivol)
    {
      //Loop over the four source dirac indexes
      for(int id_source=0;id_source<4;id_source++) //dirac index of source
	{
	  //Switch the color_spinspin into the spincolor.
	  get_spincolor_from_colorspinspin(temp1,in[ivol],id_source);
	  
	  unsafe_apply_point_chromo_operator_to_spincolor(temp2,Pmunu[ivol],temp1);
	  
	  //Switch back the spincolor into the colorspinspin
	  put_spincolor_into_colorspinspin(out[ivol],temp2,id_source);
	}
    }
  
  //invalidate borders
  set_borders_invalid(out);
}

//apply the chromo operator to the passed su3spinspin
//normalization as in ape next
void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin *out,as2t_su3 *Pmunu,su3spinspin *in)
{
  spincolor temp1,temp2;
  
  nissa_loc_vol_loop(ivol)
    {
      //Loop over the four source dirac indexes
      for(int id_source=0;id_source<4;id_source++) //dirac index of source
	for(int ic_source=0;ic_source<3;ic_source++) //color index of source
	  {
	    //Switch the su3spinspin into the spincolor.
	    get_spincolor_from_su3spinspin(temp1,in[ivol],id_source,ic_source);
	    
	    unsafe_apply_point_chromo_operator_to_spincolor(temp2,Pmunu[ivol],temp1);
	    
	    //Switch back the spincolor into the colorspinspin
	    put_spincolor_into_su3spinspin(out[ivol],temp2,id_source,ic_source);
	  }
    }

  //invalidate borders
  set_borders_invalid(out);
}

//measure the topological charge site by site
void local_topological_charge(double *charge,quad_su3 *conf)
{
  double norm_fact=1/(128*M_PI*M_PI);
  
  as2t_su3 *leaves=nissa_malloc("leaves",loc_vol,as2t_su3);
  
  vector_reset(charge);
  
  //compute the clover-shape paths
  four_leaves(leaves,conf);
  
  //list the three combinations of plans
  int plan_id[3][2]={{0,5},{1,4},{2,3}};
  int sign[3]={1,-1,1};
  
  //loop on the three different combinations of plans
  for(int iperm=0;iperm<3;iperm++)
    {
      //take the index of the two plans
      int ip0=plan_id[iperm][0];
      int ip1=plan_id[iperm][1];
      
      nissa_loc_vol_loop(ivol)
        {
	  //products
	  su3 clock,aclock;
	  unsafe_su3_prod_su3_dag(clock,leaves[ivol][ip0],leaves[ivol][ip1]);
	  unsafe_su3_prod_su3(aclock,leaves[ivol][ip0],leaves[ivol][ip1]);
	  
	  //take the trace
	  complex tclock,taclock;
	  su3_trace(tclock,clock);
	  su3_trace(taclock,aclock);
	  
	  //takes the combination with appropriate sign
	  charge[ivol]+=sign[iperm]*(tclock[RE]-taclock[RE])*norm_fact;
	}
    }
  
  set_borders_invalid(charge);
  
  nissa_free(leaves);
}

//average the topological charge
double average_topological_charge(quad_su3 *conf)
{
  double *charge=nissa_malloc("charge",loc_vol,double);
  
  //compute local charge
  local_topological_charge(charge,conf);
  
  //average over local volume
  double ave_charge=0;
  nissa_loc_vol_loop(ivol)
    ave_charge+=charge[ivol];
  
  nissa_free(charge);
  
  //return the reduction over all nodes
  return glb_reduce_double(ave_charge);
}

//wrapper for eos case
double average_topological_charge(quad_su3 **eo_conf)
{
  //convert to lex
  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
  
  double charge=average_topological_charge(lx_conf);
  
  nissa_free(lx_conf);
  
  return charge;
}

//measure the topologycal charge
void measure_topology(top_meas_pars_type &pars,quad_su3 **uncooled_conf,int iconf)
{
  FILE *file=open_file(pars.path,(iconf==0)?"w":"a");
  
  //allocate a temorary conf to be cooled
  quad_su3 *cooled_conf[2];
  for(int par=0;par<2;par++)
    {
      cooled_conf[par]=nissa_malloc("cooled_conf",loc_volh+bord_volh+edge_volh,quad_su3);
      vector_copy(cooled_conf[par],uncooled_conf[par]);
    }
  
  //print curent measure and cool
  for(int istep=0;istep<=(pars.cool_nsteps/pars.meas_each)*pars.meas_each;istep++)
    {
      if(istep%pars.meas_each==0) master_fprintf(file,"%d %d %16.16lg\n",iconf,istep,average_topological_charge(cooled_conf));
      if(istep!=pars.cool_nsteps) cool_conf(cooled_conf,pars.cool_overrelax_flag,pars.cool_overrelax_exp);
    }
  
  //discard cooled conf
  for(int par=0;par<2;par++) nissa_free(cooled_conf[par]);
  
  if(rank==0) fclose(file);
}
