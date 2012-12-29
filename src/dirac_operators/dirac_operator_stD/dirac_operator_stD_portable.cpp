#include "../../base/routines.h"
#include <omp.h>

double app_time=0;
int napp=0;

void apply_st2Doe(color *out,quad_su3 **conf,color *in)
{
  communicate_eo_quad_su3_borders(conf);
  communicate_ev_color_borders(in);
  
#pragma omp parallel for
  nissa_loc_volh_loop(io)
    {
      //printf("nthreads: %d\n",omp_get_num_threads());
      //neighbours search
      int evup0=loceo_neighup[ODD][io][0];
      int evdw0=loceo_neighdw[ODD][io][0];
      
      //derivative in the time direction - without self-summ
      unsafe_su3_prod_color(      out[io],conf[ODD][io   ][0],in[evup0]);
      su3_dag_subt_the_prod_color(out[io],conf[EVN][evdw0][0],in[evdw0]);
      
      //derivatives in the spatial direction - with self summ
      for(int mu=1;mu<4;mu++)
	{
	  int evup=loceo_neighup[ODD][io][mu];
	  int evdw=loceo_neighdw[ODD][io][mu];
	  
	  su3_summ_the_prod_color(    out[io],conf[ODD][io  ][mu],in[evup]);
	  su3_dag_subt_the_prod_color(out[io],conf[EVN][evdw][mu],in[evdw]);
	}
    }
  
  set_borders_invalid(out);
}

//put the 0.5 factor
void apply_stDoe(color *out,quad_su3 **conf,color *in)
{
  apply_st2Doe(out,conf,in);

  nissa_loc_volh_loop(io)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	out[io][ic][ri]*=0.5;
}

//multply also for an additional 1/2
void apply_stDeo_half(color *out,quad_su3 **conf,color *in)
{
  communicate_eo_quad_su3_borders(conf);
  communicate_od_color_borders(in);
  
#pragma omp parallel for
  nissa_loc_volh_loop(ie)
    {
      int odup0=loceo_neighup[EVN][ie][0];
      int oddw0=loceo_neighdw[EVN][ie][0];
      
      unsafe_su3_dag_prod_color(out[ie],conf[ODD][oddw0][0],in[oddw0]);
      su3_subt_the_prod_color(  out[ie],conf[EVN][ie   ][0],in[odup0]);
      
      for(int mu=1;mu<4;mu++)
	{
	  int odup=loceo_neighup[EVN][ie][mu];
	  int oddw=loceo_neighdw[EVN][ie][mu];
	  
	  su3_dag_summ_the_prod_color(out[ie],conf[ODD][oddw][mu],in[oddw]);
	  su3_subt_the_prod_color(    out[ie],conf[EVN][ie  ][mu],in[odup]);
	}
      
      color_prod_double(out[ie],out[ie],0.25);      
    }
  
  set_borders_invalid(out);
}

void aapply_stD2ee(color *out,quad_su3 **conf,color *temp,double mass,color *in)
{
  app_time-=take_time();
  
  communicate_eo_quad_su3_borders(conf);
  communicate_ev_color_borders(in);
  
  //check arguments
  if(out==in)   crash("out==in!");
  if(out==temp) crash("out==temp!");
  if(temp==in)  crash("temp==in!");
  
  double mass2=mass*mass;
  
  //perform the off diagonal multiplication
  communicate_eo_quad_su3_borders(conf);
  communicate_ev_color_borders(in);

#pragma omp parallel
  {
#pragma omp for
    nissa_loc_volh_loop(io)
    {
      //printf("nthreads: %d\n",omp_get_num_threads());
      //neighbours search
      int evup0=loceo_neighup[ODD][io][0];
      int evdw0=loceo_neighdw[ODD][io][0];
      
      //derivative in the time direction - without self-summ
      unsafe_su3_prod_color(      temp[io],conf[ODD][io   ][0],in[evup0]);
      su3_dag_subt_the_prod_color(temp[io],conf[EVN][evdw0][0],in[evdw0]);
      
      //derivatives in the spatial direction - with self summ
      for(int mu=1;mu<4;mu++)
	{
	  int evup=loceo_neighup[ODD][io][mu];
	  int evdw=loceo_neighdw[ODD][io][mu];
	  
	  su3_summ_the_prod_color(    temp[io],conf[ODD][io  ][mu],in[evup]);
	  su3_dag_subt_the_prod_color(temp[io],conf[EVN][evdw][mu],in[evdw]);
	}
    }
  
    #pragma omp single
    {
      set_borders_invalid(temp);
      communicate_od_color_borders(temp);
    }
    
#pragma omp for
  nissa_loc_volh_loop(ie)
    {
      int odup0=loceo_neighup[EVN][ie][0];
      int oddw0=loceo_neighdw[EVN][ie][0];
      
      unsafe_su3_dag_prod_color(out[ie],conf[ODD][oddw0][0],temp[oddw0]);
      su3_subt_the_prod_color(  out[ie],conf[EVN][ie   ][0],temp[odup0]);
      
      for(int mu=1;mu<4;mu++)
	{
	  int odup=loceo_neighup[EVN][ie][mu];
	  int oddw=loceo_neighdw[EVN][ie][mu];
	  
	  su3_dag_summ_the_prod_color(out[ie],conf[ODD][oddw][mu],temp[oddw]);
	  su3_subt_the_prod_color(    out[ie],conf[EVN][ie  ][mu],temp[odup]);
	}
      
      color_prod_double(out[ie],out[ie],0.25);      
    }
  
  //summ the mass
#pragma omp for
  nissa_loc_volh_loop(ivol)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	out[ivol][ic][ri]+=mass2*in[ivol][ic][ri];
  }
  
  set_borders_invalid(out);

  app_time+=take_time();
  napp++;
}

void apply_stD2ee(color *out,quad_su3 **conf,color *temp,double mass,color *in)
{
  app_time-=take_time();
  
  communicate_eo_quad_su3_borders(conf);
  communicate_ev_color_borders(in);
  
  //check arguments
  if(out==in)   crash("out==in!");
  if(out==temp) crash("out==temp!");
  if(temp==in)  crash("temp==in!");
  
  double mass2=mass*mass;
  
  //perform the off diagonal multiplication
  apply_st2Doe(temp,conf,in);
  apply_stDeo_half(out,conf,temp);
  
  //summ the mass
  nissa_loc_volh_loop(ivol)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	out[ivol][ic][ri]+=mass2*in[ivol][ic][ri];
  
  set_borders_invalid(out);

  app_time+=take_time();
  napp++;
}
