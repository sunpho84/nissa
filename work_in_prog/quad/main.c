#include "nissa.h"

void compute_residue_128(spincolor *source,quad_su3 *conf,double kappa,double mass,spincolor *temp,spincolor *temp_sol,spincolor *tot_source)
{

  apply_tmQ2(source,conf,kappa,mass,temp,temp_sol);
  
  float_128 res_128={0,0};
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  {
	    source[ivol][id][ic][ri]=tot_source[ivol][id][ic][ri]-source[ivol][id][ic][ri];
	    float_128 temp={source[ivol][id][ic][ri]*source[ivol][id][ic][ri],0};
	    float_128_summassign(res_128,temp);
	  }
  
  master_printf("RESIDUE: %lg\n",res_128[0]);
}


//add the solution to the total one
void solution_add_128(spincolor *tot_sol,spincolor *sol)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  tot_sol[ivol][id][ic][ri]+=sol[ivol][id][ic][ri];
}

int main()
{
  init_nissa();
  
  //init the grid
  init_grid(8,4);
  
  //read conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+loc_bord+loc_edge,quad_su3);
  read_ildg_gauge_conf(conf,"../../test/data/L4T8conf");
  
  //prepare the total source
  spincolor *tot_source=nissa_malloc("tot_source",loc_vol+loc_bord,spincolor);
  memset(tot_source,0,loc_vol*sizeof(spincolor));
  if(rank==0) tot_source[0][0][0][0]=1;
  set_borders_invalid(tot_source);
  
  //allocate the total solution
  spincolor *tot_sol=nissa_malloc("tot_sol",loc_vol+loc_bord,spincolor);
  memset(tot_sol,0,loc_vol*sizeof(spincolor));
  
  //set kappa, mass and residue
  double mass=0.50;
  double kappa=0.177;
  double residue=1.e-10;
  
  //allocate the source, the solution and a temporary vec
  spincolor *source=nissa_malloc("source",loc_vol+loc_bord,spincolor);
  spincolor *sol=nissa_malloc("sol",loc_vol+loc_bord,spincolor);
  spincolor *temp=nissa_malloc("temp",loc_vol+loc_bord,spincolor);
 
  for(int iter=0;iter<5;iter++)
    {
      //compute the source: source = tot_source - Q2 * tot_sol
      compute_residue_128(source,conf,kappa,mass,temp,tot_sol,tot_source);
      
      //compute partial sol
      inv_tmQ2_cg(sol,NULL,conf,kappa,mass,100000,1,residue,source);

      //add the solution to the total one
      solution_add_128(tot_sol,sol);
    }
  
  nissa_free(source);
  nissa_free(sol);
  nissa_free(tot_source);
  nissa_free(tot_sol);
  nissa_free(conf);
  
  return 0;
}
