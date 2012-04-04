#include "nissa.h"
#include "dirac_operator_tmQ_128_portable.c"

//compute the residue in 128 bit
double compute_tmQ2_residue_128(spincolor *residue,quad_su3 *conf,double kappa,double mass,spincolor_128 *temp,spincolor_128 *tot_sol,spincolor *tot_source)
{
  //compute D*sol in quadruple precision
  spincolor_128 *residue_128=nissa_malloc("residue_128",loc_vol+loc_bord,spincolor_128);
  apply_tmQ2_128(residue_128,conf,kappa,mass,temp,tot_sol);
  
  //compute the residue
  float_128 loc_res_summ_128={0,0};
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  {
	    float_128_subt_from_64(residue_128[ivol][id][ic][ri],tot_source[ivol][id][ic][ri],residue_128[ivol][id][ic][ri]);
	    float_128_summ_the_prod(loc_res_summ_128,residue_128[ivol][id][ic][ri],residue_128[ivol][id][ic][ri]);
	    residue[ivol][id][ic][ri]=double_from_float_128(residue_128[ivol][id][ic][ri]);
	  }
  //perform global reduction of the residue
  float_128 res_summ_128;
  MPI_Allreduce(loc_res_summ_128,res_summ_128,1,MPI_FLOAT_128,MPI_FLOAT_128_SUM,MPI_COMM_WORLD);
  
  //set borders invalid and free
  set_borders_invalid(residue);  
  nissa_free(residue_128);
  
  return res_summ_128[0];
}

//add the solution to the total one
void solution_add_128(spincolor_128 *tot_sol,spincolor *sol)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  float_128_summassign_64(tot_sol[ivol][id][ic][ri],sol[ivol][id][ic][ri]);
  
  set_borders_invalid(tot_sol);
}

void inv_tmQ2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double mass,int niter,double external_solver_residue,spincolor *source)
{
  //allocate temporary vectorthe solution in 128 bit
  spincolor *residue_source=nissa_malloc("residue_source",loc_vol+loc_bord,spincolor);
  spincolor_128 *temp_128=nissa_malloc("temp",loc_vol+loc_bord,spincolor_128); 
  spincolor_128 *sol_128=nissa_malloc("sol_128",loc_vol+loc_bord,spincolor_128);
  memset(sol_128,0,loc_vol*sizeof(spincolor_128));
  set_borders_invalid(sol_128);
  set_borders_invalid(temp_128);
  
  //internal solver stopping condition
  double inner_solver_residue=min_double(1.e-20,external_solver_residue);

  //invert until residue is lower than required
  double current_residue;
  do
    {
      //compute the residue source: residue_source = source - Q2 * sol_128
      current_residue=compute_tmQ2_residue_128(residue_source,conf,kappa,mass,temp_128,sol_128,source);
      master_printf("\nExternal loop residue: %lg\n\n",current_residue);

      //check if we residue not reached
      if(current_residue>=external_solver_residue)
	{
	  //compute partial sol
	  inv_tmQ2_cg(sol,NULL,conf,kappa,mass,100000,1,inner_solver_residue,residue_source);
	  
	  //add the solution to the total one
	  solution_add_128(sol_128,sol);
	}
    }
  while(current_residue>=external_solver_residue);
  
  //copy the solution
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  sol[ivol][id][ic][ri]=double_from_float_128(sol_128[ivol][id][ic][ri]);
  
  nissa_free(temp_128);
  nissa_free(sol_128);
  nissa_free(residue_source);
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();

  //initialize the program
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //read the input
  open_input(arg[1]);

  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Kappa
  double kappa;
  read_str_double("Kappa",&kappa);
  //Mass
  double mass;
  read_str_double("Mass",&mass);
  //Residue
  double residue;
  read_str_double("Residue",&residue);
  //Path of conf
  char gauge_conf_path[1024];
  read_str_str("GaugeConfPath",gauge_conf_path,1024);
  
  close_input();
  
  //read conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+loc_bord+loc_edge,quad_su3);
  read_ildg_gauge_conf(conf,gauge_conf_path);
  
  //prepare the total source
  spincolor *source=nissa_malloc("source",loc_vol+loc_bord,spincolor);
  memset(source,0,loc_vol*sizeof(spincolor));
  if(rank==0) source[0][0][0][0]=1;
  set_borders_invalid(source);
  
  //allocate solution
  spincolor *sol=nissa_malloc("sol",loc_vol+loc_bord,spincolor);
  
  inv_tmQ2_cg_128(sol,NULL,conf,kappa,mass,10000,residue,source);
  
  nissa_free(sol);
  nissa_free(source);
  nissa_free(conf);
  
  close_nissa();
  
  return 0;
}
