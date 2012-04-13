#include "nissa.h"

void inv_tmQ2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double mass,int niter,double external_solver_residue,spincolor *external_source)
{
  //Allocate the solution in 128 bit and initialize it. If guess passed copy it.
  spincolor_128 *sol_128=nissa_malloc("sol_128",loc_vol+bord_vol,spincolor_128);
  memset(sol_128,0,loc_vol*sizeof(spincolor_128));
  if(guess!=NULL) quadruple_vector_summassign_double_vector((float_128*)sol_128,(double*)guess,loc_vol*24);
  set_borders_invalid(sol_128);
  
  //internal inverter source
  spincolor *internal_source=nissa_malloc("internal_source",loc_vol+bord_vol,spincolor);
  
  //residue in 128 bit
  spincolor_128 *residue_128=nissa_malloc("residue_128",loc_vol+bord_vol,spincolor_128);
  
  //temporary vector for inverter
  spincolor_128 *temp_128=nissa_malloc("temp",loc_vol+bord_vol,spincolor_128); 
  
  //internal solver stopping condition
  double inner_solver_residue=min_double(1.e-20,external_solver_residue);

  //invert until residue is lower than required
  double current_residue;
  do
    {
      // 1) compute D*sol in quadruple precision
      apply_tmQ2_128(residue_128,conf,kappa,mass,temp_128,sol_128);
      
      // 2) compute the new internal_source = external_source - OP * sol_128, and residue module
      quadruple_vector_subt_from_double_vector((float_128*)residue_128,(double*)external_source,(float_128*)residue_128,loc_vol*24);
      double_vector_from_quadruple_vector((double*)internal_source,(float_128*)residue_128,loc_vol*24);
      current_residue=double_conv_quadruple_vector_glb_scalar_prod((float_128*)residue_128,(float_128*)residue_128,loc_vol*24);
      master_printf("\nExternal loop residue: %lg\n\n",current_residue);

      // 3) if residue not reached, compute the new approximated solution
      if(current_residue>=external_solver_residue)
	{
	  //compute partial sol
	  inv_tmQ2_cg(sol,NULL,conf,kappa,mass,100000,1,inner_solver_residue,internal_source);
	  
	  //add the approximated solution to the total one
	  quadruple_vector_summassign_double_vector((float_128*)sol_128,(double*)sol,loc_vol*24);
	}
    }
  while(current_residue>=external_solver_residue);
  
  //copy the solution in 128 bit to the 64 bit
  double_vector_from_quadruple_vector((double*)sol,(float_128*)sol_128,loc_vol*24);
  
  nissa_free(residue_128);
  nissa_free(temp_128);
  nissa_free(sol_128);
  nissa_free(internal_source);
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
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  read_ildg_gauge_conf(conf,gauge_conf_path);
  
  //prepare the total source
  spincolor *source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  memset(source,0,loc_vol*sizeof(spincolor));
  if(rank==0) source[0][0][0][0]=1;
  set_borders_invalid(source);
  
  //allocate solution
  spincolor *sol=nissa_malloc("sol",loc_vol+bord_vol,spincolor);
  
  inv_tmQ2_cg_128(sol,NULL,conf,kappa,mass,10000,residue,source);
  
  nissa_free(sol);
  nissa_free(source);
  nissa_free(conf);
  
  close_nissa();
  
  return 0;
}
