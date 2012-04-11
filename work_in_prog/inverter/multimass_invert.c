#include "nissa.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();

  if(narg<2) crash("Use: %s input_file",arg[0]);

  open_input(arg[1]);

  //Init the MPI grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  //masses
  int nmass;
  read_str_int("NMass",&(nmass));
  double m[nmass];
  read_str_double("masses",&(m[0]));
  for(int imass=1;imass<nmass;imass++) read_double(&(m[imass]));
  //kappa
  double kappa;
  read_str_double("kappa",&(kappa));
  //Initialize the gauge configuration and read the path
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  char gauge_file[1024];
  read_str_str("GaugeConf",gauge_file,1024);
  //read the thetas in multiples of 2*pi
  double theta[4];
  read_str_double("ThetaTXYZ",&(theta[0]));
  read_double(&(theta[1]));
  read_double(&(theta[2]));
  read_double(&(theta[3]));
  master_printf("Thetas: %f %f %f %f\n",theta[0],theta[1],theta[2],theta[3]);

  //load the configuration, put boundaries condition and communicate borders
  read_ildg_gauge_conf(conf,gauge_file);
  put_boundaries_conditions(conf,theta,1,0);
  communicate_lx_quad_su3_borders(conf);

  //initialize the source
  spincolor *source=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+bord_vol));
  spincolor *source_reco=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  
  //initialize solution
  spincolor *solution[nmass];
  
  for(int imass=0;imass<nmass;imass++) solution[imass]=nissa_malloc("solution",loc_vol+bord_vol,spincolor);
  
  double residue,minimal_residue;
  read_str_double("Residue",&residue);
  int stopping_criterion=numb_known_stopping_criterion;
  char str_stopping_criterion[1024];
  read_str_str("StoppingCriterion",str_stopping_criterion,1024);
  int isc=0;
  do
    {
      if(strcasecmp(list_known_stopping_criterion[isc],str_stopping_criterion)==0) stopping_criterion=isc;
      isc++;
    }
  while(isc<numb_known_stopping_criterion && stopping_criterion==numb_known_stopping_criterion);

  if(stopping_criterion==numb_known_stopping_criterion && rank==0)
    {
      fprintf(stderr,"Unknown stopping criterion: %s\n",str_stopping_criterion);
      fprintf(stderr,"List of known stopping criterions:\n");
      for(int isc=0;isc<numb_known_stopping_criterion;isc++) fprintf(stderr," %s\n",list_known_stopping_criterion[isc]);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  if(stopping_criterion==sc_standard) read_str_double("MinimalResidue",&minimal_residue);
  
  int nitermax;
  read_str_int("NiterMax",&nitermax);
  
  int nsource;
  read_str_int("NSource",&nsource);

  char source_file[1024];
  char output_file[1024];

  for(int isource=0;isource<nsource;isource++)
    {
      master_printf("Source #%d\n",isource);

      read_str_str("Source",source_file,1024);
      read_spincolor(source,source_file);
      
      read_str_str("Output",output_file,1024);
      
      //multiply the source by gamma5
      for(int X=0;X<loc_vol;X++)
	for(int d=2;d<4;d++)
	  for(int c=0;c<3;c++)
	    for(int r=0;r<2;r++)
	      source[X][d][c][r]=-source[X][d][c][r];
      
      communicate_lx_spincolor_borders(source);
      
      ///////////////////////////////////////////
      
      //take initial time                                                                                                        
      double tinv=-take_time;
      inv_Q2_cgmms(solution,source,NULL,conf,kappa,m,nmass,nitermax,residue,minimal_residue,stopping_criterion);
      tinv+=take_time();
      
      master_printf("\nTotal time elapsed: %f s\n",tinv);
      
      for(int imass=0;imass<nmass;imass++)
	{
	  apply_Q2(source_reco,solution[imass],conf,kappa,m[imass],NULL,NULL,NULL);
	  
	  for(int ud=0;ud<2;ud++)
	    {
	      if(ud==0) apply_Q(source_reco,solution[imass],conf,kappa,m[imass]);
	      else apply_Q(source_reco,solution[imass],conf,kappa,m[imass]);
	      
	      char output_full[1024];
	      sprintf(output_full,"%s.%02d.%d",output_file,imass,ud);
	      write_spincolor(output_full,source_reco,64);
	    }
	}
    }
  
  close_input();
  
  ///////////////////////////////////////////
  
  for(int imass=0;imass<nmass;imass++) nissa_free(solution[imass]);
  nissa_free(source);
  nissa_free(source_reco);

  nissa_free(conf);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
