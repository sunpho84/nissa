#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));

  int nmass;
  read_str_int("NMass",&(nmass));

  double m[nmass];
  double kappa;
  read_str_double("masses",&(m[0]));
  for(int imass=1;imass<nmass;imass++) read_double(&(m[imass]));

  read_str_double("kappa",&(kappa));

  //Init the MPI grid 
  init_grid();

  //Initialize the gauge configuration and read the path
  quad_su3 *conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+loc_bord));
  char gauge_file[1024];
  read_str_str("GaugeConf",gauge_file,1024);
  
  //read the thetas in multiples of 2*pi
  double theta[4];
  read_str_double("ThetaTXYZ",&(theta[0]));
  read_double(&(theta[1]));
  read_double(&(theta[2]));
  read_double(&(theta[3]));
  if(rank==0)
    printf("Thetas: %f %f %f %f\n",theta[0],theta[1],theta[2],theta[3]);

  //load the configuration, put boundaries condition and communicate borders
  read_local_gauge_conf(conf,gauge_file);
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);

  //initialize the source
  spincolor *source=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  spincolor *source_reco=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  
  //initialize solution
  spincolor *solution[nmass];
  
  for(int imass=0;imass<nmass;imass++) solution[imass]=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));

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
      printf("Source #%d\n",isource);

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
      double tic;
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
      inv_Q2_cgmms(solution,source,NULL,conf,kappa,m,nmass,nitermax,residue,minimal_residue,stopping_criterion);
      
      MPI_Barrier(cart_comm);
      double tac=MPI_Wtime();
      if(rank==0)
	printf("\nTotal time elapsed: %f s\n",tac-tic);
      
      for(int imass=0;imass<nmass;imass++)
	{
	  apply_Q2(source_reco,solution[imass],conf,kappa,m[imass],NULL,NULL,NULL);
	  
	  //printing
	  double truered,loc_truered=0;
	  for(int loc_site=0;loc_site<loc_vol;loc_site++)
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		{
		  double tempr=source[loc_site][id][ic][0]-source_reco[loc_site][id][ic][0];
		  double tempi=source[loc_site][id][ic][1]-source_reco[loc_site][id][ic][1];
		  loc_truered+=tempr*tempr+tempi*tempi;
		}
	  
	  if(rank_tot>0) MPI_Allreduce(&loc_truered,&truered,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);                                 
	  else truered=loc_truered;
	  
	  if(rank==0) printf("Residue for mass %d: %g\n",imass,truered);
	  
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
  
  for(int imass=0;imass<nmass;imass++) free(solution[imass]);
  free(source);
  free(source_reco);

  free(conf);
  
  ///////////////////////////////////////////
  
  close_appretto();
  
  return 0;
}
