#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

/*
            	 
           .....           
         ..     ..         
        .         .        
   op  X           .  S
        .         .        
         ..     ..         
           .....           
                      
*/

//Read the stopping criterion
void get_stopping_criterion(int *stopping_criterione,double *minimal_residue)
{
  int isc=0;
  stopping_criterion=numb_known_stopping_criterion;
  char str_stopping_criterion[1024];
  read_str_str("StoppingCriterion",str_stopping_criterion,1024);
  isc=0;
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
  if(stopping_criterion==sc_standard) read_str_double("MinimalResidue",minimal_residue);
}

//This function contract a source with a propagator putting the passed list of operators
void contract_with_source(complex **corr,colorspinspin *prop,int *list_op,colorspinspin *source,int ncontr)
{
  //Temporary vector for the internal matrices
  dirac_matr t1[ncontr];
  dirac_matr t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Init the second 
      t1[icontr]=base_gamma[list_op[icontr]];
      t2[icontr]=base_gamma[0];
    }

  //Call for the routine which does the real contraction
  trace_g_sdag_g_s(corr,t2,source,t1,prop,ncontr);
}

//This perform the contraction using the so-callled (by Silvano) Chris Micheal trick
//otherwise called "twisted trick"
void contract_with_chris_michael(complex **corr,int *list_op,colorspinspin *s1,colorspinspin *s2,int ncontr)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {      
      //copy the dirac matrices
      t1[icontr]=base_gamma[list_op[icontr]];
      t2[icontr]=base_gamma[0];
    }
  
  //Call for the routine which does the real contraction
  trace_g_sdag_g_s(corr,t2,s2,t1,s1,ncontr);
}
	  
//print all the passed contractions to the file
void print_bubbles_to_file(FILE *fout,int ncontr,int *op,complex **contr,int twall,const char *tag)
{
  if(rank==0)
    {
      fprintf(fout,"\n");

      for(int icontr=0;icontr<ncontr;icontr++)
	{
	  fprintf(fout," # %s%s\n\n",tag,gtag[op[icontr]]);
	  for(int tempt=0;tempt<glb_size[0];tempt++)
	    {
	      int t=tempt+twall;
	      if(t>=glb_size[0]) t-=glb_size[0];

	      fprintf(fout,"%+016.16g\t%+016.16g\n",contr[icontr][t][0],contr[icontr][t][1]);
	    }
	  fprintf(fout,"\n");
	}
      fprintf(fout,"\n");
    }
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_appretto();

  //take initial time
  double tot_time=-take_time();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  //Read the volume
  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));

  //Init the MPI grid 
  init_grid();

  //Read the number of contractions
  int ncontr;
  read_str_int("NContr",&ncontr);
  if(rank==0) printf("Number of contractions: %d\n",ncontr);

  //Initialize the list of correlations and the list of operators
  complex **contr=(complex**)malloc(sizeof(complex*)*ncontr);
  contr[0]=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
  int *op=(int*)malloc(sizeof(int)*ncontr);
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      contr[icontr]=contr[0]+icontr*glb_size[0];
      //Read the operator
      read_int(&(op[icontr]));
      
      if(rank==0 && debug) printf(" contr.%d %d\n",icontr,op[icontr]);
    }
  
  //Read the number of contractions with insertion of the chromo-magnetic operator
  int nch_contr;
  read_str_int("NChromoContr",&nch_contr);
  if(rank==0) printf("Number of chromo-contractions: %d\n",nch_contr);
  
  //Initialize the list of chromo correlations and the list of operators
  //contiguous allocation
  complex **ch_contr=(complex**)malloc(sizeof(complex*)*nch_contr);
  ch_contr[0]=(complex*)malloc(sizeof(complex)*nch_contr*glb_size[0]); 
  int *ch_op=(int*)malloc(sizeof(int)*nch_contr);
  for(int ich_contr=0;ich_contr<nch_contr;ich_contr++)
    {
      ch_contr[ich_contr]=ch_contr[0]+ich_contr*glb_size[0];
      //Read the operator
      read_int(&(ch_op[ich_contr]));
      
      if(rank==0 && debug) printf(" chromo contr.%d %d\n",ich_contr,ch_op[ich_contr]);
    }
  
  //reading of gauge conf and computation of Pmunu
  char gaugeconf_file[1024];
  read_str_str("GaugeConf",gaugeconf_file,1024);
  quad_su3 *gauge_conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+loc_bord+loc_edge));
  as2t_su3 *Pmunu=(as2t_su3*)malloc(sizeof(as2t_su3)*loc_vol);
  if((gauge_conf==NULL || Pmunu==NULL) && rank==0)
    {
      fprintf(stderr,"Error! Not enoug memory to allocate the gauge conf and the chromo term!\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  read_local_gauge_conf(gauge_conf,gaugeconf_file);
  communicate_gauge_borders(gauge_conf);
  communicate_gauge_edges(gauge_conf);
  double gplaq=global_plaquette(gauge_conf);
  if(rank==0) printf("plaq: %.10g\n",gplaq);
  Pmunu_term(Pmunu,gauge_conf);
  free(gauge_conf);
  
  //allocate the source and the chromo-prop
  colorspinspin *source=allocate_colorspinspin(loc_vol,"source");
  spincolor *inv_source=allocate_spincolor(loc_vol+loc_bord,"inv_source");
  colorspinspin *ch_prop=allocate_colorspinspin(loc_vol,"chromo-prop");
  spincolor *reco_solution[2];
  for(int r=0;r<2;r++) reco_solution[r]=llocate_spincolor(loc_vol,"reco_solution");

  //read the seed
  int seed;
  read_str_int("Seed",&seed);
  
  //read the number of sources
  int nsource;
  read_str_int("NSource",&nsource);

  //Read kappa
  double kappa;
  read_str_double("Kappa",&kappa);

  //Read the number of masses and allocate spinors for the cgmms
  int nmass;
  read_str_int("NMass",&nmass);
  double mass[nmass];
  colorspinspin **S[2]=(colorspinspin**)malloc(nmass*sizeof(colorspinspin*));
  spincolor **QQ=(spincolor**)malloc(nmass*sizeof(spincolor*));
  for(int imass=0;imass<nmass;imass++)
    {
      for(int r=0;r<2;r++) S[r]=allocate_colorspinspin(loc_vol,"S");
      QQ[imass]=allocate_spincolor(loc_vol+loc_bord,"QQ");
      read_double(&(mass[imass]));
    }
    
  //Residue
  double residue,minimal_residue;
  read_str_double("Residue",&stopping_residue);
  //Stopping criterion
  int stopping_criterion;
  get_stopping_criterion(&stopping_criterion,&minimal_residue);
  //Number of iterations
  read_str_int("NiterMax",&niter_max);

  //Inizialization of output
  char outfile[1024];
  read_str_str("Output",outfile,1024);
  FILE *fout=open_text_file_for_output(outfile);
  
  //loop over the sources
  for(int isource=0;isource<nsource;isource++)
    {
      generate_volume_source(source);
  
      //Invert
      for(int idso=0;idso<4;idso++)
	{
	  
	  //prepare the source
	  for(int ivol=0;ivol<loc_vol;ivol++)
	    {
	      get_spincolor_from_colorspinspin(inv_source[ivol],source[ivol],idso);

	      //put the g5
	      for(int idsi=0;idsi<2;idsi++)
		for(int ri=0;ri<2;ri++)
		  inv_source[ivol][idsi][ri]=-inv_source[ivol][idsi][ri];
	    }

	  inv_Q2_cgmms(QQ,inv_source,NULL,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	  //put the solution inside the S vector
	  for(int imass=0;imass<nmass;imass++) 
	    {
	      reconstruct_doublet(reco_solution[0],reco_solution[1],QQ[imass],conf,kappa,mass[imass]);
	      for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
                for(int ivol=0;ivol<loc_vol;ivol++)
                  put_spincolor_into_colorspinspin(S0[r][imass][ivol],reco_solution[r][ivol],idso);
	    }
	}
      
      

      //apply the chromo magnetic operator to the second spinor
      unsafe_apply_chromo_operator_to_colorspinspin(ch_prop,Pmunu,prop);
      
      if(rank==0) fprintf(fout," # %s r=%d\n",tag,r);
      
      contract_with_source(contr,prop,op,source,ncontr,phys,r);
      print_bubbles_to_file(fout,ncontr,op,contr,twall,"");

      contract_with_source(ch_contr,ch_prop,ch_op,source,nch_contr,phys,r);
      print_bubbles_to_file(fout,nch_contr,ch_op,ch_contr,twall,"CHROMO-");

      if(r==0)
	{
	  contract_with_chris_michael(contr,op,prop,prop,ncontr,phys,r,phys,r);
	  print_bubbles_to_file(fout,ncontr,op,contr,twall,"CHRIS-");
	  
	  contract_with_chris_michael(ch_contr,ch_op,ch_prop,prop,nch_contr,phys,r,phys,r);
	  print_bubbles_to_file(fout,nch_contr,ch_op,ch_contr,twall,"CHROMO-CHRIS-");
	}
    }

  close_input();

  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();
      if(rank==0) printf("Time elapsed in doing %d contractions: %g s\n",nprop*(ncontr+nch_contr),tac-tic);
    }

  ///////////////////////////////////////////
  
  if(rank==0)
    {
      printf("\nEverything ok, exiting!\n");
      fclose(fout);
    }
  
  free(ch_prop);
  free(Pmunu);
  
  free(prop);
  free(source);

  free(op);
  free(contr[0]);
  free(contr);

  close_appretto();

  return 0;
}
