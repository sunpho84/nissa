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

//This function contract a source with a propagator putting the passed list of operators
void contract_with_source(complex **corr,colorspinspin *prop,int *list_op,colorspinspin *source,int ncontr,int fprop,int rprop)
{
  //Temporary vector for the internal matrices
  dirac_matr t1[ncontr];
  dirac_matr t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Init the second 
      t1[icontr]=base_gamma[list_op[icontr]];
      t2[icontr]=base_gamma[0];

      //Remind that D- rotates as 1+ig5, but D-^-1 rotates as 1-ig5,

      //f1 < 1: do rotation for a quark propagator
      //f1 = 1: do nothing (physical basis already)
      //f1 > 1: do rotation for a (charged) sequential propagator

      if(fprop<1)
	switch(rprop)
	  {
	  case 0: //This is D-^-1
	    dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pminus);
	    dirac_prod(&(t2[icontr]), &Pminus,&(t2[icontr]));
	    break;
	  case 1: //This is D+^-1
	    dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pplus);
	    dirac_prod(&(t2[icontr]), &Pplus,&(t2[icontr]));
	    break;
	  }

      if(fprop>1)
        switch(rprop)
          {
          case 0: //This is D-^-1 on the sequential. t1 is on the sink
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pminus);
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pplus);
            break;
          case 1: //This is D+^-1 on the sequential. t1 is on the sink
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pplus);
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pminus);
            break;
          }
    }

  //Call for the routine which does the real contraction
  trace_g_sdag_g_s(corr,t2,source,t1,prop,ncontr);
}

//This perform the contraction using the so-callled (by Silvano) Chris Micheal trick
//otherwise called "twisted trick"
void contract_with_chris_michael(complex **corr,int *list_op,colorspinspin *s1,colorspinspin *s2,int ncontr,int f1,int r1,int f2,int r2)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      
      //copy the dirac matrices
      t1[icontr]=base_gamma[list_op[icontr]];
      t2[icontr]=base_gamma[0];
      
      //Remind that D- rotates as 1+ig5, but D-^-1 rotates as 1-ig5,
      //moreover (D^-1)^dagger rotates again as 1+ig5 (pweee!!!)

      //f1 < 1: do rotation for a quark propagator
      //f1 = 1: do nothing (physical basis already)
      //f1 > 1: do rotation for a (charged) sequential propagator

      if(f1<1)
        switch(r1)
          {
          case 0: //This is (D-^-1)^dag
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pplus);
            dirac_prod(&(t2[icontr]), &Pplus,&(t2[icontr]));
            break;
          case 1: //This is (D+^-1)^dag
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pminus);
            dirac_prod(&(t2[icontr]), &Pminus,&(t2[icontr]));
            break;
          }

      if(f1>1)
        switch(r1)
          {
          case 0: //This is (D-^-1)^dag
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pplus);
            dirac_prod(&(t2[icontr]), &Pminus,&(t2[icontr]));
            break;
          case 1: //This is (D+^-1)^dag
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pminus);
            dirac_prod(&(t2[icontr]), &Pplus,&(t2[icontr]));
            break;
          }

      if(f2<1)
        switch(r2)
          {
          case 0: //This is D-^-1
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pminus);
            dirac_prod(&(t1[icontr]), &Pminus,&(t1[icontr]));
            break;
          case 1: //This is D+^-1
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pplus);
            dirac_prod(&(t1[icontr]), &Pplus,&(t1[icontr]));
            break;
          }

      if(f2>1)
        switch(r2)
          {
          case 0: //This is D-^-1
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pminus);
            dirac_prod(&(t1[icontr]), &Pplus,&(t1[icontr]));
            break;
          case 1: //This is D+^-1
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pplus);
            dirac_prod(&(t1[icontr]), &Pminus,&(t1[icontr]));
            break;
          }
    }
  
  //Call for the routine which does the real contraction
  trace_g_sdag_g_s(corr,t2,s2,t1,s1,ncontr);
}
	  
//print all the passed contractions to the file
void print_contractions_to_file(FILE *fout,int ncontr,int *op,complex **contr,int twall,const char *tag)
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
  double tic;
  if(debug)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      tic=MPI_Wtime();
    }

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

  //Read the time location of the source
  int twall;
  read_str_int("TWall",&twall);

  //Read the source
  char path[1024];
  read_str_str("Source",path,1024);

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
  
  //allocate the source and the prop
  colorspinspin *source=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  colorspinspin *prop=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  if((source==NULL || prop==NULL ) && rank==0)
    {
      fprintf(stderr,"Error! Not enoug memory to allocate source and spinor!\n");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  //read the source
  read_colorspinspin(source,path,NULL);
  
  //allocation of the chromo-spinor
  colorspinspin *ch_prop=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  if(ch_prop==NULL && rank==0)
    {
      fprintf(stderr,"Error! Not enoug memory to allocate the chromo spinor!\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  //Inizialization of output
  char outfile[1024];
  FILE *fout=NULL;
  read_str_str("Output",outfile,1024);
  if(rank==0) fout=fopen(outfile,"w");
  
  //Read the number of the props
  int nprop;
  read_str_int("NProp",&nprop);
  
  //Loop over propagators
  int phys,r;
  char tag[1024];
  for(int iprop=0;iprop<nprop;iprop++)
    {
      //Read the name, mass, theta and r for the list of the propagators
      read_str(path,1024);
      read_str(tag,1024);
      read_int(&phys);
      read_int(&r);
      
      if(debug && rank==0)
	printf(" prop.%d %s tagged as '%s', phys=%d r=%d\n",iprop,path,tag,phys,r);
      
      //Read the propagator
      read_colorspinspin(prop,path,NULL);
      
      //apply the chromo magnetic operator to the second spinor
      unsafe_apply_chromo_operator_to_colorspinspin(ch_prop,Pmunu,prop);
      
      if(rank==0) fprintf(fout," # %s r=%d\n",tag,r);
      
      contract_with_source(contr,prop,op,source,ncontr,phys,r);
      print_contractions_to_file(fout,ncontr,op,contr,twall,"");

      contract_with_source(ch_contr,ch_prop,ch_op,source,nch_contr,phys,r);
      print_contractions_to_file(fout,nch_contr,ch_op,ch_contr,twall,"CHROMO-");

      if(r==0)
	{
	  contract_with_chris_michael(contr,op,prop,prop,ncontr,phys,r,phys,r);
	  print_contractions_to_file(fout,ncontr,op,contr,twall,"CHRIS-");
	  
	  contract_with_chris_michael(ch_contr,ch_op,ch_prop,prop,nch_contr,phys,r,phys,r);
	  print_contractions_to_file(fout,nch_contr,ch_op,ch_contr,twall,"CHROMO-CHRIS-");
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
