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

int ncontr,nch_contr;
complex **contr,**ch_contr;
int *op,*ch_op;

char gaugeconf_file[1024];
quad_su3 *conf;
as2t_su3 *Pmunu;

colorspinspin *source,*ch_prop,***S;
spincolor *inv_source,*reco_solution[2],**QQ;

int seed,nsource,noise_type;

int nmass;
double kappa,*mass;

double residue,minimal_residue;
int stopping_criterion,niter_max;
char outfile[1024];

double tot_time=0,inv_time=0,contr_time=0;
int ncontr_tot,ninv_tot;

const double rad2=1.414213562373095048801688724209;

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
void print_bubbles_to_file(FILE *fout,int ncontr,int *op,complex **contr,const char *tag)
{
  if(rank==0)
    {
      fprintf(fout,"\n");

      for(int icontr=0;icontr<ncontr;icontr++)
	{
	  fprintf(fout," # %s%s\n\n",tag,gtag[op[icontr]]);
	  for(int t=0;t<glb_size[0];t++)
	    fprintf(fout,"%+016.16g\t%+016.16g\n",contr[icontr][t][0],contr[icontr][t][1]);
	  fprintf(fout,"\n");
	}
      fprintf(fout,"\n");
    }
}

void initialize_bubbles(char *input_path)
{
  //take initial time
  tot_time=-take_time();

  open_input(input_path);

  //Read the volume
  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));

  //Init the MPI grid 
  init_grid();

  //Read the number of contractions
  read_str_int("NContr",&ncontr);
  if(rank==0) printf("Number of contractions: %d\n",ncontr);

  //Initialize the list of correlations and the list of operators
  contr=(complex**)malloc(sizeof(complex*)*ncontr);
  contr[0]=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);
  op=(int*)malloc(sizeof(int)*ncontr);
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      contr[icontr]=contr[0]+icontr*glb_size[0];
      //Read the operator
      read_int(&(op[icontr]));
      
      if(rank==0 && debug) printf(" contr.%d %d\n",icontr,op[icontr]);
    }
  
  //Read the number of contractions with insertion of the chromo-magnetic operator
  read_str_int("NChromoContr",&nch_contr);
  if(rank==0) printf("Number of chromo-contractions: %d\n",nch_contr);
  
  //Initialize the list of chromo correlations and the list of operators
  //contiguous allocation
  ch_contr=(complex**)malloc(sizeof(complex*)*nch_contr);
  ch_contr[0]=(complex*)malloc(sizeof(complex)*nch_contr*glb_size[0]); 
  ch_op=(int*)malloc(sizeof(int)*nch_contr);
  for(int ich_contr=0;ich_contr<nch_contr;ich_contr++)
    {
      ch_contr[ich_contr]=ch_contr[0]+ich_contr*glb_size[0];
      //Read the operator
      read_int(&(ch_op[ich_contr]));
      
      if(rank==0 && debug) printf(" chromo contr.%d %d\n",ich_contr,ch_op[ich_contr]);
    }
  
  //reading of gauge conf and computation of Pmunu
  read_str_str("GaugeConf",gaugeconf_file,1024);
  conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"gauge conf");
  read_local_gauge_conf(conf,gaugeconf_file);
  communicate_gauge_borders(conf);
  communicate_gauge_edges(conf);

  //calculate plaquette, Pmunu
  double gplaq=global_plaquette(conf);
  if(rank==0) printf("plaq: %.10g\n",gplaq);
  Pmunu=allocate_as2t_su3(loc_vol,"clover term");
  Pmunu_term(Pmunu,conf);
  
  //allocate the source and the chromo-prop
  source=allocate_colorspinspin(loc_vol,"source");
  inv_source=allocate_spincolor(loc_vol+loc_bord,"inv_source");
  ch_prop=allocate_colorspinspin(loc_vol,"chromo-prop");
  for(int r=0;r<2;r++) reco_solution[r]=allocate_spincolor(loc_vol,"reco_solution");

  //read the seed
  read_str_int("Seed",&seed);
  init_random(seed);
  //read the number of sources
  read_str_int("NSource",&nsource);
  //read the noise type
  read_str_int("NoiseType",&noise_type);

  //Read kappa
  read_str_double("Kappa",&kappa);

  //Read the number of masses and allocate spinors for the cgmms
  read_str_int("NMass",&nmass);
  mass=(double*)malloc(sizeof(double)*nmass);
  S=(colorspinspin***)malloc(nmass*sizeof(colorspinspin**));
  QQ=(spincolor**)malloc(nmass*sizeof(spincolor*));
  for(int imass=0;imass<nmass;imass++)
    {
      S[imass]=(colorspinspin**)malloc(2*sizeof(colorspinspin*));
      for(int r=0;r<2;r++) S[imass][r]=allocate_colorspinspin(loc_vol,"S");
      QQ[imass]=allocate_spincolor(loc_vol+loc_bord,"QQ");
      read_double(&(mass[imass]));
    }
    
  //Residue
  read_str_double("Residue",&residue);
  //Stopping criterion
  get_stopping_criterion(&stopping_criterion,&minimal_residue);
  //Number of iterations
  read_str_int("NiterMax",&niter_max);

  //Read path of output
  read_str_str("Output",outfile,1024);
  
  close_input();
}

void generate_volume_source()
{ //reset
  memset(source,0,sizeof(colorspinspin)*loc_vol);
  
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    for(int ic=0;ic<3;ic++)
      { //real part
	if(noise_type>=2) source[loc_site][ic][0][0][0]=pm_one(loc_site)/rad2;
	else source[loc_site][ic][0][0][0]=noise_type;
	//imaginary part
	if(noise_type==4) source[loc_site][ic][0][0][1]=pm_one(loc_site)/rad2;
        
	for(int d=1;d<4;d++) //copy the other three dirac indexes
	  memcpy(source[loc_site][ic][d][d],source[loc_site][ic][0][0],sizeof(complex));
      }
}

void calculate_S()
{
  inv_time-=take_time();

  //Invert
  for(int idso=0;idso<4;idso++)
    {
      //prepare the source
      for(int ivol=0;ivol<loc_vol;ivol++)
	{
	  get_spincolor_from_colorspinspin(inv_source[ivol],source[ivol],idso);
	  
	  //put the g5
	  for(int idsi=0;idsi<2;idsi++)
	    for(int icol=0;icol<3;icol++)
	      for(int ri=0;ri<2;ri++)
		inv_source[ivol][idsi][icol][ri]=-inv_source[ivol][idsi][icol][ri];
	}
      
      inv_Q2_cgmms(QQ,inv_source,NULL,conf,kappa,mass,nmass,niter_max,residue,minimal_residue,stopping_criterion);
      ninv_tot++;
      //put the solution inside the S vector
      for(int imass=0;imass<nmass;imass++) 
	{
	  reconstruct_doublet(reco_solution[0],reco_solution[1],QQ[imass],conf,kappa,mass[imass]);
	  for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
	    for(int ivol=0;ivol<loc_vol;ivol++)
	      put_spincolor_into_colorspinspin(S[imass][r][ivol],reco_solution[r][ivol],idso);
	}
    }
  
  //rotate to physical basis
  for(int imass=0;imass<nmass;imass++) //put the (1+ig5)/sqrt(2) factor
    for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
      rotate_vol_colorspinspin_to_physical_basis(S[imass][r],!r,!r);

  inv_time+=take_time();
}

void calculate_all_contractions(int isource)
{
  contr_time-=take_time();

  char outp[1024];
  sprintf(outp,"%s%d",outfile,isource);
  FILE *fout=open_text_file_for_output(outfile);;
  
  for(int imass=0;imass<nmass;imass++)
    for(int r=0;r<2;r++)
      {
	//apply the chromo magnetic operator to the second spinor
	unsafe_apply_chromo_operator_to_colorspinspin(ch_prop,Pmunu,S[imass][r]);
	
	if(rank==0) fprintf(fout," # mass=%g r=%d\n",mass[imass],r);
	
	contract_with_source(contr,S[imass][r],op,source,ncontr);
	print_bubbles_to_file(fout,ncontr,op,contr,"");

	contract_with_source(ch_contr,ch_prop,ch_op,source,nch_contr);
	print_bubbles_to_file(fout,nch_contr,ch_op,ch_contr,"CHROMO-");
	
	ncontr_tot+=nch_contr+ncontr;
	
	if(r==0)
	  {
	    contract_with_chris_michael(contr,op,S[imass][r],S[imass][r],ncontr);
	    print_bubbles_to_file(fout,ncontr,op,contr,"CHRIS-");
	    
	    contract_with_chris_michael(ch_contr,ch_op,ch_prop,S[imass][r],nch_contr);
	    print_bubbles_to_file(fout,nch_contr,ch_op,ch_contr,"CHROMO-CHRIS-");
	    
	    ncontr_tot+=nch_contr+ncontr;
	  }
      }

  contr_time+=take_time();

  if(rank==0) fclose(fout);
}

void close_bubbles()
{

  //take final time
  tot_time+=take_time();

  if(rank==0)
    {
     printf("\n");
     printf("Total time: %g, of which:\n",tot_time);
     printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
     printf(" - %02.2f%s to perform %d contr. (%2.2gs avg)\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
    }

  ///////////////////////////////////////////
  
  if(rank==0) printf("\nEverything ok, exiting!\n");
  
  close_appretto();
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  initialize_bubbles(arg[1]);

  //loop over the sources
  for(int isource=0;isource<nsource;isource++)
    {
      generate_volume_source(source);
      calculate_S();
      calculate_all_contractions(isource);
    }

  close_bubbles();

  return 0;
}
