#include "appretto.h"

//gauge info
int ngauge_conf;
char conf_path[1024];
quad_su3 *conf,*sme_conf;
double kappa;

int nmass;
double *mass;
double put_theta[4],old_theta[4]={0,0,0,0};

//source data
int seed,noise_type,twall;
double rad2=1.414213562373095048801688724209;
spincolor *source;
colorspinspin *original_source;

//smearing parameters
double jacobi_kappa,ape_alpha;
int *jacobi_niter,ape_niter;
int nsm_lev;

//vectors for the spinor data
int npropS0;
colorspinspin **S0[2];
spincolor **cgmms_solution,*reco_solution[2];

//cgmms inverter parameters
double stopping_residue;
double minimal_residue;
int stopping_criterion;
int niter_max;

//two points contractions
int ncontr_2pts;
complex *contr_2pts;
int *op1_2pts,*op2_2pts;
char outfile_2pts[1024];

//timings
int ninv_tot=0,ncontr_tot=0;
double tot_time=0,inv_time=0,contr_time=0;

//This function takes care to make the revert on the FIRST spinor, putting the needed gamma5
void meson_two_points(complex *corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Put the two gamma5 needed for the revert of the first spinor
      dirac_prod(&(t1[icontr]), &(base_gamma[list_op1[icontr]]),&(base_gamma[5]));
      dirac_prod(&(t2[icontr]), &(base_gamma[5]),&(base_gamma[list_op2[icontr]]));
    }
  
  //Call the routine which does the real contraction
  trace_g_sdag_g_s(corr,t1,s1,t2,s2,ncontr);
}

//Generate the source for the dirac index
void generate_source()
{ //reset
  memset(original_source,0,sizeof(colorspinspin)*loc_vol);
  
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    if(glb_coord_of_loclx[loc_site][0]==twall)
      for(int ic=0;ic<3;ic++)
	{ //real part
	  if(noise_type>=2) original_source[loc_site][ic][0][0][0]=pm_one(loc_site)/rad2;
	  else original_source[loc_site][ic][0][0][0]=noise_type;
	  //imaginary part
	  if(noise_type==4) original_source[loc_site][ic][0][0][1]=pm_one(loc_site)/rad2;
	  
	  for(int d=1;d<4;d++) //copy the other three dirac indexes
	    memcpy(original_source[loc_site][ic][d][d],original_source[loc_site][ic][0][0],sizeof(complex));
	}
}

//Parse all the input file
void initialize_semileptonic(char *input_path)
{
  open_input(input_path);

  // 1) Read information about the gauge conf
  
  //Read the volume
  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  //Init the MPI grid 
  init_grid(); 
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the source and masses
  
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  init_random(seed);
  //Read the noise type
  read_str_int("NoiseType",&noise_type);
  //Read the masses
  read_list_of_doubles("NMass",&nmass,&mass);
  
  // 3) Smearing parameters
  
  //Smearing parameters
  read_str_double("ApeAlpha",&ape_alpha);
  read_str_int("ApeNiter",&ape_niter);
  read_str_double("JacobiKappa",&jacobi_kappa);
  read_list_of_ints("JacobiNiters",&nsm_lev,&jacobi_niter);
  if(jacobi_niter[0]!=0)
    crash("Error, jacobi level of smearing 0 have to be null, obtained: %d!",jacobi_niter[0]);
  for(int iter=1;iter<nsm_lev;iter++)
    if(jacobi_niter[iter]<jacobi_niter[iter-1])
      crash("Error, jacobi level %d minor than %d (%d, %d)!",iter,iter-1,jacobi_niter[iter],jacobi_niter[iter-1]);
  
  // 4) Info about the inverter

  //Residue
  read_str_double("Residue",&stopping_residue);
  //Stopping criterion
  stopping_criterion=numb_known_stopping_criterion;
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
      
  //Number of iterations
  read_str_int("NiterMax",&niter_max);
  
  // 5) contraction list for two points
  
  read_str_int("NContrTwoPoints",&ncontr_2pts);
  contr_2pts=appretto_malloc("contr_2pts",ncontr_2pts*glb_size[0],complex);
  op1_2pts=appretto_malloc("op1_2pts",ncontr_2pts,int);
  op2_2pts=appretto_malloc("op2_2pts",ncontr_2pts,int);
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1_2pts[icontr]));
      read_int(&(op2_2pts[icontr]));

      if(rank==0 && debug) printf(" contr.%d %d %d\n",icontr,op1_2pts[icontr],op2_2pts[icontr]);
    }
  
  read_str_int("NGaugeConf",&ngauge_conf);  
  
  ////////////////////////////////////// end of input reading/////////////////////////////////
  
  //allocate gauge conf and all the needed spincolor and colorspinspin
  conf=appretto_malloc("or_conf",loc_vol+loc_bord,quad_su3);
  sme_conf=appretto_malloc("sm_conf",loc_vol+loc_bord,quad_su3);
  
  //Allocate all the S0 colorspinspin vectors
  npropS0=nmass;
  S0[0]=appretto_malloc("S0[0]",npropS0,colorspinspin*);
  S0[1]=appretto_malloc("S0[1]",npropS0,colorspinspin*);
  for(int iprop=0;iprop<npropS0;iprop++)
    {
      S0[0][iprop]=appretto_malloc("S0[0][iprop]",loc_vol,colorspinspin);
      S0[1][iprop]=appretto_malloc("S0[1][iprop]",loc_vol,colorspinspin);
    }
  
  //Allocate nmass spincolors, for the cgmms solutions
  cgmms_solution=appretto_malloc("cgmms_solution",nmass,spincolor*);
  for(int imass=0;imass<nmass;imass++) cgmms_solution[imass]=appretto_malloc("cgmms_solution[imass]",loc_vol+loc_bord,spincolor);
  reco_solution[0]=appretto_malloc("reco_sol[0]",loc_vol+loc_bord,spincolor);
  reco_solution[1]=appretto_malloc("reco_sol[1]",loc_vol+loc_bord,spincolor);
  
  //Allocate one spincolor for the source
  source=appretto_malloc("source",loc_vol+loc_bord,spincolor);
  original_source=appretto_malloc("orig_source",loc_vol,colorspinspin);
}

//load the conf, smear it and put boundary cond
void load_gauge_conf()
{
  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_gauge_conf(conf,conf_path);
  //prepare the smerded version and calculate plaquette
  ape_smearing(sme_conf,conf,ape_alpha,ape_niter);
  communicate_gauge_borders(conf);
  communicate_gauge_borders(sme_conf);
  
  double gplaq=global_plaquette(conf);
  if(rank==0) printf("plaq: %.18g\n",gplaq);
  gplaq=global_plaquette(sme_conf);
  if(rank==0) printf("smerded plaq: %.18g\n",gplaq);

  //Put the anti-periodic condition on the temporal border
  old_theta[0]=0;
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,0);
}

//Finalization
void close_semileptonic()
{
  appretto_free(original_source);appretto_free(source);
  appretto_free(reco_solution[1]);appretto_free(reco_solution[0]);
  for(int imass=0;imass<nmass;imass++) appretto_free(cgmms_solution[imass]);
  appretto_free(cgmms_solution);
  for(int r=0;r<2;r++)
    {
      for(int iprop=0;iprop<npropS0;iprop++)
	appretto_free(S0[r][iprop]);
      appretto_free(S0[r]);
    }
  appretto_free(sme_conf);
  appretto_free(conf);
  appretto_free(contr_2pts);
  appretto_free(op1_2pts);
  appretto_free(op2_2pts);

  if(rank==0)
    {
      printf("\n");
      printf("Total time: %g, of which:\n",tot_time);
      printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
      printf(" - %02.2f%s to perform %d contr. (%2.2gs avg)\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
    }
  
  close_appretto();
}

//calculate the standard propagators
void calculate_S0(int sm_lev_sour)
{
  for(int id=0;id<4;id++)
    { //loop over the source dirac index
      for(int ivol=0;ivol<loc_vol;ivol++)
	{
	  get_spincolor_from_colorspinspin(source[ivol],original_source[ivol],id);
	  //put the g5
	  for(int id1=2;id1<4;id1++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic][ri]*=-1;
	}
      
      //smerd the source if needed
      if(sm_lev_sour!=0) jacobi_smearing(source,source,sme_conf,jacobi_kappa,jacobi_niter[sm_lev_sour]);
      
      double part_time=-take_time();
      communicate_lx_spincolor_borders(source);
      inv_Q2_cgmms(cgmms_solution,source,NULL,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
      part_time+=take_time();ninv_tot++;inv_time+=part_time;
      if(rank==0) printf("Finished the inversion of S0, dirac index %d in %g sec\n",id,part_time);
      
      for(int imass=0;imass<nmass;imass++)
	{ //reconstruct the doublet
	  reconstruct_doublet(reco_solution[0],reco_solution[1],cgmms_solution[imass],conf,kappa,mass[imass]);
	  if(rank==0) printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	  for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
	    for(int i=0;i<loc_vol;i++)
	      put_spincolor_into_colorspinspin(S0[r][imass][i],reco_solution[r][i],id);
	}
    }
  
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
      rotate_vol_colorspinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
}

//Calculate and print to file the 2pts
void calculate_all_2pts(int sm_lev_sour)
{
  //loop over smearing of the sink
  for(int sm_lev_sink=0;sm_lev_sink<nsm_lev;sm_lev_sink++)
    {
      char path[1024];
      sprintf(path,"%s_%d%d",outfile_2pts,sm_lev_sour,sm_lev_sink);
      FILE *fout=open_text_file_for_output(path);
      
      //in the case smear the sink
      if(sm_lev_sink!=0)
	for(int r=0;r<2;r++)
	  for(int imass=0;imass<nmass;imass++)
	    for(int id=0;id<4;id++)
	      {
		for(int ivol=0;ivol<loc_vol;ivol++) get_spincolor_from_colorspinspin(reco_solution[0][ivol],S0[r][imass][ivol],id);
		jacobi_smearing(reco_solution[1],reco_solution[0],sme_conf,jacobi_kappa,jacobi_niter[sm_lev_sink]-jacobi_niter[sm_lev_sink-1]);
		for(int ivol=0;ivol<loc_vol;ivol++) put_spincolor_into_colorspinspin(S0[r][imass][ivol],reco_solution[1][ivol],id);
	      }
      
      //perform the contractions
      contr_time-=take_time();      
      for(int im2=0;im2<nmass;im2++)
	for(int r2=0;r2<2;r2++)
	  for(int im1=0;im1<nmass;im1++)
	    for(int r1=0;r1<2;r1++)
	      {
		if(rank==0)
		  {
		    fprintf(fout," # m1=%f r1=%d , m2=%f r2=%d ,",mass[im1],r1,mass[im2],r2);
		    fprintf(fout," smear_source=%d smear_sink=%d\n",jacobi_niter[sm_lev_sour],jacobi_niter[sm_lev_sink]);
		  }
		
		meson_two_points(contr_2pts,op1_2pts,S0[r1][im1],op2_2pts,S0[r2][im2],ncontr_2pts);
		ncontr_tot+=ncontr_2pts;
		print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,twall,"");
		
		if(rank==0) fprintf(fout,"\n");
	      }
      contr_time+=take_time();
      
      if(rank==0) fclose(fout);
    }
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_appretto();
  
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  tot_time-=take_time();
  initialize_semileptonic(arg[1]);
  
  for(int iconf=0;iconf<ngauge_conf;iconf++)
  {
      //Gauge path
      read_str(conf_path,1024);
      read_int(&twall);
      read_str(outfile_2pts,1024);
      
      load_gauge_conf();
      generate_source();
      
      for(int sm_lev_sour=0;sm_lev_sour<nsm_lev;sm_lev_sour++)
      {
	  calculate_S0(sm_lev_sour);
	  calculate_all_2pts(sm_lev_sour);
      }
  }
  
  tot_time+=take_time();
  close_semileptonic();
  
  return 0;
}
