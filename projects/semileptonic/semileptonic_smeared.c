#include "nissa.h"

/*

  Program to compute two and three points correlation functions 
  relevant for semileptonic decays.
  
  Nomenclature
  
  |       *        |          Q2          |
  |   Q1 / \ Q2    |      .--------.      |
  |     /   \      |  SO *          * SI  |
  |  SO*_____*SE   |      '--------'      |
  |      spec      |          Q1          |
  
  SO = source
  SI = sink
  SE = sequential slice
  
  TSep specifiy the separation between sequential slice and source.
  
  For the SO, SE and SI a different number of smearing level can be specified
  Take into account that or each SO and SE level of smearing, an inversion is required,
  while different SI levels comes almost for free.
  
  Spectator for three points are specified in the input, according to
   iThetaMassR [iTheta] [iMass] [r]
  where iTheta and iMass are the indexes of the theta and mass inside their list
  and r=0,1 is the flavor in the twisted doublet.
  
  Two points are computed for each combination of masses and r of Q1 and Q2,
  and for each theta of Q2, but only for the values of theta Q1 listed in the "spec" list.
  
  Three points are computed for all the mass and theta values of Q1 and Q2, only
  for the charged twisted mass combination.
  
  To avoid the calculation of three points and at the same time the inversion, specify:
   NContrThreePoints 0
   NChContrThreePoints 0
   
  The list of configurations must be specified in the format:
   [confname] [tsource] [outfolder]
   
  The output folder must *not* be present and will be created by the program.
  If the directory is present, the configuration will not be analyzed.
  
  The program will run up to finish the list of configurations, or when 
  the WallTime (in seconds) is passed.
  
*/

#ifdef POINT_SOURCE_VERSION
#define prop_type su3spinspin
#else
#define prop_type colorspinspin
#endif

//gauge info
int ngauge_conf,nanalyzed_conf=0;
char conf_path[1024],outfolder[1024];
quad_su3 *conf,*sme_conf;
as2t_su3 *Pmunu;
double kappa;
double put_theta[4],old_theta[4]={0,0,0,0};

//list of masses and theta
int nmassS0,nthetaS0;
int nmassS1,nthetaS1;
double *massS0,*thetaS0;
double *massS1,*thetaS1;

//source data
int seed,noise_type;
coords source_coord;
spincolor *source;
prop_type *original_source;

//smearing parameters
double jacobi_kappa,ape_alpha;
int ape_niter;
int *jacobi_niter_so,nsm_lev_so;
int *jacobi_niter_si,nsm_lev_si;
int *jacobi_niter_se,nsm_lev_se;

//vectors for the spinor data
int npropS0;
prop_type **S0[2];
int ncgmms_solution;
spincolor **cgmms_solution,*temp_vec[2];
int ithetaS0_min,ithetaS0_max,nthetaS0;

//cgmms inverter parameters
double stopping_residue;
double minimal_residue;
int stopping_criterion;
int niter_max;

//two points contractions
int ncontr_2pts;
complex *contr_2pts;
int *op1_2pts,*op2_2pts;

//two points chromo contractions
int nch_contr_2pts;
complex *ch_contr_2pts;
int *ch_op1_2pts,*ch_op2_2pts;
prop_type *ch_prop;

//sequential props info
int nspec;
int tsep;
int *imass_spec,*r_spec,*ith_spec;
prop_type *sequential_source;

//sequential propagators
int npropS1;
prop_type **S1;

//three points contractions
int ncontr_3pts;
complex *contr_3pts;
int *op1_3pts,*op2_3pts;

//two points chromo contractions
int nch_contr_3pts;
complex *ch_contr_3pts;
int *ch_op1_3pts,*ch_op2_3pts;

//timings
int ninv_tot=0,ncontr_tot=0;
int wall_time;
double tot_time=0,inv_time=0;
double smear_time=0;
double load_time=0,contr_time;
double contr_save_time=0;
double contr_2pts_time=0;
double contr_3pts_time=0;

//return the position of the propagator of theta and mass
int ipropS0(int itheta,int imass){return itheta*nmassS0+imass;}
int ipropS1(int itheta,int imass){return itheta*nmassS1+imass;}

//This function takes care to make the revert on the FIRST spinor, putting the needed gamma5
void meson_two_points(complex *corr,int *list_op1,prop_type *s1,int *list_op2,prop_type *s2,int ncontr)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr],t2[ncontr];
  
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Put the two gamma5 needed for the revert of the first spinor
      dirac_prod(&(t1[icontr]), &(base_gamma[list_op1[icontr]]),&(base_gamma[5]));
      dirac_prod(&(t2[icontr]), &(base_gamma[5]),&(base_gamma[list_op2[icontr]]));
    }
  
  //Call the routine which perform the contraction
#ifdef POINT_SOURCE_VERSION
  trace_g_ccss_dag_g_ccss(corr,t1,s1,t2,s2,ncontr);
#else
  trace_g_sdag_g_s(corr,t1,s1,t2,s2,ncontr);
#endif
}

//This function contract a source with a sequential spinor putting the passed list of operators
void contract_with_source(complex *corr,prop_type *S1,int *list_op,prop_type *source)
{
  //Temporary vector for the internal matrices
  dirac_matr t1[ncontr_2pts],t2[ncontr_2pts];

  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Init the second 
      t1[icontr]=base_gamma[0];
      t2[icontr]=base_gamma[list_op[icontr]];
    }

  //Call the routine which does the real contraction
#ifdef POINT_SOURCE_VERSION
  trace_g_ccss_dag_g_ccss(corr,t1,S1,t2,source,ncontr_2pts);
#else
  trace_g_sdag_g_s(corr,t1,S1,t2,source,ncontr_2pts);
#endif
}

//generate the source
void generate_source()
{
#ifdef POINT_SOURCE_VERSION
  generate_delta_source(original_source,source_coord);
#else
  enum rnd_type type[5]={RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_Z2,RND_Z2,RND_Z4};
  generate_spindiluted_source(original_source,type[noise_type],source_coord[0]);
#endif
}

//Generate a sequential source for S1
void generate_sequential_source(int ispec)
{
  int r=r_spec[ispec];
  
  master_printf("\nCreating the sequential source for spectator %d\n",ispec);
  nissa_loc_vol_loop(ivol)
    {
      //put to zero everything but the slice
      if(glb_coord_of_loclx[ivol][0]!=(source_coord[0]+tsep)%glb_size[0])
	memset(sequential_source[ivol],0,sizeof(prop_type));
      else
	{
	  //avoid to put g5, beacuse commutate with (i+-g5)/sqrt(2) and cancel with those of the QQ
	  memcpy(sequential_source[ivol],S0[r][ipropS0(ith_spec[ispec],imass_spec[ispec])][ivol],sizeof(prop_type));
	  for(int c=0;c<3;c++) //rotate as r because it's D^-1
	    {
#ifdef POINT_SOURCE_VERSION
	      for(int c1=0;c1<3;c1++) rotate_spinspin_to_physical_basis(sequential_source[ivol][c][c1],r,r);
#else
	      rotate_spinspin_to_physical_basis(sequential_source[ivol][c],r,r);
#endif
	    }
	}
    }
  master_printf("Sequential source created\n\n");
}  

//Parse all the input file
void initialize_semileptonic(char *input_path)
{
  open_input(input_path);

  // 1) Read information about the gauge conf
  
  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Wall time
  read_str_int("WallTime",&wall_time);
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the source
  
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  //Read the noise type
  read_str_int("NoiseType",&noise_type);

  // 3) Smearing parameters
  
  //Smearing parameters
  read_str_double("ApeAlpha",&ape_alpha);
  read_str_int("ApeNiter",&ape_niter);
  read_str_double("JacobiKappa",&jacobi_kappa);
  read_list_of_ints("JacobiNiterSo",&nsm_lev_so,&jacobi_niter_so);
  for(int iter=1;iter<nsm_lev_so;iter++)
    if(jacobi_niter_so[iter]<jacobi_niter_so[iter-1])
      crash("Error, jacobi lev sou %d minor than %d (%d, %d)!\n",iter,iter-1,jacobi_niter_so[iter],jacobi_niter_so[iter-1]);
  read_list_of_ints("JacobiNiterSe",&nsm_lev_se,&jacobi_niter_se);
  for(int iter=1;iter<nsm_lev_se;iter++)
    if(jacobi_niter_se[iter]<jacobi_niter_se[iter-1])
      crash("Error, jacobi lev seq %d minor than %d (%d, %d)!\n",iter,iter-1,jacobi_niter_se[iter],jacobi_niter_se[iter-1]);
  read_list_of_ints("JacobiNiterSi",&nsm_lev_si,&jacobi_niter_si);
  for(int iter=1;iter<nsm_lev_si;iter++)
    if(jacobi_niter_si[iter]<jacobi_niter_si[iter-1])
      crash("Error, jacobi lev seq %d minor than %d (%d, %d)!\n",iter,iter-1,jacobi_niter_si[iter],jacobi_niter_si[iter-1]);
  
  // 4) Read list of masses and of thetas

  read_list_of_doubles("NMassS0",&nmassS0,&massS0);
  read_list_of_doubles("NThetaS0",&nthetaS0,&thetaS0);
  
  // 5) Info about the inverter
  
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
  
  if(stopping_criterion==numb_known_stopping_criterion)
    {
      master_fprintf(stderr,"Unknown stopping criterion: %s\n",str_stopping_criterion);
      master_fprintf(stderr,"List of known stopping criterions:\n");
      for(int isc=0;isc<numb_known_stopping_criterion;isc++) master_fprintf(stderr," %s\n",list_known_stopping_criterion[isc]);
      crash("check the input file");
    }
  if(stopping_criterion==sc_standard) read_str_double("MinimalResidue",&minimal_residue);
      
  //Number of iterations
  read_str_int("NiterMax",&niter_max);
  
  // 6) contraction list for two points
  
  read_str_int("NContrTwoPoints",&ncontr_2pts);
  contr_2pts=nissa_malloc("contr_2pts",ncontr_2pts*glb_size[0],complex);
  op1_2pts=nissa_malloc("op1_2pts",ncontr_2pts,int);
  op2_2pts=nissa_malloc("op2_2pts",ncontr_2pts,int);
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1_2pts[icontr]));
      read_int(&(op2_2pts[icontr]));
      
      master_printf(" contr.%d %d %d\n",icontr,op1_2pts[icontr],op2_2pts[icontr]);
    }
  
  read_str_int("NChromoContrTwoPoints",&nch_contr_2pts);
  ch_contr_2pts=nissa_malloc("ch_contr_2pts",nch_contr_2pts*glb_size[0],complex);
  ch_op1_2pts=nissa_malloc("ch_op1_2pts",ncontr_2pts,int);
  ch_op2_2pts=nissa_malloc("ch_op2_2pts",ncontr_2pts,int);
  for(int icontr=0;icontr<nch_contr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(ch_op1_2pts[icontr]));
      read_int(&(ch_op2_2pts[icontr]));
      
      master_printf(" ch-contr.%d %d %d\n",icontr,ch_op1_2pts[icontr],ch_op2_2pts[icontr]);
    }
  
  // 7) Read list of masses and of thetas for S1
  
  read_list_of_doubles("NMassS1",&nmassS1,&massS1);
  read_list_of_doubles("NThetaS1",&nthetaS1,&thetaS1);
  
  // 8) three points functions
  
  sequential_source=nissa_malloc("Sequential source",loc_vol,prop_type);
  read_str_int("TSep",&tsep);
  read_str_int("NSpec",&nspec);
  if(nspec==0) crash("it has no meaning to specify 0 spectators");
  ith_spec=nissa_malloc("ith_spec",nspec,int);
  imass_spec=nissa_malloc("imass_spec",nspec,int);
  r_spec=nissa_malloc("r_spec",nspec,int);
  for(int ispec=0;ispec<nspec;ispec++)
    {
      expect_str("iThetaMassR");
      read_int(&(ith_spec[ispec]));
      read_int(&(imass_spec[ispec]));
      read_int(&(r_spec[ispec]));
      
      if(ith_spec[ispec]<0||ith_spec[ispec]>nthetaS0)    crash("theta for ispec %d out of bounds",ispec);
      if(imass_spec[ispec]<0||imass_spec[ispec]>nmassS0) crash("mass for ispec %d out of bounds",ispec);
      if(r_spec[ispec]<0||r_spec[ispec]>1)               crash("r for ispec %d out of bounds",ispec);
      
      master_printf(" spec %d: th=%g, m=%g, r=%d\n",ispec,thetaS0[ith_spec[ispec]],massS0[imass_spec[ispec]],r_spec[ispec]);
    }
  read_str_int("NContrThreePoints",&ncontr_3pts);
  contr_3pts=nissa_malloc("contr_3pts",ncontr_3pts*glb_size[0],complex); 
  op1_3pts=nissa_malloc("op1_3pts",ncontr_3pts,int);
  op2_3pts=nissa_malloc("op2_3pts",ncontr_3pts,int);
  for(int icontr=0;icontr<ncontr_3pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1_3pts[icontr]));
      read_int(&(op2_3pts[icontr]));
      
      master_printf(" contr.%d %d %d\n",icontr,op1_3pts[icontr],op2_3pts[icontr]);
    }
  
  read_str_int("NChromoContrThreePoints",&nch_contr_3pts);
  ch_contr_3pts=nissa_malloc("ch_contr_3pts",nch_contr_3pts*glb_size[0],complex);
  ch_op1_3pts=nissa_malloc("ch_op1_3pts",nch_contr_3pts,int);
  ch_op2_3pts=nissa_malloc("ch_op2_3pts",nch_contr_3pts,int);
  for(int icontr=0;icontr<nch_contr_3pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(ch_op1_3pts[icontr]));
      read_int(&(ch_op2_3pts[icontr]));
      
      master_printf(" ch-contr.%d %d %d\n",icontr,ch_op1_3pts[icontr],ch_op2_3pts[icontr]);
    }
  
  read_str_int("NGaugeConf",&ngauge_conf);
  
  master_printf("\n");
  
  ////////////////////////////////////// end of input reading/////////////////////////////////
  
  //allocate gauge conf, Pmunu and all the needed spincolor and propagators
  conf=nissa_malloc("or_conf",loc_vol+loc_bord+loc_edge,quad_su3);
  sme_conf=nissa_malloc("sm_conf",loc_vol+loc_bord,quad_su3);
  Pmunu=nissa_malloc("Pmunu",loc_vol,as2t_su3);
  
  //Allocate all the S0 prop_type vectors
  npropS0=nthetaS0*nmassS0;
  S0[0]=nissa_malloc("S0[0]",npropS0,prop_type*);
  S0[1]=nissa_malloc("S0[1]",npropS0,prop_type*);
  for(int iprop=0;iprop<npropS0;iprop++)
    {
      S0[0][iprop]=nissa_malloc("S0[0]",loc_vol,prop_type);
      S0[1][iprop]=nissa_malloc("S0[1]",loc_vol,prop_type);
    }
  
  //Allocate nmass spincolors, for the cgmms solutions
  ncgmms_solution=max_int(nmassS0,nmassS1);
  cgmms_solution=nissa_malloc("cgmms_solution",ncgmms_solution,spincolor*);
  for(int imass=0;imass<ncgmms_solution;imass++) cgmms_solution[imass]=nissa_malloc("cgmms_solution",loc_vol+loc_bord,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol,spincolor);
  
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+loc_bord,spincolor);
  original_source=nissa_malloc("original_source",loc_vol,prop_type);
  
  //Allocate one prop_type for the chromo-contractions
  ch_prop=nissa_malloc("chromo-prop",loc_vol,prop_type);
  
  //Allocate all the S1 prop_type vectors
  npropS1=nthetaS1*nmassS1;
  S1=nissa_malloc("S1",loc_vol,prop_type*);
  for(int iprop=0;iprop<npropS0;iprop++) S1[iprop]=nissa_malloc("S1[i]",loc_vol,prop_type);
}

//find a new conf
int read_conf_parameters(int *iconf)
{
  int ok_conf;

  do
    {
      //Gauge path
      read_str(conf_path,1024);
      
      //Source coord
      read_int(&(source_coord[0]));
#ifdef POINT_SOURCE_VERSION
      read_int(&(source_coord[1]));
      read_int(&(source_coord[2]));
      read_int(&(source_coord[3]));
#endif
      
      //Out folder
      read_str(outfolder,1024);
      
      //Check if the conf exist
      master_printf("Considering configuration \"%s\" with output path \"%s\".\n",outfolder,conf_path);
      ok_conf=!(dir_exists(outfolder));
      if(ok_conf)
	{
	  int ris=create_dir(outfolder);
	  if(ris==0) master_printf(" Output path \"%s\" not present: configuration \"%s\" not yet analyzed, starting.\n",outfolder,conf_path);
	  else
	    {
	      ok_conf=0;
	      master_printf(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
	    }
	}
      else
	master_printf(" Output path \"%s\" already present: configuration \"%s\" already analyzed, skipping.\n",outfolder,conf_path);
      (*iconf)++;
    }
  while(!ok_conf && (*iconf)<ngauge_conf);
  
  master_printf("\n");
  
  return ok_conf;
}

//read the conf and setup it
void setup_conf()
{
  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_ildg_gauge_conf(conf,conf_path);
  Pmunu_term(Pmunu,conf);
  
  //prepare the smerded version
  ape_spatial_smear_conf(sme_conf,conf,ape_alpha,ape_niter);
  
  //compute plaquette
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  master_printf("smerded plaq: %.18g\n",global_plaquette_lx_conf(sme_conf));
  
  //put the anti-periodic condition on the temporal border
  old_theta[0]=old_theta[1]=old_theta[2]=old_theta[3]=0;
  put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
  adapt_theta(conf,old_theta,put_theta,1,1);
}

//Finalization
void close_semileptonic()
{
  contr_time=contr_2pts_time+contr_3pts_time+contr_save_time;
  
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g s, of which:\n",tot_time);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf(" - %02.2f%s to sink-smear propagators\n",smear_time*100.0/tot_time,"%");
  master_printf(" - %02.2f%s to perform %d contr. (%2.2gs avg) of which:\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
  master_printf("   * %02.2f%s to compute two points\n",contr_2pts_time*100.0/contr_time,"%");
  master_printf("   * %02.2f%s to compute three points\n",contr_3pts_time*100.0/contr_time,"%");
  master_printf("   * %02.2f%s to save correlations\n",contr_save_time*100.0/contr_time,"%");
  
  nissa_free(Pmunu);nissa_free(conf);nissa_free(sme_conf);
  for(int iprop=0;iprop<npropS0;iprop++){nissa_free(S0[0][iprop]);nissa_free(S0[1][iprop]);nissa_free(S1[iprop]);}
  nissa_free(S0[0]);nissa_free(S0[1]);nissa_free(S1);
  nissa_free(temp_vec[0]);nissa_free(temp_vec[1]);
  nissa_free(ch_prop);nissa_free(sequential_source);
  nissa_free(contr_2pts);nissa_free(ch_contr_2pts);
  nissa_free(contr_3pts);nissa_free(ch_contr_3pts);
  nissa_free(op1_2pts);nissa_free(op2_2pts);
  nissa_free(op1_3pts);nissa_free(op2_3pts);
  nissa_free(ch_op1_2pts);nissa_free(ch_op2_2pts);
  nissa_free(ch_op1_3pts);nissa_free(ch_op2_3pts);
  nissa_free(ith_spec);nissa_free(r_spec);nissa_free(imass_spec);
  for(int imass=0;imass<ncgmms_solution;imass++) nissa_free(cgmms_solution[imass]);
  nissa_free(cgmms_solution);
  nissa_free(source);nissa_free(original_source);
  close_nissa();
}

//smear addditivily a propagator
void smear_additive_propagator(prop_type *out,prop_type *in,int ism_lev,int *jacobi_niter)
{
  spincolor *temp=nissa_malloc("temp",loc_vol+loc_bord,spincolor);
  
  int nsme=jacobi_niter[ism_lev];
  if(ism_lev>0) nsme-=jacobi_niter[ism_lev-1];
  
  //loop over dirac index
#ifdef POINT_SOURCE_VERSION
  for(int ic=0;ic<3;ic++)
#endif
    for(int id=0;id<4;id++)
      {
	nissa_loc_vol_loop(ivol)
	  {
#ifdef POINT_SOURCE_VERSION
	    get_spincolor_from_su3spinspin(temp[ivol],in[ivol],id,ic);
#else
	    get_spincolor_from_colorspinspin(temp[ivol],in[ivol],id);
#endif
	  }
	set_borders_invalid(temp);
	
	jacobi_smearing(temp,temp,sme_conf,jacobi_kappa,nsme);
	
	nissa_loc_vol_loop(ivol)
	  {
#ifdef POINT_SOURCE_VERSION
	    put_spincolor_into_su3spinspin(out[ivol],temp[ivol],id,ic);
#else
	    put_spincolor_into_colorspinspin(out[ivol],temp[ivol],id);
#endif
	  }
      }
  
  nissa_free(temp);
}

//calculate the standard propagators
void calculate_S0(int ism_lev_so)
{
  //smear additively the source
  master_printf("\nSource Smearing level: %d\n",ism_lev_so);
  smear_additive_propagator(original_source,original_source,ism_lev_so,jacobi_niter_so);
  master_printf("\n");
  
 //loop over the source dirac index
#ifdef POINT_SOURCE_VERSION
  for(int ic=0;ic<3;ic++)
#endif
    for(int id=0;id<4;id++)
      { 
	//put the g5
	nissa_loc_vol_loop(ivol)
	  {
#ifdef POINT_SOURCE_VERSION
	    get_spincolor_from_su3spinspin(source[ivol],original_source[ivol],id,ic);
#else
	    get_spincolor_from_colorspinspin(source[ivol],original_source[ivol],id);
#endif
	    safe_dirac_prod_spincolor(source[ivol],&(base_gamma[5]),source[ivol]);
	  }
	set_borders_invalid(source);
	
	for(int itheta=0;itheta<nthetaS0;itheta++)
	  {
	    //adapt the border condition
	    put_theta[1]=put_theta[2]=put_theta[3]=thetaS0[itheta];
	    adapt_theta(conf,old_theta,put_theta,1,1);
	    
	    //inverting
	    double part_time=-take_time();
	    inv_tmQ2_cgmms(cgmms_solution,conf,kappa,massS0,nmassS0,niter_max,stopping_residue,minimal_residue,stopping_criterion,source);
	    part_time+=take_time();ninv_tot++;inv_time+=part_time;
#ifdef POINT_SOURCE_VERSION
	    master_printf("Finished the inversion of S0 theta %d, color index %d, dirac index %d in %g sec\n",itheta,ic,id,part_time);
#else
	    master_printf("Finished the inversion of S0 theta %d, dirac index %d in %g sec\n",itheta,id,part_time);
#endif	    
	    //reconstruct the doublet
	    for(int imass=0;imass<nmassS0;imass++)
	      {
		reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,massS0[imass],cgmms_solution[imass]);
		master_printf("Mass %d (%g) reconstructed \n",imass,massS0[imass]);
		for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
		  nissa_loc_vol_loop(ivol)
		    {
#ifdef POINT_SOURCE_VERSION
		      put_spincolor_into_su3spinspin(S0[r][ipropS0(itheta,imass)][ivol],temp_vec[r][ivol],id,ic);
#else
		      put_spincolor_into_colorspinspin(S0[r][ipropS0(itheta,imass)][ivol],temp_vec[r][ivol],id);
#endif
		    }
	      }
	  }
      }
  
  //rotate to physical basis
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
      {
#ifdef POINT_SOURCE_VERSION
	rotate_vol_su3spinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
#else
	rotate_vol_colorspinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
#endif
      }
  master_printf("Propagators rotated\n");

  master_printf("\n");
}

//calculate the sequential propagators
void calculate_S1(int ispec,int ism_lev_se)
{
  //smear additively the seq
  master_printf("\nSeq Smearing level: %d (will be applied twice!)\n",ism_lev_se);
  smear_additive_propagator(sequential_source,sequential_source,ism_lev_se,jacobi_niter_se);
  smear_additive_propagator(sequential_source,sequential_source,ism_lev_se,jacobi_niter_se);
  master_printf("\n");
  
  //loop over se
  for(int ic=0;ic<3;ic++)
    for(int id=0;id<4;id++)
      { 
	nissa_loc_vol_loop(ivol) //avoid the g5 insertion
	  {
#ifdef POINT_SOURCE_VERSION
	    get_spincolor_from_su3spinspin(source[ivol],sequential_source[ivol],id,ic);
#else
	    get_spincolor_from_colorspinspin(source[ivol],sequential_source[ivol],id);
#endif
	  }
	set_borders_invalid(source);
	
	for(int itheta=0;itheta<nthetaS1;itheta++)
	  {
	    //adapt the border condition
	    put_theta[1]=put_theta[2]=put_theta[3]=thetaS1[itheta];
	    adapt_theta(conf,old_theta,put_theta,1,1);
	    
	    double part_time=-take_time();
	    inv_tmQ2_cgmms(cgmms_solution,conf,kappa,massS1,nmassS1,niter_max,stopping_residue,minimal_residue,stopping_criterion,source);
	    part_time+=take_time();ninv_tot++;inv_time+=part_time;
	    master_printf("Finished the inversion of S1 theta %d, seq sme lev %d, dirac index %d in %g sec\n",itheta,ism_lev_se,id,part_time);
	    
	    for(int imass=0;imass<nmassS1;imass++)
	      {
		//reconstruct the doublet: r(S1)=!r(spec), so we have to multiply by Q+ if r(spec)==1 and Q- if 0
		double reco_mass=-massS1[imass];
		if(r_spec[ispec]==1) reco_mass=-reco_mass;
		//use temp_vec[0] as temporary storage
		apply_tmQ(temp_vec[0],conf,kappa,reco_mass,cgmms_solution[imass]);
		master_printf("Mass %d (%g) reconstructed \n",imass,massS1[imass]);
		nissa_loc_vol_loop(ivol)
		  {
#ifdef POINT_SOURCE_VERSION
		    put_spincolor_into_su3spinspin(S1[ipropS1(itheta,imass)][ivol],temp_vec[0][ivol],id,ic);
#else
		    put_spincolor_into_colorspinspin(S1[ipropS1(itheta,imass)][ivol],temp_vec[0][ivol],id);
#endif
		  }
	      }
	  }
      }

  //put the (1+-ig5)/sqrt(2) factor. On the source rotate as r_spec, on the sink as !r_spec
  for(int ipropS1=0;ipropS1<npropS1;ipropS1++) //but, being D^-1, everything is swapped
    {
#ifdef POINT_SOURCE_VERSION
      rotate_vol_su3spinspin_to_physical_basis(S1[ipropS1],!(r_spec[ispec]),(r_spec[ispec]));
#else
      rotate_vol_colorspinspin_to_physical_basis(S1[ipropS1],!(r_spec[ispec]),(r_spec[ispec]));
#endif
    }

  master_printf("Propagators rotated\n");

  master_printf("\n");
}

//Calculate and print to file the 2pts
void calculate_all_2pts(int ism_lev_so,int ism_lev_si)
{
  smear_time-=take_time();
    
  for(int r=0;r<2;r++)
    for(int iprop=0;iprop<npropS0;iprop++)
      smear_additive_propagator(S0[r][iprop],S0[r][iprop],ism_lev_si,jacobi_niter_si);
  
  double temp_time=take_time();
  smear_time+=temp_time;
  contr_2pts_time-=temp_time;
  
  char path[1024];
  sprintf(path,"%s/2pts_%02d_%02d",outfolder,jacobi_niter_so[ism_lev_so],jacobi_niter_si[ism_lev_si]);
  FILE *fout=open_text_file_for_output(path);
  
  for(int ispec=0;ispec<nspec;ispec++)
    {
      int ith1=ith_spec[ispec];
      for(int ith2=0;ith2<nthetaS0;ith2++)
	for(int im2=0;im2<nmassS0;im2++)
  	  {
	    int ip2=ipropS0(ith2,im2);
	    for(int r2=0;r2<2;r2++)
	      {
		if(nch_contr_2pts>0)
		  {
#ifdef POINT_SOURCE_VERSION
		    unsafe_apply_chromo_operator_to_su3spinspin(ch_prop,Pmunu,S0[r2][ip2]);
#else
		    unsafe_apply_chromo_operator_to_colorspinspin(ch_prop,Pmunu,S0[r2][ip2]);
#endif
		  }
		
		for(int im1=0;im1<nmassS0;im1++)
		  {
		    int ip1=ipropS0(ith1,im1);
		    
		    for(int r1=0;r1<2;r1++)
		      {
			//header
			master_fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d",massS0[im1],thetaS0[ith1],r1,massS1[im2],thetaS1[ith2],r2);
			master_fprintf(fout," smear_source=%d smear_sink=%d\n",jacobi_niter_so[ism_lev_so],jacobi_niter_si[ism_lev_si]);
			
			//compute contractions
			meson_two_points(contr_2pts,op1_2pts,S0[r1][ip1],op2_2pts,S0[r2][ip2],ncontr_2pts);
			ncontr_tot+=ncontr_2pts;
			
			//write 
			contr_save_time-=take_time();
			print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,source_coord[0],"",1.0);
			contr_save_time+=take_time();
			
			//if chromo contractions
			if(nch_contr_2pts>0)
			  {
			    //compute them
			    meson_two_points(ch_contr_2pts,ch_op1_2pts,S0[r1][ip1],ch_op2_2pts,ch_prop,nch_contr_2pts);
			    ncontr_tot+=nch_contr_2pts;
			    
			    //print them
			    contr_save_time-=take_time();
			    print_contractions_to_file(fout,nch_contr_2pts,ch_op1_2pts,ch_op2_2pts,ch_contr_2pts,source_coord[0],"CHROMO-",1.0);
			    contr_save_time+=take_time();
			  }
			master_fprintf(fout,"\n");
		      }
		    
		    ncontr_tot+=nch_contr_2pts;
		  }
	      }
	  }
    }
  
  contr_2pts_time+=take_time();
  if(rank==0) fclose(fout);
}

//Calculate and print to file the 3pts
void calculate_all_3pts(int ispec,int ism_lev_so,int ism_lev_se)
{
  char path[1024];
  
  sprintf(path,"%s/3pts_sp%d_%02d_%02d",outfolder,ispec,jacobi_niter_so[ism_lev_so],jacobi_niter_se[ism_lev_se]);
  
  FILE *fout=open_text_file_for_output(path);
  
  contr_3pts_time-=take_time();
  
  int r1=r_spec[ispec];
  int r2=!r1;
  for(int ith2=0;ith2<nthetaS1;ith2++)
    for(int im2=0;im2<nmassS1;im2++)
      {
	int ip2=ipropS1(ith2,im2);
	if(nch_contr_3pts>0)
	  {
#ifdef POINT_SOURCE_VERSION
	    unsafe_apply_chromo_operator_to_su3spinspin(ch_prop,Pmunu,S1[ip2]);
#else
	    unsafe_apply_chromo_operator_to_colorspinspin(ch_prop,Pmunu,S1[ip2]);
#endif
	  }
	
	for(int ith1=0;ith1<nthetaS0;ith1++)
	    for(int im1=0;im1<nmassS0;im1++)
	      {
		int ip1=ipropS0(ith1,im1);
		
		//header
		master_fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d,",massS0[im1],thetaS0[ith1],r1,massS1[im2],thetaS1[ith2],r2);
		master_fprintf(fout," smear_source=%d smear_seq=%d\n",jacobi_niter_so[ism_lev_so],jacobi_niter_se[ism_lev_se]);
		
		//compute contractions
		meson_two_points(contr_3pts,op1_3pts,S0[r1][ip1],op2_3pts,S1[ip2],ncontr_3pts);
		ncontr_tot+=ncontr_3pts;
		
		//write them
		contr_save_time-=take_time();
		print_contractions_to_file(fout,ncontr_3pts,op1_3pts,op2_3pts,contr_3pts,source_coord[0],"",1.0);
		contr_save_time+=take_time();
		
		//if chromo contractions
		if(nch_contr_3pts>0)
		  {
		    //compute them
		    meson_two_points(ch_contr_3pts,ch_op1_3pts,S0[r1][ip1],ch_op2_3pts,ch_prop,nch_contr_3pts);
		    ncontr_tot+=nch_contr_3pts;
		    
		    //and write them
		    contr_save_time-=take_time();
		    print_contractions_to_file(fout,nch_contr_3pts,ch_op1_3pts,ch_op2_3pts,ch_contr_3pts,source_coord[0],"CHROMO-",1.0);
		    contr_save_time+=take_time();
		  }
		
		master_fprintf(fout,"\n");
	    }
    }
  
  contr_3pts_time+=take_time();
  
  if(rank==0) fclose(fout);
}

//check all the two points
void check_two_points(int ispec,int ism_lev_so,int ism_lev_se)
{
  char path[1024];
  sprintf(path,"%s/2pts_check_sp%d_%02d_%02d",outfolder,ispec,jacobi_niter_so[ism_lev_so],jacobi_niter_se[ism_lev_se]);
  FILE *fout=open_text_file_for_output(path);

  for(int ith2=0;ith2<nthetaS1;ith2++)
    for(int im2=0;im2<nmassS1;im2++)
      {
	int ip2=ipropS1(ith2,im2);
	contract_with_source(contr_2pts,S1[ip2],op2_2pts,original_source);
	
	master_fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",
		       massS0[imass_spec[ispec]],thetaS0[ith_spec[ispec]],r_spec[ispec],massS1[im2],thetaS1[ith2],r_spec[ispec]);
	
	master_fprintf(fout,"\n");
	
	for(int icontr=0;icontr<ncontr_2pts;icontr++)
	  if(op1_2pts[icontr]==5)
	    master_fprintf(fout," # P5%s\t%+016.16g\t%+016.16g\n",gtag[op2_2pts[icontr]],(contr_2pts+icontr*glb_size[0])[source_coord[0]][0],(contr_2pts+icontr*glb_size[0])[source_coord[0]][1]);
	master_fprintf(fout,"\n");
      }
  
  if(rank==0) fclose(fout);
}

//check if the time is enough
int check_remaining_time()
{
  int enough_time;

  //check remaining time                                                                                                                                                                        
  double temp_time=take_time()+tot_time;
  double ave_time=temp_time/nanalyzed_conf;
  double left_time=wall_time-temp_time;
  enough_time=left_time>(ave_time*1.1);

  master_printf("Remaining time: %lg sec\n",left_time);
  master_printf("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
  if(enough_time) master_printf("Continuing with next conf!\n");
  else master_printf("Not enough time, exiting!\n");
  
  return enough_time;
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();
  tot_time-=take_time();
  
  //initialize the program
  if(narg<2) crash("Use: %s input_file",arg[0]);
  initialize_semileptonic(arg[1]);
  
  //loop over the configs
  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf && enough_time && read_conf_parameters(&iconf))
    {
      //smear the conf and generate the source
      setup_conf();
      generate_source();
      
      //loop on smearing of the source
      for(int sm_lev_so=0;sm_lev_so<nsm_lev_so;sm_lev_so++)
	{
	  //compute S0 propagator smearing the config the appropriate number of time the source
	  calculate_S0(sm_lev_so);
	  
	  //loop on spectator
	  if(ncontr_3pts!=0 || nch_contr_3pts!=0)
	    for(int ispec=0;ispec<nspec;ispec++)
	      {
		//select a timeslice and multiply by gamma5
		generate_sequential_source(ispec);
		
		//loop on smearing of the sequential prop
		for(int sm_lev_se=0;sm_lev_se<nsm_lev_se;sm_lev_se++)
		  {
		    calculate_S1(ispec,sm_lev_se);
		    check_two_points(ispec,sm_lev_so,sm_lev_se);
		    calculate_all_3pts(ispec,sm_lev_so,sm_lev_se);
		  }
	      }
	  
	  //loop on the smearing of the sink
	  for(int sm_lev_si=0;sm_lev_si<nsm_lev_si;sm_lev_si++) calculate_all_2pts(sm_lev_so,sm_lev_si);
	}
      
      //pass to the next conf if there is enough time
      nanalyzed_conf++;
      enough_time=check_remaining_time();
    }
  
  if(iconf==ngauge_conf) master_printf("Finished all the conf!\n");
  
  tot_time+=take_time();
  close_semileptonic();
  
  return 0;
}
