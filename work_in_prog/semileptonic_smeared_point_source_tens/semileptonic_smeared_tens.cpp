#include "nissa.h"

/*

  Program to compute two and three points correlation functions 
  relevant for semileptonic decays.
  
  Nomenclature
  
  |         *        |         Q2          |
  |     Q1 / \ Q2    |     .--------.      |
  |       /   \      | SO *          *  SI |
  |    SO*_____*SE   |     '--------'      |
  |        spec      |         Q1          |
  
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
  

//gauge info
int ngauge_conf,nanalized_conf=0;
char conf_path[1024];
quad_su3 *conf,*sme_conf;
as2t_su3 *Pmunu;
double kappa;
char outfolder[1024];

//list of masses and theta
int nmass,ntheta;
double *mass,*theta;
double put_theta[4],old_theta[4]={0,0,0,0};

//source data
int seed,noise_type,twall;
spincolor *source;
colorspinspin *original_source;
spincolor *der_source;

//smearing parameters
double jacobi_kappa,ape_alpha;
int ape_niter;
int *jacobi_niter_so,nsm_lev_so;
int *jacobi_niter_si,nsm_lev_si;
int *jacobi_niter_se,nsm_lev_se;

//vectors for the spinor data
int npropS0;
colorspinspin **S0[2];
spincolor **cgm_solution,*temp_vec[2];
colorspinspin **S0_der[2][3][3];

spincolor *temp_der_vec[2];
//spincolor **cgm_der_solution;

//cgm inverter parameters
double *stopping_residues;
int niter_max=100000;

//two points contractions
int ncontr_2pts;
complex *contr_2pts;
int *op1_2pts,*op2_2pts;

//two points chromo contractions
int nch_contr_2pts;
complex *ch_contr_2pts;
int *ch_op1_2pts,*ch_op2_2pts;
colorspinspin *ch_colorspinspin;

//sequential props info
int nspec;
int tsep;
int *imass_spec,*r_spec,*ith_spec;
colorspinspin *sequential_source;

//sequential propagators
int npropS1;
colorspinspin **S1;

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
double load_time=0,contr_time;
double contr_save_time=0;
double contr_2pts_time=0;
double contr_3pts_time=0;


//Apply the covariant derivative having direction i on a spincolor

void apply_nabla_i_spincolor(spincolor *out, spincolor *in, quad_su3 *conf, int i)
{
	
	memset(out,0,loc_vol*sizeof(spincolor));
	communicate_lx_spincolor_borders(in);
	
	for (int ix=0; ix< loc_vol; ix++) 
	{
		int Xup,Xdw;
		Xup=loclx_neighup[ix][i];
		Xdw=loclx_neighdw[ix][i];
		
		unsafe_su3_prod_spincolor(out[ix], conf[ix][i],in[Xup]);
		unsafe_su3_dag_subt_the_prod_spincolor(out[ix],conf[Xdw][i],in[Xdw]);
	}
	set_borders_invalid(out);
}

/*
void apply_nabla_dag_i_spincolor(spincolor *out, spincolor *in, quad_su3 *conf, int i)
{
	
	memset(out,0,loc_vol*sizeof(spincolor));
	communicate_lx_spincolor_borders(in);
	
	for (int ix=0; ix< loc_vol; ix++) 
	{
		int Xup,Xdw;
		Xup=loclx_neighup[ix][i];
		Xdw=loclx_neighdw[ix][i];
		
		safe_su3_dag_prod_spincolor(out[ix], conf[ix][i],in[Xup]);
		su3_subt_the_prod_spincolor(out[ix],conf[Xdw][i],in[Xdw]);
		
	}
	set_borders_invalid(out);
}
*/

//return the position of the propagator of theta and mass
int iprop_of(int itheta,int imass){return itheta*nmass+imass;}

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

//This function contract a source with a sequential spinor putting the passed list of operators
void contract_with_source(complex *corr,colorspinspin *S1,int *list_op,colorspinspin *source)
{
  //Temporary vector for the internal matrices
  dirac_matr t1[ncontr_2pts],t2[ncontr_2pts];

  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    { //Init the second 
      t1[icontr]=base_gamma[0];
      t2[icontr]=base_gamma[list_op[icontr]];
    }

  //Call the routine which does the real contraction
  trace_g_sdag_g_s(corr,t1,S1,t2,source,ncontr_2pts);
}

//generate the source
void generate_source()
{
  enum rnd_type type[5]={RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_Z2,RND_Z2,RND_Z4};
  generate_spindiluted_source(original_source,type[noise_type],twall);
}

//Generate a sequential source for S1
void generate_sequential_source(int ispec)
{
  int r=r_spec[ispec];
  
  master_printf("\nCreating the sequential source for spectator %d\n",ispec);
  for(int ivol=0;ivol<loc_vol;ivol++)
    { //put to zero everything but the slice
      if(glb_coord_of_loclx[ivol][0]!=(twall+tsep)%glb_size[0])
	memset(sequential_source[ivol],0,sizeof(colorspinspin));
      else
	{ //avoid to put g5, beacuse commutate with (i+-g5)/sqrt(2) and cancel with those of the QQ
	  memcpy(sequential_source[ivol],S0[r][iprop_of(ith_spec[ispec],imass_spec[ispec])][ivol],sizeof(colorspinspin));
	  for(int c=0;c<3;c++) //rotate as r because it's D^-1
	    rotate_spinspin_to_physical_basis(sequential_source[ivol][c],r,r);
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

  read_list_of_double_pairs("MassResidues",&nmass,&mass,&stopping_residues);
  read_list_of_doubles("NTheta",&ntheta,&theta);

  // 5) contraction list for two points
  
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

  // 6) three points functions
  
  sequential_source=nissa_malloc("Sequential source",loc_vol,colorspinspin);
  read_str_int("TSep",&tsep);
  read_str_int("NSpec",&nspec);
  ith_spec=nissa_malloc("ith_spec",nspec,int);
  imass_spec=nissa_malloc("imass_spec",nspec,int);
  r_spec=nissa_malloc("r_spec",nspec,int);
  for(int ispec=0;ispec<nspec;ispec++)
    {
      expect_str("iThetaMassR");
      read_int(&(ith_spec[ispec]));
      read_int(&(imass_spec[ispec]));
      read_int(&(r_spec[ispec]));
      master_printf(" spec %d: th=%g, m=%g, r=%d\n",ispec,theta[ith_spec[ispec]],mass[imass_spec[ispec]],r_spec[ispec]);
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

  //allocate gauge conf, Pmunu and all the needed spincolor and colorspinspin
  conf=nissa_malloc("or_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  sme_conf=nissa_malloc("sm_conf",loc_vol+bord_vol,quad_su3);
  Pmunu=nissa_malloc("Pmunu",loc_vol,as2t_su3);

  //Allocate all the S0 colorspinspin vectors
  npropS0=ntheta*nmass;
  S0[0]=nissa_malloc("S0[0]",npropS0,colorspinspin*);
  S0[1]=nissa_malloc("S0[1]",npropS0,colorspinspin*);
  for(int iprop=0;iprop<npropS0;iprop++)
    {
      S0[0][iprop]=nissa_malloc("S0[0]",loc_vol,colorspinspin);
      S0[1][iprop]=nissa_malloc("S0[1]",loc_vol,colorspinspin);
    }
	
	
	//Allocate all the S0_der colorspinspin vectors
	
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++) 
		{
			S0_der[0][i][j]=nissa_malloc("S0_der[0]",npropS0,colorspinspin*);
			S0_der[1][i][j]=nissa_malloc("S0_der[1]",npropS0,colorspinspin*);
			
			for(int iprop=0;iprop<npropS0;iprop++)
			{   S0_der[0][i][j][iprop]=nissa_malloc("S0_der[0]",loc_vol,colorspinspin);
				S0_der[1][i][j][iprop]=nissa_malloc("S0_der[1]",loc_vol,colorspinspin);							
			}
		}
		
	}


  //Allocate nmass spincolors, for the cgm solutions
  cgm_solution=nissa_malloc("cgm_solution",nmass,spincolor*);
  for(int imass=0;imass<nmass;imass++) cgm_solution[imass]=nissa_malloc("cgm_solution",loc_vol+bord_vol,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol+bord_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol+bord_vol,spincolor);
  
	//Allocate nmass spincolors, for the cgm "derivative" solutions
	/*cgm_der_solution=nissa_malloc("cgm_der_solution",nmass,spincolor*);
	for(int imass=0;imass<nmass;imass++) cgm_der_solution[imass]=nissa_malloc("cgm_der_solution",loc_vol+bord_vol,spincolor);*/
	temp_der_vec[0]=nissa_malloc("temp_der_vec[0]",loc_vol+bord_vol,spincolor);
	temp_der_vec[1]=nissa_malloc("temp_der_vec[1]",loc_vol+bord_vol,spincolor);
	
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  original_source=nissa_malloc("original_source",loc_vol,colorspinspin);

	//Allocate one spincolor for the der. source
	der_source=nissa_malloc("der_source",loc_vol+bord_vol,spincolor);
	

  //Allocate one colorspinspin for the chromo-contractions
  ch_colorspinspin=nissa_malloc("chromo-colorspinspin",loc_vol,colorspinspin);

  //Allocate all the S1 colorspinspin vectors
  npropS1=ntheta*nmass;
  S1=nissa_malloc("S1",loc_vol,colorspinspin*);
  for(int iprop=0;iprop<npropS0;iprop++) S1[iprop]=nissa_malloc("S1[i]",loc_vol,colorspinspin);
	
}
	
	//find a new conf
int read_conf_parameters(int *iconf)
{
  int ok_conf;

  do
    {
      //Gauge path
      read_str(conf_path,1024);
      
      //Twall
      read_int(&(twall));
      
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
  memset(old_theta,0,4*sizeof(double));
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,1);
}

//Finalization
void close_semileptonic()
{
  contr_time=contr_2pts_time+contr_3pts_time+contr_save_time;
  
  master_printf("\n");
  master_printf("Total time: %g, of which:\n",tot_time);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d contr. (%2.2gs avg) of which:\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
  master_printf("   * %02.2f%s to compute two points\n",contr_2pts_time*100.0/contr_time,"%");
  master_printf("   * %02.2f%s to compute three points\n",contr_3pts_time*100.0/contr_time,"%");
  master_printf("   * %02.2f%s to save correlations\n",contr_save_time*100.0/contr_time,"%");
  
  nissa_free(Pmunu);nissa_free(conf);nissa_free(sme_conf);
	
for(int iprop=0;iprop<npropS0;iprop++){nissa_free(S0[0][iprop]);nissa_free(S0[1][iprop]);nissa_free(S1[iprop]);} 
for (int i=0; i<3; i++)
for (int j=0; j<3; j++) 
for(int iprop=0;iprop<npropS0;iprop++)
{nissa_free(S0_der[0][i][j][iprop]);nissa_free(S0_der[1][i][j][iprop]);}
	nissa_free(S0[0]);nissa_free(S0[1]);nissa_free(S1);
for (int i=0; i<3; i++)
		for (int j=0; j<3; j++) 
		{nissa_free(S0_der[0][i][j]); nissa_free(S0_der[1][i][j]);}
nissa_free(temp_der_vec[0]);nissa_free(temp_der_vec[1]);
nissa_free(temp_vec[0]);nissa_free(temp_vec[1]);
  nissa_free(ch_colorspinspin);nissa_free(sequential_source);
  nissa_free(contr_2pts);nissa_free(ch_contr_2pts);
  nissa_free(contr_3pts);nissa_free(ch_contr_3pts);
  nissa_free(op1_2pts);nissa_free(op2_2pts);
  nissa_free(op1_3pts);nissa_free(op2_3pts);
  nissa_free(ch_op1_2pts);nissa_free(ch_op2_2pts);
  nissa_free(ch_op1_3pts);nissa_free(ch_op2_3pts);
  nissa_free(ith_spec);nissa_free(r_spec);nissa_free(imass_spec);
  for(int imass=0;imass<nmass;imass++) {nissa_free(cgm_solution[imass]);}
  nissa_free(cgm_solution);
  nissa_free(source);nissa_free(original_source);nissa_free(der_source);
  close_nissa();
}

//smear addditivily a colorspinspin
void smear_additive_colorspinspin(colorspinspin *out,colorspinspin *in,int ism_lev,int *jacobi_niter)
{
  spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  
  int nsme=jacobi_niter[ism_lev];
  if(ism_lev>0) nsme-=jacobi_niter[ism_lev-1];

  for(int id=0;id<4;id++)
    {
      for(int ivol=0;ivol<loc_vol;ivol++)
	get_spincolor_from_colorspinspin(temp[ivol],in[ivol],id);
      
      jacobi_smearing(temp,temp,sme_conf,jacobi_kappa,nsme);

      for(int ivol=0;ivol<loc_vol;ivol++)
	put_spincolor_into_colorspinspin(out[ivol],temp[ivol],id);
    }
  
  nissa_free(temp);
}

//calculate the standard propagators
void calculate_S0(int ism_lev_so)
{
  //smear additively the source
  master_printf("\nSource Smearing level: %d\n",ism_lev_so);
  smear_additive_colorspinspin(original_source,original_source,ism_lev_so,jacobi_niter_so);
  master_printf("\n");
  
 //loop over the source dirac index
  for(int id=0;id<4;id++)
    { 
      //put the g5
      for(int ivol=0;ivol<loc_vol;ivol++)
	{
	  get_spincolor_from_colorspinspin(source[ivol],original_source[ivol],id);
	  for(int id1=2;id1<4;id1++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic][ri]*=-1;
	}
      set_borders_invalid(source);
      
      for(int itheta=0;itheta<ntheta;itheta++)
	{
	  //adapt the border condition
	  put_theta[1]=put_theta[2]=put_theta[3]=theta[itheta];
	  adapt_theta(conf,old_theta,put_theta,1,1);
	  
	  //inverting
	  double part_time=-take_time();
	  inv_tmQ2_cgm(cgm_solution,conf,kappa,mass,nmass,niter_max,stopping_residues,source);
	  part_time+=take_time();ninv_tot++;inv_time+=part_time;
	  master_printf("Finished the inversion of S0 theta %d, dirac index %d in %g sec\n",itheta,id,part_time);

	  //reconstruct the doublet
	  for(int imass=0;imass<nmass;imass++)
	    {
	      reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass[imass],cgm_solution[imass]);
	      master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	      for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
		for(int i=0;i<loc_vol;i++)
		  put_spincolor_into_colorspinspin(S0[r][iprop_of(itheta,imass)][i],temp_vec[r][i],id);
	    }
	}
    }
  
  //rotate to physical basis
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
      rotate_vol_colorspinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
  master_printf("Propagators rotated\n");

  master_printf("\n");
}


// calculate the derivative standard propagators
void calculate_S0_derivative(int ism_lev_so)

{
	//smear additively the source
	master_printf("\nSource Smearing level: %d\n",ism_lev_so);
	smear_additive_colorspinspin(original_source,original_source,ism_lev_so,jacobi_niter_so);
	master_printf("\n");
	
	//loop over the source dirac index
	for(int id=0;id<4;id++)
    { 
		//put the g5
		for(int ivol=0;ivol<loc_vol;ivol++)
		{
			get_spincolor_from_colorspinspin(source[ivol],original_source[ivol],id);
			for(int id1=2;id1<4;id1++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic][ri]*=-1;
		}
		set_borders_invalid(source);
		
		//derivative operator applied D_i on the smeared source
		for (int idir=1;idir<4;idir++)
		{
			apply_nabla_i_spincolor(der_source, source, conf, idir);
					
		for(int itheta=0;itheta<ntheta;itheta++)
		{
			//adapt the border condition
			put_theta[1]=put_theta[2]=put_theta[3]=theta[itheta];
			adapt_theta(conf,old_theta,put_theta,1,1);
			
			//inverting
			double part_time=-take_time();
			//master_printf("calculate the solution with derivative source \n");
			inv_tmQ2_cgm(cgm_solution,conf,kappa,mass,nmass,niter_max,stopping_residues,der_source);
			part_time+=take_time();ninv_tot++;inv_time+=part_time;
			master_printf("Finished the inversion of S0_der theta %d, dirac index %d in %g sec\n",itheta,id,part_time);
			
			//reconstruct the doublet
			for(int imass=0;imass<nmass;imass++)
			{
			  reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass[imass],cgm_solution[imass]);
				master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
				//derivative_operator applied on the solution
				for (int jdir=1;jdir<4;jdir++) 
				{
					apply_nabla_i_spincolor(temp_der_vec[0], temp_vec[0], conf, jdir);
					apply_nabla_i_spincolor(temp_der_vec[1], temp_vec[1], conf, jdir);
					
					for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
					for(int i=0;i<loc_vol;i++)
						put_spincolor_into_colorspinspin(S0_der[r][idir-1][jdir-1][iprop_of(itheta,imass)][i],temp_der_vec[r][i],id);
			}
		}
    }
		}
	}
	
	//rotate to physical basis
	for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
		for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
			for (int idir=1;idir<4;idir++)
				for (int jdir=1;jdir<4;jdir++)
			rotate_vol_colorspinspin_to_physical_basis(S0_der[r][idir-1][jdir-1][ipropS0],!r,!r);
	master_printf("Propagators rotated\n");
	
	master_printf("\n");
}


//calculate the sequential propagators
void calculate_S1(int ispec,int ism_lev_se)
{
  //smear additively the seq
  master_printf("\nSeq Smearing level: %d (will be applied twice!)\n",ism_lev_se);
  smear_additive_colorspinspin(sequential_source,sequential_source,ism_lev_se,jacobi_niter_se);
  smear_additive_colorspinspin(sequential_source,sequential_source,ism_lev_se,jacobi_niter_se);
  master_printf("\n");
  
  //loop over se
  for(int id=0;id<4;id++)
    { 
      for(int ivol=0;ivol<loc_vol;ivol++) //avoid the g5 insertion
	get_spincolor_from_colorspinspin(source[ivol],sequential_source[ivol],id);
      set_borders_invalid(source);

      for(int itheta=0;itheta<ntheta;itheta++)
	{ //adapt the border condition
	  put_theta[1]=put_theta[2]=put_theta[3]=theta[itheta];
	  adapt_theta(conf,old_theta,put_theta,1,1);

	  double part_time=-take_time();
	  inv_tmQ2_cgm(cgm_solution,conf,kappa,mass,nmass,niter_max,stopping_residues,source);
	  part_time+=take_time();ninv_tot++;inv_time+=part_time;
	  master_printf("Finished the inversion of S1 theta %d, seq sme lev %d, dirac index %d in %g sec\n",itheta,ism_lev_se,id,part_time);
	  
	  for(int imass=0;imass<nmass;imass++)
	    { //reconstruct the doublet: r(S1)=!r(spec), so we have to multiply by Q+ if r(spec)==1 and Q- if 0
	      double reco_mass=-mass[imass];
	      if(r_spec[ispec]==1) reco_mass=-reco_mass;
	      //use temp_vec[0] as temporary storage
	      apply_tmQ(temp_vec[0],conf,kappa,reco_mass,cgm_solution[imass]);
	      master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	      for(int i=0;i<loc_vol;i++) put_spincolor_into_colorspinspin(S1[iprop_of(itheta,imass)][i],temp_vec[0][i],id);
	    }
	}
    }

  //put the (1+-ig5)/sqrt(2) factor. On the source rotate as r_spec, on the sink as !r_spec
  for(int ipropS1=0;ipropS1<npropS1;ipropS1++) //but, being D^-1, everything is swapped
    rotate_vol_colorspinspin_to_physical_basis(S1[ipropS1],!(r_spec[ispec]),(r_spec[ispec]));
  master_printf("Propagators rotated\n");

  master_printf("\n");
}

//Calculate and print to file the 2pts
 void calculate_all_2pts_ten(int ism_lev_so,int ism_lev_si)
{
  contr_2pts_time-=take_time();

  for(int r=0;r<2;r++)
    for(int iprop=0;iprop<npropS0;iprop++)
      smear_additive_colorspinspin(S0[r][iprop],S0[r][iprop],ism_lev_si,jacobi_niter_si);
	
	for(int r=0;r<2;r++)
		for(int iprop=0;iprop<npropS0;iprop++)
			for (int idir=1;idir<4;idir++)
				for (int jdir=1;jdir<4;jdir++)
					smear_additive_colorspinspin(S0_der[r][idir-1][jdir-1][iprop],S0_der[r][idir-1][jdir-1][iprop],ism_lev_si,jacobi_niter_si);
	
	char path[1024];
  sprintf(path,"%s/2pts_%02d_%02d",outfolder,jacobi_niter_so[ism_lev_so],jacobi_niter_si[ism_lev_si]);
  FILE *fout=open_text_file_for_output(path);

  for(int ispec=0;ispec<nspec;ispec++)
    {
      int ith1=ith_spec[ispec];
		
		for (int idir=1;idir<4;idir++)			
		for (int jdir=1;jdir<4;jdir++)
			
      for(int ith2=0;ith2<ntheta;ith2++)
	for(int im2=0;im2<nmass;im2++)
  	  {
	   int ip2=iprop_of(ith2,im2);
	    for(int r2=0;r2<2;r2++)
	      {
		  if(nch_contr_2pts>0)
		      unsafe_apply_chromo_operator_to_colorspinspin(ch_colorspinspin,Pmunu,S0[r2][ip2]);
		  
		  for(int im1=0;im1<nmass;im1++)
		    {
		      int ip1=iprop_of(ith1,im1);
		      
		      for(int r1=0;r1<2;r1++)
		        {
			    
			    if(rank==0)
			     {
					 fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d idir=%d jdir=%d \n",
							 mass[im1],theta[ith1],r1,mass[im2],theta[ith2],r2, idir, jdir);
					 fprintf(fout," smear_source=%d smear_sink=%d\n",jacobi_niter_so[ism_lev_so],jacobi_niter_si[ism_lev_si]);	
				 }
		  
			   meson_two_points(contr_2pts,op1_2pts,S0_der[r1][idir-1][jdir-1][ip1],op2_2pts,S0[r2][ip2],ncontr_2pts);
			    ncontr_tot+=ncontr_2pts;
			    
			    contr_save_time-=take_time();
			    print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,twall,"",1.0);
			    contr_save_time+=take_time();
			    
			    if(nch_contr_2pts>0)
			    {
				meson_two_points(ch_contr_2pts,ch_op1_2pts,S0[r1][ip1],ch_op2_2pts,ch_colorspinspin,nch_contr_2pts);
				ncontr_tot+=nch_contr_2pts;
				
				contr_save_time-=take_time();
				print_contractions_to_file(fout,nch_contr_2pts,ch_op1_2pts,ch_op2_2pts,ch_contr_2pts,twall,"CHROMO-",1.0);
				contr_save_time+=take_time();
			    }
			    if(rank==0) fprintf(fout,"\n");
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
  for(int ith2=0;ith2<ntheta;ith2++)
    for(int im2=0;im2<nmass;im2++)
    {
	int ip2=iprop_of(ith2,im2);
	if(nch_contr_3pts>0)
	    unsafe_apply_chromo_operator_to_colorspinspin(ch_colorspinspin,Pmunu,S1[ip2]);
	
	for(int ith1=0;ith1<ntheta;ith1++)
	    for(int im1=0;im1<nmass;im1++)
	    {
		int ip1=iprop_of(ith1,im1);
		
		if(rank==0)
		{
		    fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d,",mass[im1],theta[ith1],r1,mass[im2],theta[ith2],r2);
		    fprintf(fout," smear_source=%d smear_seq=%d\n",jacobi_niter_so[ism_lev_so],jacobi_niter_se[ism_lev_se]);
		}
		
		meson_two_points(contr_3pts,op1_3pts,S0[r1][ip1],op2_3pts,S1[ip2],ncontr_3pts);
		ncontr_tot+=ncontr_3pts;
		
		contr_save_time-=take_time();
		print_contractions_to_file(fout,ncontr_3pts,op1_3pts,op2_3pts,contr_3pts,twall,"",1.0);
		contr_save_time+=take_time();
		
		if(nch_contr_3pts>0)
		{
		    meson_two_points(ch_contr_3pts,ch_op1_3pts,S0[r1][ip1],ch_op2_3pts,ch_colorspinspin,nch_contr_3pts);
		    ncontr_tot+=nch_contr_3pts;
		    
		    contr_save_time-=take_time();
		    print_contractions_to_file(fout,nch_contr_3pts,ch_op1_3pts,ch_op2_3pts,ch_contr_3pts,twall,"CHROMO-",1.0);
		    contr_save_time+=take_time();
		}
		if(rank==0) fprintf(fout,"\n");
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

  for(int ith2=0;ith2<ntheta;ith2++)
    for(int im2=0;im2<nmass;im2++)
      {
	int ip2=iprop_of(ith2,im2);
	contract_with_source(contr_2pts,S1[ip2],op2_2pts,original_source);
	
	if(rank==0)
	  {
	    fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",
		    mass[imass_spec[ispec]],theta[ith_spec[ispec]],r_spec[ispec],mass[im2],theta[ith2],r_spec[ispec]);

	    fprintf(fout,"\n");
	    
	    for(int icontr=0;icontr<ncontr_2pts;icontr++)
	      if(op1_2pts[icontr]==5)
		fprintf(fout," # P5%s\t%+016.16g\t%+016.16g\n",gtag[op2_2pts[icontr]],(contr_2pts+icontr*glb_size[0])[twall][0],(contr_2pts+icontr*glb_size[0])[twall][1]);
	    fprintf(fout,"\n");
	  }
      }
  
  if(rank==0) fclose(fout);
}

//check if the time is enough
int check_remaining_time()
{
  int enough_time;

  //check remaining time                                                                                                                                                                        
  double temp_time=take_time()+tot_time;
  double ave_time=temp_time/nanalized_conf;
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

  if(narg<2) crash("Use: %s input_file",arg[0]);

  tot_time-=take_time();
  initialize_semileptonic(arg[1]);
  
  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf && enough_time && read_conf_parameters(&iconf))
    {
      setup_conf();
      generate_source();
      
      //loop on smearing of the source
      for(int sm_lev_so=0;sm_lev_so<nsm_lev_so;sm_lev_so++)
	{
	  
		calculate_S0(sm_lev_so);
		calculate_S0_derivative(sm_lev_so);
	  
	  //loop on spectator
	  if(ncontr_3pts!=0 && nch_contr_3pts!=0)
	  for(int ispec=0;ispec<nspec;ispec++)
	    {
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
	  for(int sm_lev_si=0;sm_lev_si<nsm_lev_si;sm_lev_si++) calculate_all_2pts_ten(sm_lev_so,sm_lev_si);
	}

      nanalized_conf++;

      enough_time=check_remaining_time();
    }
  
  if(iconf==ngauge_conf) master_printf("Finished all the conf!\n");
  
  tot_time+=take_time();
  close_semileptonic();
  
  return 0;
}
