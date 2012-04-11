#include "nissa.h"

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

//smearing parameters
double jacobi_kappa,ape_alpha;
int ape_niter;
int *jacobi_niter_so,nsm_lev_so;
int *jacobi_niter_si,nsm_lev_si;
int *jacobi_niter_se,nsm_lev_se;

//vectors for the spinor data
int npropS0;
colorspinspin **S0[2];

//cgmms inverter parameters
double stopping_residue;
double minimal_residue;
int stopping_criterion;
int niter_max;

//two points contractions
int ncontr_2pts;
int *op_2pts[2];

//two points chromo contractions
int nch_contr_2pts;
int *ch_op_2pts[2];

//sequential props info
int nspec;
int tsep;
int *imass_spec,*r_spec,*ith_spec;
int *start_imass3pts;
colorspinspin *sequential_source;

//sequential propagators
int npropS1;
colorspinspin **S1;

//three points contractions
int ncontr_3pts;
int *op_3pts[2];

//two points chromo contractions
int nch_contr_3pts;
int *ch_op_3pts[2];

//timings
int ninv_tot=0,ncontr_tot=0;
int wall_time;
double tot_time=0,inv_time=0;
double load_time=0;
double contr_2pts_time=0;
double contr_3pts_time=0;
double contr_save_time=0;

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
  nissa_loc_vol_loop(ivol)
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
  {
    int spat_ranks=rank_tot/nrank_dir[0];
    if(rank%spat_ranks!=plan_rank[0]) crash("plan and proj rank do not agree");
  }
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

  read_list_of_doubles("NMass",&nmass,&mass);
  read_list_of_doubles("NTheta",&ntheta,&theta);

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
  
  // 6) contraction list for two points
  
  read_str_int("NContrTwoPoints",&ncontr_2pts);
  op_2pts[0]=nissa_malloc("op0_2pts",ncontr_2pts,int);
  op_2pts[1]=nissa_malloc("op1_2pts",ncontr_2pts,int);
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op_2pts[0][icontr]));
      read_int(&(op_2pts[1][icontr]));
      
      master_printf(" contr.%d %d %d\n",icontr,op_2pts[0][icontr],op_2pts[1][icontr]);
    }
  
  read_str_int("NChromoContrTwoPoints",&nch_contr_2pts);
  ch_op_2pts[0]=nissa_malloc("ch_op0_2pts",ncontr_2pts,int);
  ch_op_2pts[1]=nissa_malloc("ch_op1_2pts",ncontr_2pts,int);
  for(int icontr=0;icontr<nch_contr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(ch_op_2pts[0][icontr]));
      read_int(&(ch_op_2pts[1][icontr]));
      
      master_printf(" ch-contr.%d %d %d\n",icontr,ch_op_2pts[0][icontr],ch_op_2pts[1][icontr]);
    }

  // 7) three points functions
  
  sequential_source=nissa_malloc("Sequential source",loc_vol,colorspinspin);
  read_str_int("TSep",&tsep);
  if(tsep<0||tsep>=glb_size[0]) crash("tsep=%d while T=%d",twall,glb_size[0]);
  
  read_str_int("NSpec",&nspec);
  ith_spec=nissa_malloc("ith_spec",nspec,int);
  imass_spec=nissa_malloc("imass_spec",nspec,int);
  r_spec=nissa_malloc("r_spec",nspec,int);
  start_imass3pts=nissa_malloc("start_imass3pts",nspec,int);
  for(int ispec=0;ispec<nspec;ispec++)
    {
      expect_str("iThetaMassR");
      read_int(&(ith_spec[ispec]));
      read_int(&(imass_spec[ispec]));
      read_int(&(r_spec[ispec]));
      if(ith_spec[ispec]>=ntheta||ith_spec[ispec]<0||imass_spec[ispec]>=nmass||imass_spec[ispec]<0||r_spec[ispec]>=2||r_spec[ispec]<0)
	  crash("requiring uncomputed propagator as sequential source");
      master_printf(" spec %d: th=%lg, m=%lg, r=%d\n",ispec,theta[ith_spec[ispec]],mass[imass_spec[ispec]],r_spec[ispec]);
      
      read_str_int("StartIMass",&(start_imass3pts[ispec]));
      if(start_imass3pts[ispec]<0||start_imass3pts[ispec]>=nmass)
	crash("Requiring to start from uncomputed mass");
      master_printf(" mass range for sequential: [%d,%d]\n",start_imass3pts[ispec],nmass-1);
    }
  read_str_int("NContrThreePoints",&ncontr_3pts);
  op_3pts[0]=nissa_malloc("op0_3pts",ncontr_3pts,int);
  op_3pts[1]=nissa_malloc("op1_3pts",ncontr_3pts,int);
  for(int icontr=0;icontr<ncontr_3pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op_3pts[0][icontr]));
      read_int(&(op_3pts[1][icontr]));
      
      master_printf(" contr.%d %d %d\n",icontr,op_3pts[0][icontr],op_3pts[1][icontr]);
    }

  read_str_int("NChromoContrThreePoints",&nch_contr_3pts);
  ch_op_3pts[0]=nissa_malloc("ch_op0_3pts",nch_contr_3pts,int);
  ch_op_3pts[1]=nissa_malloc("ch_op1_3pts",nch_contr_3pts,int);
  for(int icontr=0;icontr<nch_contr_3pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(ch_op_3pts[0][icontr]));
      read_int(&(ch_op_3pts[1][icontr]));
      
      master_printf(" ch-contr.%d %d %d\n",icontr,ch_op_3pts[0][icontr],ch_op_3pts[1][icontr]);
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

  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  original_source=nissa_malloc("original_source",loc_vol,colorspinspin);

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
      if(twall<0||twall>=glb_size[0]) crash("twall=%d while T=%d",twall,glb_size[0]);
      
      //Out folder
      read_str(outfolder,1024);
      
      //Check if the conf exist
      master_printf("Considering configuration %s\n",conf_path);
      ok_conf=!(dir_exists(outfolder));
      if(ok_conf)
	{
	  int ris=create_dir(outfolder);
	  if(ris==0) master_printf("Configuration %s not already analized, starting.\n",conf_path);
	  else
	    {
	      ok_conf=0;
	      master_printf("Configuration %s taken by someone else.\n",conf_path);
	    }
	}
      else
	master_printf("Configuration %s already analized, skipping.\n",conf_path);
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
  communicate_lx_quad_su3_borders(conf);
  communicate_lx_gauge_edges(conf);
  Pmunu_term(Pmunu,conf);
  
  //prepare the smerded version
  ape_spatial_smear_conf(sme_conf,conf,ape_alpha,ape_niter);
  communicate_lx_quad_su3_borders(conf);
  communicate_lx_quad_su3_borders(sme_conf);
  
  //compute plaquette
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  master_printf("smerded plaq: %.18g\n",global_plaquette_lx_conf(sme_conf));
  
  //put the anti-periodic condition on the temporal border
  memset(old_theta,0,4*sizeof(double));
  memset(put_theta,0,4*sizeof(double));
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,1);
}

//Finalization
void close_semileptonic()
{
  double contr_time=contr_save_time+contr_2pts_time+contr_3pts_time;

  master_printf("\n");
  master_printf("Total time: %g",tot_time);
  if(ninv_tot>0)
    {
      master_printf(", of which:\n");
      master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
      master_printf(" - %02.2f%s to perform %d contractions (%2.2gs avg), of which:\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
      master_printf("   * %02.2f%s to compute 2pts\n",contr_2pts_time*100.0/contr_time,"%");
      master_printf("   * %02.2f%s to compute 3pts\n",contr_3pts_time*100.0/contr_time,"%");
      master_printf("   * %02.2f%s to save correlations\n",contr_save_time*100.0/contr_time,"%");
    }
  else master_printf(".\n");
  
  nissa_free(Pmunu);nissa_free(conf);nissa_free(sme_conf);
  for(int iprop=0;iprop<npropS0;iprop++){nissa_free(S0[0][iprop]);nissa_free(S0[1][iprop]);nissa_free(S1[iprop]);}
  nissa_free(S0[0]);nissa_free(S0[1]);nissa_free(S1);
  nissa_free(sequential_source);
  nissa_free(op_2pts[0]);nissa_free(op_2pts[1]);
  nissa_free(op_3pts[0]);nissa_free(op_3pts[1]);
  nissa_free(ch_op_2pts[0]);nissa_free(ch_op_2pts[1]);
  nissa_free(ch_op_3pts[0]);nissa_free(ch_op_3pts[1]);
  nissa_free(ith_spec);nissa_free(r_spec);nissa_free(imass_spec);
  nissa_free(start_imass3pts);
  nissa_free(source);nissa_free(original_source);
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
      nissa_loc_vol_loop(ivol)
	get_spincolor_from_colorspinspin(temp[ivol],in[ivol],id);
      
      jacobi_smearing(temp,temp,sme_conf,jacobi_kappa,nsme);

      nissa_loc_vol_loop(ivol)
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
  
  //Allocate nmass spincolors, for the cgmms solutions
  spincolor **cgmms_solution,*temp_vec[2];
  cgmms_solution=nissa_malloc("cgmms_solution",nmass,spincolor*);
  for(int imass=0;imass<nmass;imass++) cgmms_solution[imass]=nissa_malloc("cgmms_solution",loc_vol+bord_vol,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol,spincolor);
  
 //loop over the source dirac index
  for(int id=0;id<4;id++)
    { 
      //put the g5
      nissa_loc_vol_loop(ivol)
	{
	  get_spincolor_from_colorspinspin(source[ivol],original_source[ivol],id);
	  for(int id1=2;id1<4;id1++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic][ri]*=-1;
	}
      
      communicate_lx_spincolor_borders(source);
      
      for(int itheta=0;itheta<ntheta;itheta++)
	{
	  //adapt the border condition
	  put_theta[1]=put_theta[2]=put_theta[3]=theta[itheta];
	  adapt_theta(conf,old_theta,put_theta,1,1);
	  
	  //inverting
	  double part_time=-take_time();
	  inv_tmQ2_cgmms(cgmms_solution,source,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	  part_time+=take_time();ninv_tot++;inv_time+=part_time;
	  master_printf("Finished the inversion of S0 theta %d, dirac index %d in %g sec\n",itheta,id,part_time);

	  //reconstruct the doublet
	  for(int imass=0;imass<nmass;imass++)
	    {
	      reconstruct_tm_doublet(temp_vec[0],temp_vec[1],cgmms_solution[imass],conf,kappa,mass[imass]);
	      master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	      for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
		nissa_loc_vol_loop(ivol)
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
  
  //free vectors
  for(int imass=0;imass<nmass;imass++) nissa_free(cgmms_solution[imass]);
  nissa_free(cgmms_solution);
  nissa_free(temp_vec[0]);nissa_free(temp_vec[1]);
}

//calculate the sequential propagators
void calculate_S1(int ispec,int ism_lev_se)
{
  int nmass_3pts=nmass-start_imass3pts[ispec];
  double *mass_3pts=mass+start_imass3pts[ispec];
  
  //smear additively the seq
  master_printf("\nSeq Smearing level: %d (will be applied twice!)\n",ism_lev_se);
  smear_additive_colorspinspin(sequential_source,sequential_source,ism_lev_se,jacobi_niter_se);
  smear_additive_colorspinspin(sequential_source,sequential_source,ism_lev_se,jacobi_niter_se);
  master_printf("\n");
  
  //Allocate nmass spincolors, for the cgmms solutions
  spincolor **cgmms_solution,*temp_vec[2];
  cgmms_solution=nissa_malloc("cgmms_solution",nmass_3pts,spincolor*);
  for(int imass=0;imass<nmass_3pts;imass++) cgmms_solution[imass]=nissa_malloc("cgmms_solution",loc_vol+bord_vol,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol,spincolor);
  
  //loop over se
  for(int id=0;id<4;id++)
    { 
      nissa_loc_vol_loop(ivol) //avoid the g5 insertion
	get_spincolor_from_colorspinspin(source[ivol],sequential_source[ivol],id);
      communicate_lx_spincolor_borders(source);

      for(int itheta=0;itheta<ntheta;itheta++)
	{ //adapt the border condition
	  put_theta[1]=put_theta[2]=put_theta[3]=theta[itheta];
	  adapt_theta(conf,old_theta,put_theta,1,1);

	  double part_time=-take_time();
	  inv_tmQ2_cgmms(cgmms_solution,source,conf,kappa,mass_3pts,nmass_3pts,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	  part_time+=take_time();ninv_tot++;inv_time+=part_time;
	  master_printf("Finished the inversion of S1 theta %d, seq sme lev %d, dirac index %d in %g sec\n",itheta,ism_lev_se,id,part_time);
	  
	  for(int imass=0;imass<nmass_3pts;imass++)
	    { //reconstruct the doublet: r(S1)=!r(spec), so we have to multiply by Q+ if r(spec)==1 and Q- if 0
	      double reco_mass=-mass[imass+start_imass3pts[ispec]];
	      if(r_spec[ispec]==1) reco_mass=-reco_mass;
	      //use temp_vec[0] as temporary storage
	      apply_tmQ(temp_vec[0],cgmms_solution[imass],conf,kappa,reco_mass);
	      master_printf("Mass %d (%g) reconstructed \n",imass,mass_3pts[imass]);
	      nissa_loc_vol_loop(ivol) put_spincolor_into_colorspinspin(S1[iprop_of(itheta,imass+start_imass3pts[ispec])][ivol],temp_vec[0][ivol],id);
	    }
	}
    }

  //put the (1+-ig5)/sqrt(2) factor. On the source rotate as r_spec, on the sink as !r_spec
  for(int imass=0;imass<nmass_3pts;imass++) //but, being D^-1, everything is swapped
    for(int itheta=0;itheta<ntheta;itheta++)
      rotate_vol_colorspinspin_to_physical_basis(S1[iprop_of(itheta,imass+start_imass3pts[ispec])],!(r_spec[ispec]),(r_spec[ispec]));
  master_printf("Propagators rotated\n");

  master_printf("\n");

  //free vectors
  for(int imass=0;imass<nmass_3pts;imass++) nissa_free(cgmms_solution[imass]);
  nissa_free(cgmms_solution);
  nissa_free(temp_vec[0]);nissa_free(temp_vec[1]);
}

//Calculate all two points contractions
void two_points(int ism_lev_so,int ism_lev_si)
{
  //semar
  for(int r=0;r<2;r++)
    for(int iprop=0;iprop<npropS0;iprop++)
      smear_additive_colorspinspin(S0[r][iprop],S0[r][iprop],ism_lev_si,jacobi_niter_si);
  
  //take initial time
  contr_2pts_time-=take_time();
  
  //find the number of propagator combinations to be hold on each x0=0 rank
  int nprop_combo=nspec*ntheta*nmass*nmass*2*2;
  int nrank_x0=rank_tot/nrank_dir[0];
  int nprop_combo_per_rank_x0=(int)ceil((double)nprop_combo/nrank_x0);
  
  //allocate space for contractions
  complex *glb_contr   =nissa_malloc(   "glb_2pts_contr",nprop_combo_per_rank_x0*glb_size[0]*ncontr_2pts,complex);
  complex *glb_ch_contr=nissa_malloc("glb_ch_2pts_contr",nprop_combo_per_rank_x0*glb_size[0]*nch_contr_2pts,complex);
  
  //if needed take into account chromo contractions
  colorspinspin *S2[nmass*ntheta];
  if(nch_contr_2pts>0)
    {
      //Allocate one colorspinspin for the chromo-contractions
      colorspinspin *ch_colorspinspin=nissa_malloc("chromo-colorspinspin",loc_vol,colorspinspin);
      
      //apply Pmunu      
      for(int iprop=0;iprop<nmass*ntheta;iprop++)
	{
	  S2[iprop]=nissa_malloc("S2",loc_vol,colorspinspin);
	  
	  //apply to r=0
	  unsafe_apply_chromo_operator_to_colorspinspin(ch_colorspinspin,Pmunu,S0[0][iprop]);
	  memcpy(S1[iprop],ch_colorspinspin,loc_vol*sizeof(colorspinspin));
	  
	  //apply to r=1
	  unsafe_apply_chromo_operator_to_colorspinspin(ch_colorspinspin,Pmunu,S0[1][iprop]);
	  memcpy(S2[iprop],ch_colorspinspin,loc_vol*sizeof(colorspinspin));
	}
      
      //free
      nissa_free(ch_colorspinspin);
    }
  
  //create the list of propagators
  colorspinspin *SA[nspec*nmass*2],*SB[ntheta*nmass*2],*SC[ntheta*nmass*2];
  int iA=0,iB=0,iC=0;
  for(int ispec=0;ispec<nspec;ispec++)
    for(int im1=0;im1<nmass;im1++)
      for(int r1=0;r1<2;r1++)
	SA[iA++]=S0[r1][iprop_of(ith_spec[ispec],im1)];
  for(int ith2=0;ith2<ntheta;ith2++)
    for(int im2=0;im2<nmass;im2++)
      {
	int ip2=iprop_of(ith2,im2);
	for(int r2=0;r2<2;r2++)
	  SB[iB++]=S0[r2][ip2];
	if(nch_contr_2pts>0)
	  {
	    SC[iC++]=S1[ip2]; //will use S1 for r=0 ch-prop
	    SC[iC++]=S2[ip2]; // and use S2 for r=1
	  }
      }
  
  //create the list of propagator combinations
  intpair nprop={nspec*nmass*2,2*ntheta*nmass};
  intpair *prop_combo=nissa_malloc("prop_combo",nprop_combo,intpair);
  int icombo=0;
  for(int ispec=0;ispec<nspec;ispec++)
    for(int ith2=0;ith2<ntheta;ith2++)
      for(int im2=0;im2<nmass;im2++)
	for(int r2=0;r2<2;r2++)
	  for(int im1=0;im1<nmass;im1++)
	    for(int r1=0;r1<2;r1++)
	      {
		prop_combo[icombo][0]=r1+2*(im1+nmass*ispec);
		prop_combo[icombo][1]=r2+2*(im2+nmass*ith2);
		icombo++;
	      }
  
  //perform all the contractions
  lot_of_mesonic_contractions(glb_contr,op_2pts,ncontr_2pts,SA,SB,nprop,prop_combo,nprop_combo,twall);

  //perform chromo contractions
  if(nch_contr_2pts>0)
    lot_of_mesonic_contractions(glb_ch_contr,ch_op_2pts,nch_contr_2pts,SA,SC,nprop,prop_combo,nprop_combo,twall);
  
  contr_2pts_time+=take_time();
  
  //save all the output
  contr_save_time-=take_time();
  
  //loop over corelation functions
  int ioff=0;
  int ich_off=0;
  
  //create output, opening it on rank 0 where first bunch of corr is stored
  char path[1024];
  sprintf(path,"%s/2pts_%02d_%02d",outfolder,jacobi_niter_so[ism_lev_so],jacobi_niter_si[ism_lev_si]);
  FILE *fout=open_text_file_for_output(path);
  int ist_combo=0,irank=0;
  
  if(rank_coord[0]==0)
    {
      for(int ispec=0;ispec<nspec;ispec++)
	{
	  int ith1=ith_spec[ispec];
	  for(int ith2=0;ith2<ntheta;ith2++)
	    for(int im2=0;im2<nmass;im2++)
	      for(int r2=0;r2<2;r2++)
		for(int im1=0;im1<nmass;im1++)
		  for(int r1=0;r1<2;r1++)
		    {
		      //if next rank must write, close, pass to next rank and reopen 
		      if(ist_combo==nprop_combo_per_rank_x0)
			{
			  //close and reopen
			  if(plan_rank[0]==irank) fclose(fout);
			  irank++;
			  
			  //wait
			  MPI_Barrier(plan_comm[0]);
			  
			  if(plan_rank[0]==irank) fout=fopen(path,"a");
			  //reset ncombo
			  ist_combo=ioff=ich_off=0;
			}
		      
		      //on appropriate rank, write
		      if(plan_rank[0]==irank)
			{
			  //print combination header
			  fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d",mass[im1],theta[ith1],r1,mass[im2],theta[ith2],r2);
			  fprintf(fout," smear_source=%d smear_sink=%d\n\n",jacobi_niter_so[ism_lev_so],jacobi_niter_si[ism_lev_si]);
			  
			  //loop over contraction of the combo
			  for(int icontr=0;icontr<ncontr_2pts;icontr++)
			    {
			      //print the contraction header
			      fprintf(fout," # %s%s\n",gtag[op_2pts[1][icontr]],gtag[op_2pts[0][icontr]]);
			      
			      //print the contraction
			      for(int t=0;t<glb_size[0];t++)
				{
				  fprintf(fout,"%+016.16g\t%+016.16g\n",glb_contr[ioff][0],glb_contr[ioff][1]);
				  ioff++;
				}
			      
			      //paragraph end
			      fprintf(fout,"\n");
			    }
			  
			  //loop over ch-contraction of the combo
			  for(int icontr=0;icontr<nch_contr_2pts;icontr++)
			    {
			      fprintf(fout," # CHROMO-%s%s\n",gtag[ch_op_2pts[1][icontr]],gtag[ch_op_2pts[0][icontr]]);
			      
			      //print the contraction
			      for(int t=0;t<glb_size[0];t++)
				{
				  fprintf(fout,"%+016.16g\t%+016.16g\n",glb_ch_contr[ich_off][0],glb_ch_contr[ich_off][1]);
				  ich_off++;
				}
			      
			      //paragraph end
			      fprintf(fout,"\n");
			    }
			}
		      //increment the number of stored combo
		      ist_combo++;
		    }
	}
      
      //close the output
      if(plan_rank[0]==irank) fclose(fout);
    }
  
  contr_save_time+=take_time();
  
  //free all the vectors
  
  if(nch_contr_2pts>0)
      for(int iprop=0;iprop<nmass*ntheta;iprop++)
	nissa_free(S2[iprop]);
  
  nissa_free(prop_combo);
  
  nissa_free(glb_contr);
  nissa_free(glb_ch_contr);

  ncontr_tot+=(ncontr_2pts+nch_contr_2pts)*nprop_combo;
}

//Calculate all three points contractions
void three_points(int ispec,int ism_lev_so,int ism_lev_se)
{
  //take initial time
  contr_3pts_time-=take_time();
  
  //Allocate one colorspinspin for the chromo-contractions
  colorspinspin *ch_colorspinspin=nissa_malloc("chromo-colorspinspin",loc_vol,colorspinspin);
  
  //find the number of propagator combinations to be hold on each x0=0 rank
  double *mass_3pts=mass+start_imass3pts[ispec];
  int nmass_3pts=nmass-start_imass3pts[ispec];
  int nprop_combo=ntheta*ntheta*nmass*nmass_3pts;
  int nrank_x0=rank_tot/nrank_dir[0];
  int nprop_combo_per_rank_x0=(int)ceil((double)nprop_combo/nrank_x0);
  
  //allocate space for contractions
  complex *glb_contr   =nissa_malloc(   "glb_3pts_contr",nprop_combo_per_rank_x0*glb_size[0]*ncontr_3pts,complex);
  complex *glb_ch_contr=nissa_malloc("glb_ch_3pts_contr",nprop_combo_per_rank_x0*glb_size[0]*nch_contr_3pts,complex);
  
  //create the list of propagator combinations
  intpair nprop={ntheta*nmass,ntheta*nmass_3pts};
  intpair *prop_combo=nissa_malloc("prop_combo",nprop_combo,intpair);
  int icombo=0;
  colorspinspin *S1S[nprop[1]];
  int ipropA=0;
  for(int itheta1=0;itheta1<ntheta;itheta1++)
    for(int imass1=0;imass1<nmass_3pts;imass1++)
      S1S[ipropA++]=S1[iprop_of(itheta1,imass1+start_imass3pts[ispec])];
  for(int iprop1=0;iprop1<ntheta*nmass_3pts;iprop1++)
    for(int iprop0=0;iprop0<nprop[0];iprop0++)
      {
	prop_combo[icombo][0]=iprop0;
	prop_combo[icombo][1]=iprop1;
	icombo++;
      }
  
  //perform all the contractions
  int r1=r_spec[ispec];
  lot_of_mesonic_contractions(glb_contr,op_3pts,ncontr_3pts,S0[r1],S1S,nprop,prop_combo,nprop_combo,twall);

  //perform chromo contractions
  if(nch_contr_3pts>0)
    {
      //apply Pmunu to the prop if needed  
      for(int itheta1=0;itheta1<ntheta;itheta1++)
	for(int imass1=0;imass1<nmass_3pts;imass1++)
	{
	  int iprop=iprop_of(itheta1,imass1+start_imass3pts[ispec]);
	  unsafe_apply_chromo_operator_to_colorspinspin(ch_colorspinspin,Pmunu,S1[iprop]);
	  memcpy(S1[iprop],ch_colorspinspin,loc_vol*sizeof(colorspinspin));
	}
      //perform all the contractions
      lot_of_mesonic_contractions(glb_ch_contr,ch_op_3pts,nch_contr_3pts,S0[r1],S1S,nprop,prop_combo,nprop_combo,twall);
    }
  
  contr_3pts_time+=take_time();
  
  //save all the output
  contr_save_time-=take_time();
  
  //loop over corelation functions
  int ioff=0;
  int ich_off=0;
  
  //create output, opening it on rank 0 where first bunch of corr is stored
  char path[1024];
  sprintf(path,"%s/3pts_sp%d_%02d_%02d",outfolder,ispec,jacobi_niter_so[ism_lev_so],jacobi_niter_se[ism_lev_se]);
  FILE *fout=open_text_file_for_output(path);
  int ist_combo=0,irank=0;
  
  //fix r2
  int r2=!r1;
 
  if(rank_coord[0]==0)
    {
      //loop over 2nd prop
      for(int ith2=0;ith2<ntheta;ith2++)
	for(int im2=0;im2<nmass_3pts;im2++)
	  //loop over 1st prop
	  for(int ith1=0;ith1<ntheta;ith1++)
	    for(int im1=0;im1<nmass;im1++)
	      {
		//if next rank must write, close, pass to next rank and reopen 
		if(ist_combo==nprop_combo_per_rank_x0)
		  {
		    //close and reopen
		    if(plan_rank[0]==irank) fclose(fout);
		    irank++;
		    
		    //wait
		    MPI_Barrier(plan_comm[0]);
		    
		    if(plan_rank[0]==irank) fout=fopen(path,"a");
		    //reset ncombo
		    ist_combo=ioff=ich_off=0;
		  }
		
		//on appropriate rank, write
		if(plan_rank[0]==irank)
		  {
		    //print combination header
		    fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d,",mass[im1],theta[ith1],r1,mass_3pts[im2],theta[ith2],r2);
		    fprintf(fout," smear_source=%d smear_seq=%d\n\n",jacobi_niter_so[ism_lev_so],jacobi_niter_se[ism_lev_se]);
		    
		    //loop over contraction of the combo
		    for(int icontr=0;icontr<ncontr_3pts;icontr++)
		      {
			//print the contraction header
			fprintf(fout," # %s%s\n",gtag[op_3pts[1][icontr]],gtag[op_3pts[0][icontr]]);
			
			//print the contraction
			for(int t=0;t<glb_size[0];t++)
			  {
			    fprintf(fout,"%+016.16g\t%+016.16g\n",glb_contr[ioff][0],glb_contr[ioff][1]);
			    ioff++;
			  }
			
			//paragraph end
			fprintf(fout,"\n");
		      }
		    
		    //loop over ch-contraction of the combo
		    for(int icontr=0;icontr<nch_contr_3pts;icontr++)
		      {
			fprintf(fout," # CHROMO-%s%s\n",gtag[ch_op_3pts[1][icontr]],gtag[ch_op_3pts[0][icontr]]);
			
			//print the contraction
			for(int t=0;t<glb_size[0];t++)
			  {
			    fprintf(fout,"%+016.16g\t%+016.16g\n",glb_ch_contr[ich_off][0],glb_ch_contr[ich_off][1]);
			    ich_off++;
			  }
			
			//paragraph end
			fprintf(fout,"\n");
		      }
		  }
		//increment the number of stored combo
		ist_combo++;
	      }
      
      //close the output
      if(plan_rank[0]==irank) fclose(fout);
    }
  
  contr_save_time+=take_time();
  
  //free all the vectors
  
  nissa_free(prop_combo);
  
  nissa_free(glb_contr);
  nissa_free(glb_ch_contr);

  nissa_free(ch_colorspinspin);

  ncontr_tot+=(ncontr_3pts+nch_contr_3pts)*nprop_combo;
}

//check all the two points
void check_two_points(int ispec,int ism_lev_so,int ism_lev_se)
{
  double *mass_3pts=mass+start_imass3pts[ispec];
  int nmass_3pts=nmass-start_imass3pts[ispec];

  complex *contr_2pts=nissa_malloc("contr_2pts",ncontr_2pts*glb_size[0],complex);
  
  char path[1024];
  sprintf(path,"%s/2pts_check_sp%d_%02d_%02d",outfolder,ispec,jacobi_niter_so[ism_lev_so],jacobi_niter_se[ism_lev_se]);
  FILE *fout=open_text_file_for_output(path);

  for(int ith2=0;ith2<ntheta;ith2++)
    for(int im2=0;im2<nmass_3pts;im2++)
      {
	int ip2=iprop_of(ith2,im2+start_imass3pts[ispec]);
	contract_with_source(contr_2pts,S1[ip2],op_2pts[1],original_source);
	
	if(rank==0)
	  {
	    fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",
		    mass[imass_spec[ispec]],theta[ith_spec[ispec]],r_spec[ispec],mass_3pts[im2],theta[ith2],r_spec[ispec]);

	    fprintf(fout,"\n");
	    
	    for(int icontr=0;icontr<ncontr_2pts;icontr++)
	      if(op_2pts[0][icontr]==5)
		fprintf(fout," # P5%s\t%+016.16g\t%+016.16g\n",gtag[op_2pts[1][icontr]],(contr_2pts+icontr*glb_size[0])[twall][0],(contr_2pts+icontr*glb_size[0])[twall][1]);
	    fprintf(fout,"\n");
	  }
      }
  
  if(rank==0) fclose(fout);
  
  nissa_free(contr_2pts);
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
	  
	  //loop on spectator
	  for(int ispec=0;ispec<nspec;ispec++)
	    {
	      generate_sequential_source(ispec);
	      
	      //loop on smearing of the sequential prop
	      for(int sm_lev_se=0;sm_lev_se<nsm_lev_se;sm_lev_se++)
		{
		  calculate_S1(ispec,sm_lev_se);
		  check_two_points(ispec,sm_lev_so,sm_lev_se);
		  three_points(ispec,sm_lev_so,sm_lev_se);
		}
	    }	  
	  
	  //loop on the smearing of the sink
	  for(int sm_lev_si=0;sm_lev_si<nsm_lev_si;sm_lev_si++)
	    two_points(sm_lev_so,sm_lev_si);
	}
      
      nanalized_conf++;
      
      enough_time=check_remaining_time();
    }
  
  if(iconf==ngauge_conf) master_printf("Finished all the conf!\n");
  
  tot_time+=take_time();
  close_semileptonic();
  
  return 0;
}
