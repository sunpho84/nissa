#include "nissa.h"

//gauge info
char conf_path[1024];
quad_su3 *conf,*sme_conf;
as2t_su3 *Pmunu;
double kappa;

int nmass,ntheta;
double *mass,*theta;
double put_theta[4],old_theta[4]={0,0,0,0};

//source
int source_coord[4];
spincolor *source;
su3spinspin *original_source;

//smearing parameters
double jacobi_kappa,ape_alpha;
int ape_niter;
int *jacobi_niter_so,nsm_lev_so;
int *jacobi_niter_si,nsm_lev_si;
int *jacobi_niter_se,nsm_lev_se;

//vectors for the spinor data
int npropS0;
su3spinspin **S0[2];
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

//two points chromo contractions
int nch_contr_2pts;
complex *ch_contr_2pts;
int *ch_op1_2pts,*ch_op2_2pts;
su3spinspin *ch_su3spinspin;

//sequential props info
int nspec;
int *imass_spec,*r_spec,*ith_spec;
su3spinspin *sequential_source;

//sequential propagators
int npropS1;
su3spinspin **S1;

//three points contractions
int ncontr_3pts;
complex *contr_3pts;
int *op1_3pts,*op2_3pts;
char outfile_3pts[1024];

//two points chromo contractions
int nch_contr_3pts;
complex *ch_contr_3pts;
int *ch_op1_3pts,*ch_op2_3pts;

//timings
int ninv_tot=0,ncontr_tot=0;
double tot_time=0,inv_time=0,contr_time=0;

//return the position of the propagator of theta and mass
int iprop_of(int itheta,int imass){return itheta*nmass+imass;}

//This function takes care to make the revert on the FIRST spinor, putting the needed gamma5
void meson_two_points(complex *corr,int *list_op1,su3spinspin *s1,int *list_op2,su3spinspin *s2,int ncontr)
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
  trace_g_ccss_dag_g_ccss(corr,t1,s1,t2,s2,ncontr);
}

//This function contract a source with a sequential spinor putting the passed list of operators
void contract_with_source(complex *corr,su3spinspin *S1,int *list_op,su3spinspin *source)
{
  //Temporary vector for the internal matrices
  dirac_matr t1[ncontr_2pts],t2[ncontr_2pts];

  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    { //Init the second 
      t1[icontr]=base_gamma[0];
      t2[icontr]=base_gamma[list_op[icontr]];
    }

  //Call the routine which does the real contraction
  trace_g_ccss_dag_g_ccss(corr,t1,S1,t2,source,ncontr_2pts);
}

//generate the source
void generate_source()
{
  memset(original_source,0,sizeof(su3spinspin)*loc_vol);
  
  int islocal=1,lx[4];
  for(int idir=0;idir<4;idir++)
    {
      lx[idir]=source_coord[idir]-rank_coord[idir]*loc_size[idir];
      islocal&=(lx[idir]>=0);
      islocal&=(lx[idir]<loc_size[idir]);
    }
  
  if(islocal)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	original_source[loclx_of_coord(lx)][ic][ic][id][id][0]=1;
}

//Generate a sequential source for S1
void generate_sequential_source(int ispec)
{
  int r=r_spec[ispec];

  if(rank==0) printf("Creating the sequential source\n");
  for(int ivol=0;ivol<loc_vol;ivol++)
    { //put to zero everything but the slice
      if(glb_coord_of_loclx[ivol][0]!=(source_coord[0]+glb_size[0]/2)%glb_size[0])
	memset(sequential_source[ivol],0,sizeof(su3spinspin));
      else
	{ //avoid to put g5, beacuse commutate with (i+-g5)/sqrt(2) and cancel with those of the QQ
	  memcpy(sequential_source[ivol],S0[r][iprop_of(ith_spec[ispec],imass_spec[ispec])][ivol],sizeof(su3spinspin));
	  for(int c1=0;c1<3;c1++) //rotate as r because it's D^-1
	    for(int c2=0;c2<3;c2++)
	      rotate_spinspin_to_physical_basis(sequential_source[ivol][c1][c2],r,r);
	}
    }
  if(rank==0) printf("Sequential source created\n");
}  

//Parse all the input file
void initialize_semileptonic(char *input_path)
{
  open_input(input_path);

  // 1) Read information about the gauge conf
  
  //Read the volume
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Gauge path
  read_str_str("GaugeConfPath",conf_path,1024);
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the source
  
  //Read the location of the wall
  read_str_int("SourceCoordTXYZ",&(source_coord[0]));
  read_int(&(source_coord[1]));
  read_int(&(source_coord[2]));
  read_int(&(source_coord[3]));

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
  contr_2pts=nissa_malloc("contr_2pts",ncontr_2pts*glb_size[0],complex);
  op1_2pts=nissa_malloc("op1_2pts",ncontr_2pts,int);
  op2_2pts=nissa_malloc("op2_2pts",ncontr_2pts,int);
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1_2pts[icontr]));
      read_int(&(op2_2pts[icontr]));

      if(rank==0 && debug_lvl) printf(" contr.%d %d %d\n",icontr,op1_2pts[icontr],op2_2pts[icontr]);
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

      if(rank==0 && debug_lvl) printf(" ch-contr.%d %d %d\n",icontr,ch_op1_2pts[icontr],ch_op2_2pts[icontr]);
    }

  read_str_str("OutfileTwoPoints",outfile_2pts,1024);

  // 7) three points functions
  
  sequential_source=nissa_malloc("Sequential source",loc_vol,su3spinspin);
  read_str_int("NSpec",&nspec);
  ith_spec=nissa_malloc("ith_spec",nspec,int);
  imass_spec=nissa_malloc("imass_spec",nspec,int);
  r_spec=nissa_malloc("r_spec",nspec,int);
  expect_str("iThetaMassR");
  for(int ispec=0;ispec<nspec;ispec++)
    {
      read_int(&(ith_spec[ispec]));
      read_int(&(imass_spec[ispec]));
      read_int(&(r_spec[ispec]));
      if(rank==0)printf(" spec %d: th=%g, m=%g, r=%d\n",ispec,theta[ith_spec[ispec]],mass[imass_spec[ispec]],r_spec[ispec]);
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

      if(rank==0 && debug_lvl) printf(" contr.%d %d %d\n",icontr,op1_3pts[icontr],op2_3pts[icontr]);
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

      if(rank==0 && debug_lvl) printf(" ch-contr.%d %d %d\n",icontr,ch_op1_3pts[icontr],ch_op2_3pts[icontr]);
    }

  read_str_str("OutfileThreePoints",outfile_3pts,1024);

  close_input();

  ////////////////////////////////////// end of input reading/////////////////////////////////

  //allocate gauge conf, Pmunu and all the needed spincolor and su3spinspin
  conf=nissa_malloc("or_conf",loc_vol+loc_bord+loc_edge,quad_su3);
  sme_conf=nissa_malloc("sm_conf",loc_vol+loc_bord,quad_su3);
  Pmunu=nissa_malloc("Pmunu",loc_vol,as2t_su3);

  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_ildg_gauge_conf(conf,conf_path);
  communicate_lx_gauge_borders(conf);
  communicate_lx_gauge_edges(conf);
  Pmunu_term(Pmunu,conf);
  //prepare the smerded version and calculate plaquette
  ape_smearing(sme_conf,conf,ape_alpha,ape_niter);
  communicate_lx_gauge_borders(conf);
  communicate_lx_gauge_borders(sme_conf);

  double gplaq=global_plaquette_lx_conf(conf);
  if(rank==0) printf("plaq: %.18g\n",gplaq);
  gplaq=global_plaquette_lx_conf(sme_conf);
  if(rank==0) printf("smerded plaq: %.18g\n",gplaq);
    
  //Put the anti-periodic condition on the temporal border
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,1);

  //Allocate all the S0 su3spinspin vectors
  npropS0=ntheta*nmass;
  S0[0]=nissa_malloc("S0[0]",npropS0,su3spinspin*);
  S0[1]=nissa_malloc("S0[1]",npropS0,su3spinspin*);
  for(int iprop=0;iprop<npropS0;iprop++)
    {
      S0[0][iprop]=nissa_malloc("S0[0]",loc_vol,su3spinspin);
      S0[1][iprop]=nissa_malloc("S0[1]",loc_vol,su3spinspin);
    }

  //Allocate nmass spincolors, for the cgmms solutions
  cgmms_solution=nissa_malloc("cgmms_solution",nmass,spincolor*);
  for(int imass=0;imass<nmass;imass++) cgmms_solution[imass]=nissa_malloc("cgmms_solution",loc_vol+loc_bord,spincolor);
  reco_solution[0]=nissa_malloc("reco_solution[0]",loc_vol,spincolor);
  reco_solution[1]=nissa_malloc("reco_solution[1]",loc_vol,spincolor);
  
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+loc_bord,spincolor);
  original_source=nissa_malloc("original_source",loc_vol,su3spinspin);

  //Allocate one su3pinspin for the chromo-contractions
  ch_su3spinspin=nissa_malloc("chromo-su3spinspin",loc_vol,su3spinspin);

  //Allocate all the S1 su3spinspin vectors
  npropS1=ntheta*nmass;
  S1=nissa_malloc("S1",loc_vol,su3spinspin*);
  for(int iprop=0;iprop<npropS0;iprop++) S1[iprop]=nissa_malloc("S1[i]",loc_vol,su3spinspin);
}

//Finalization
void close_semileptonic()
{
  if(rank==0)
    {
      printf("\n");
      printf("Total time: %g, of which:\n",tot_time);
      printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
      printf(" - %02.2f%s to perform %d contr. (%2.2gs avg)\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
    }

  nissa_free(Pmunu);nissa_free(conf);nissa_free(sme_conf);
  for(int iprop=0;iprop<npropS0;iprop++){nissa_free(S0[0][iprop]);nissa_free(S0[1][iprop]);nissa_free(S1[iprop]);}
  nissa_free(S0[0]);nissa_free(S0[1]);nissa_free(S1);
  nissa_free(reco_solution[0]);nissa_free(reco_solution[1]);
  nissa_free(ch_su3spinspin);nissa_free(sequential_source);
  nissa_free(contr_2pts);nissa_free(ch_contr_2pts);
  nissa_free(contr_3pts);nissa_free(ch_contr_3pts);
  nissa_free(op1_2pts);nissa_free(op2_2pts);
  nissa_free(op1_3pts);nissa_free(op2_3pts);
  nissa_free(ch_op1_2pts);nissa_free(ch_op2_2pts);
  nissa_free(ch_op1_3pts);nissa_free(ch_op2_3pts);
  nissa_free(ith_spec);nissa_free(r_spec);nissa_free(imass_spec);
  for(int imass=0;imass<nmass;imass++) nissa_free(cgmms_solution[imass]);
  nissa_free(cgmms_solution);
  nissa_free(source);nissa_free(original_source);
  close_nissa();
}

//smear addditivily a su3spinspin
void smear_additive_su3spinspin(su3spinspin *out,su3spinspin *in,int ism_lev,int *jacobi_niter)
{
  spincolor *temp=nissa_malloc("temp",loc_vol+loc_bord,spincolor);
  
  int nsme=jacobi_niter[ism_lev];
  if(ism_lev>0) nsme-=jacobi_niter[ism_lev-1];

  for(int ic=0;ic<3;ic++)
    for(int id=0;id<4;id++)
      {
	for(int ivol=0;ivol<loc_vol;ivol++)
	  get_spincolor_from_su3spinspin(temp[ivol],in[ivol],id,ic);
	
	jacobi_smearing(temp,temp,sme_conf,jacobi_kappa,nsme);
	
	for(int ivol=0;ivol<loc_vol;ivol++)
	  put_spincolor_into_su3spinspin(out[ivol],temp[ivol],id,ic);
      }
  
  nissa_free(temp);
}

//calculate the standard propagators
void calculate_S0(int ism_lev_so)
{
  smear_additive_su3spinspin(original_source,original_source,ism_lev_so,jacobi_niter_so);
  
  for(int ic=0;ic<3;ic++)
    for(int id=0;id<4;id++)
      { //loop over the source dirac index
	for(int ivol=0;ivol<loc_vol;ivol++)
	  {
	    get_spincolor_from_su3spinspin(source[ivol],original_source[ivol],id,ic);
	    //put the g5
	    for(int id1=2;id1<4;id1++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic][ri]*=-1;
	  }
	
	communicate_lx_spincolor_borders(source);
	
	for(int itheta=0;itheta<ntheta;itheta++)
	  { //adapt the border condition
	    put_theta[1]=put_theta[2]=put_theta[3]=theta[itheta];
	    adapt_theta(conf,old_theta,put_theta,1,1);
	    
	    double part_time=-take_time();
	    inv_tmQ2_cgmms(cgmms_solution,source,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	    part_time+=take_time();ninv_tot++;inv_time+=part_time;
	    if(rank==0) printf("Finished the inversion of S0 theta %d, dirac index %d in %g sec\n",itheta,id,part_time);
	    
	    for(int imass=0;imass<nmass;imass++)
	      { //reconstruct the doublet
		reconstruct_tm_doublet(reco_solution[0],reco_solution[1],cgmms_solution[imass],conf,kappa,mass[imass]);
		if(rank==0) printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
		for(int r=0;r<2;r++) //convert the id-th spincolor into the su3spinspin
		  for(int i=0;i<loc_vol;i++)
		    put_spincolor_into_su3spinspin(S0[r][iprop_of(itheta,imass)][i],reco_solution[r][i],id,ic);
	      }
	  }
      }
  
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
      rotate_vol_su3spinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
}

//calculate the sequential propagators
void calculate_S1(int ispec,int ism_lev_se)
{
  smear_additive_su3spinspin(sequential_source,sequential_source,ism_lev_se,jacobi_niter_se);

  for(int ic=0;ic<3;ic++)
    for(int id=0;id<4;id++)
      { 
	for(int ivol=0;ivol<loc_vol;ivol++) //avoid the g5 insertion
	  get_spincolor_from_su3spinspin(source[ivol],sequential_source[ivol],id,ic);
	communicate_lx_spincolor_borders(source);
	
	for(int itheta=0;itheta<ntheta;itheta++)
	  { //adapt the border condition
	    put_theta[1]=put_theta[2]=put_theta[3]=theta[itheta];
	    adapt_theta(conf,old_theta,put_theta,1,1);
	    
	    double part_time=-take_time();
	    inv_tmQ2_cgmms(cgmms_solution,source,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	    part_time+=take_time();ninv_tot++;inv_time+=part_time;
	    if(rank==0) printf("Finished the inversion of S1 theta %d, seq sme lev %d, dirac index %d in %g sec\n",itheta,ism_lev_se,id,part_time);
	    
	    for(int imass=0;imass<nmass;imass++)
	      { //reconstruct the doublet: r(S1)=!r(spec), so we have to multiply by Q+ if r(spec)==1 and Q- if 0
		double reco_mass=-mass[imass];
		if(r_spec[ispec]==1) reco_mass=-reco_mass;
		//use reco_solution[0] as temporary storage
		apply_tmQ(reco_solution[0],cgmms_solution[imass],conf,kappa,reco_mass);
		if(rank==0) printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
		for(int i=0;i<loc_vol;i++) put_spincolor_into_su3spinspin(S1[iprop_of(itheta,imass)][i],reco_solution[0][i],id,ic);
	      }
	  }
      }
  
  //put the (1+-ig5)/sqrt(2) factor. On the source rotate as r_spec, on the sink as !r_spec
  for(int ipropS1=0;ipropS1<npropS1;ipropS1++) //but, being D^-1, everything is swapped
    rotate_vol_su3spinspin_to_physical_basis(S1[ipropS1],!(r_spec[ispec]),(r_spec[ispec]));
}

//Calculate and print to file the 2pts
void calculate_all_2pts(int ism_lev_so,int ism_lev_si)
{
  for(int r=0;r<2;r++)
    for(int iprop=0;iprop<npropS0;iprop++)
      smear_additive_su3spinspin(S0[r][iprop],S0[r][iprop],ism_lev_si,jacobi_niter_si);
  
  char path[1024];
  sprintf(path,"%s_%02d_%02d",outfile_2pts,jacobi_niter_so[ism_lev_so],jacobi_niter_si[ism_lev_si]);
  FILE *fout=open_text_file_for_output(path);

  contr_time-=take_time();

  for(int ispec=0;ispec<nspec;ispec++)
    {
      int ith1=ith_spec[ispec];
      for(int ith2=0;ith2<ntheta;ith2++)
	for(int im2=0;im2<nmass;im2++)
	  for(int r2=0;r2<2;r2++)
	    for(int im1=0;im1<nmass;im1++)
	      for(int r1=0;r1<2;r1++)
		{
		  int ip1=iprop_of(ith1,im1),ip2=iprop_of(ith2,im2);
		  
		  if(rank==0)
		    {
		      fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d",mass[im1],theta[ith1],r1,mass[im2],theta[ith2],r2);
		      fprintf(fout," smear_source=%d smear_sink=%d\n",jacobi_niter_so[ism_lev_so],jacobi_niter_si[ism_lev_si]);
		    }
		  
		  meson_two_points(contr_2pts,op1_2pts,S0[r1][ip1],op2_2pts,S0[r2][ip2],ncontr_2pts);
		  ncontr_tot+=ncontr_2pts;
		  print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,source_coord[0],"",1);
		  
		  if(nch_contr_2pts>0)
		    {
		      unsafe_apply_chromo_operator_to_su3spinspin(ch_su3spinspin,Pmunu,S0[r2][ip2]);
		      meson_two_points(ch_contr_2pts,ch_op1_2pts,S0[r1][ip1],ch_op2_2pts,ch_su3spinspin,nch_contr_2pts);
		      ncontr_tot+=nch_contr_2pts;
		      print_contractions_to_file(fout,nch_contr_2pts,ch_op1_2pts,ch_op2_2pts,ch_contr_2pts,source_coord[0],"CHROMO-",1);
		    }
		  if(rank==0) fprintf(fout,"\n");
		}
    }
  contr_time+=take_time();

  if(rank==0) fclose(fout);
}

//Calculate and print to file the 3pts
void calculate_all_3pts(int ispec,int ism_lev_so,int ism_lev_se)
{
  char path[1024];
  
  sprintf(path,"%s_sp%d_%02d_%02d",outfile_3pts,ispec,jacobi_niter_so[ism_lev_so],jacobi_niter_se[ism_lev_se]);

  FILE *fout=open_text_file_for_output(path);
 
  contr_time-=take_time();
  
  int r1=r_spec[ispec];
  int r2=!r1;
  for(int ith2=0;ith2<ntheta;ith2++)
    for(int im2=0;im2<nmass;im2++)
      for(int ith1=0;ith1<ntheta;ith1++)
	for(int im1=0;im1<nmass;im1++)
	  {
	    int ip1=iprop_of(ith1,im1),ip2=iprop_of(ith2,im2);
	    
	    if(rank==0)
	      {
		fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d,",mass[im1],theta[ith1],r1,mass[im2],theta[ith2],r2);
		fprintf(fout," smear_source=%d smear_seq=%d\n",jacobi_niter_so[ism_lev_so],jacobi_niter_se[ism_lev_se]);
	      }
	    
	    meson_two_points(contr_3pts,op1_3pts,S0[r1][ip1],op2_3pts,S1[ip2],ncontr_3pts);
	    ncontr_tot+=ncontr_3pts;
	    print_contractions_to_file(fout,ncontr_3pts,op1_3pts,op2_3pts,contr_3pts,source_coord[0],"",1);
	    
	    if(nch_contr_3pts>0)
	      {
		unsafe_apply_chromo_operator_to_su3spinspin(ch_su3spinspin,Pmunu,S1[ip2]);
		meson_two_points(ch_contr_3pts,ch_op1_3pts,S0[r1][ip1],ch_op2_3pts,ch_su3spinspin,nch_contr_3pts);
		ncontr_tot+=nch_contr_3pts;
		print_contractions_to_file(fout,nch_contr_3pts,ch_op1_3pts,ch_op2_3pts,ch_contr_3pts,source_coord[0],"CHROMO-",1);
	      }
	    if(rank==0) fprintf(fout,"\n");
	  }

  contr_time+=take_time();
  
  if(rank==0) fclose(fout);
}

//check all the two points
void check_two_points(int ispec,int ism_lev_so,int ism_lev_se)
{
  char path[1024];
  sprintf(path,"2pts_check_sp%d_%02d_%02d",ispec,jacobi_niter_so[ism_lev_so],jacobi_niter_se[ism_lev_se]);
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
		fprintf(fout," # P5%s\t%+016.16g\t%+016.16g\n",gtag[op2_2pts[icontr]],(contr_2pts+icontr*glb_size[0])[source_coord[0]][0],(contr_2pts+icontr*glb_size[0])[source_coord[0]][1]);
	    fprintf(fout,"\n");
	  }
      }
  
  if(rank==0) fclose(fout);
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();

  if(narg<2) crash("Use: %s input_file",arg[0]);

  tot_time-=take_time();
  initialize_semileptonic(arg[1]);
  
  generate_source();
  
  for(int sm_lev_so=0;sm_lev_so<nsm_lev_so;sm_lev_so++)
    {
      calculate_S0(sm_lev_so);
      
      for(int ispec=0;ispec<nspec;ispec++)
	{
	  generate_sequential_source(ispec);

	  for(int sm_lev_se=0;sm_lev_se<nsm_lev_se;sm_lev_se++)
	    {
	      calculate_S1(ispec,sm_lev_se);
	      check_two_points(ispec,sm_lev_so,sm_lev_se);
	      calculate_all_3pts(ispec,sm_lev_so,sm_lev_se);
	    }
	}	  
      
      for(int sm_lev_si=0;sm_lev_si<nsm_lev_si;sm_lev_si++) calculate_all_2pts(sm_lev_so,sm_lev_si);
    }

  tot_time+=take_time();
  close_semileptonic();
  
  return 0;
}
