#include "appretto.h"

char conf_path[1024];
quad_su3 *conf;
as2t_su3 *Pmunu;
double kappa;

int nmass,ntheta;
double *mass,*theta;
double put_theta[4],old_theta[4]={0,0,0,0};

//source data
int seed,noise_type,twall;
double rad2=1.414213562373095048801688724209;
spincolor *source;
colorspinspin *original_source;

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

//two points chromo contractions
int nch_contr_2pts;
complex *ch_contr_2pts;
int *ch_op1_2pts,*ch_op2_2pts;
colorspinspin *ch_colorspinspin;

//sequential props info
int nspec;
int *imass_spec,*r_spec,*ith_spec;
colorspinspin *sequential_source;

//sequential propagators
int npropS1;
colorspinspin **S1;

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

//Generate a sequential source for S1
void generate_sequential_source(int ispec)
{
  int r=r_spec[ispec];

  if(rank==0) printf("Creating the sequential source\n");
  for(int ivol=0;ivol<loc_vol;ivol++)
    { //put to zero everything but the slice
      if(glb_coord_of_loclx[ivol][0]!=(twall+glb_size[0]/2)%glb_size[0])
	memset(sequential_source[ivol],0,sizeof(colorspinspin));
      else
	{ //avoid to put g5, beacuse commutate with (i+-g5)/sqrt(2) and cancel with those of the QQ
	  memcpy(sequential_source[ivol],S0[r][iprop_of(ith_spec[ispec],imass_spec[ispec])][ivol],sizeof(colorspinspin));
	  for(int c=0;c<3;c++) //rotate as r because it's D^-1
	    rotate_spinspin_to_physical_basis(sequential_source[ivol][c],r,r);
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
  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  //Init the MPI grid 
  init_grid(); 
  //Gauge path
  read_str_str("GaugeConfPath",conf_path,1024);
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the source
  
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  init_random(seed);
  //Read the location of the wall
  read_str_int("TWall",&twall);
  //Read the noise type
  read_str_int("NoiseType",&noise_type);

  // 3) Read list of masses and of thetas

  read_list_of_doubles("NMass",&nmass,&mass);
  read_list_of_doubles("NTheta",&ntheta,&theta);

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
  contr_2pts=(complex*)malloc(sizeof(complex)*ncontr_2pts*glb_size[0]); 
  op1_2pts=(int*)malloc(sizeof(int)*ncontr_2pts);
  op2_2pts=(int*)malloc(sizeof(int)*ncontr_2pts);
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1_2pts[icontr]));
      read_int(&(op2_2pts[icontr]));

      if(rank==0 && debug) printf(" contr.%d %d %d\n",icontr,op1_2pts[icontr],op2_2pts[icontr]);
    }
  
  read_str_int("NChromoContrTwoPoints",&nch_contr_2pts);
  ch_contr_2pts=(complex*)malloc(sizeof(complex)*nch_contr_2pts*glb_size[0]);
  ch_op1_2pts=(int*)malloc(sizeof(int)*nch_contr_2pts);
  ch_op2_2pts=(int*)malloc(sizeof(int)*nch_contr_2pts);
  for(int icontr=0;icontr<nch_contr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(ch_op1_2pts[icontr]));
      read_int(&(ch_op2_2pts[icontr]));

      if(rank==0 && debug) printf(" ch-contr.%d %d %d\n",icontr,ch_op1_2pts[icontr],ch_op2_2pts[icontr]);
    }

  read_str_str("OutfileTwoPoints",outfile_2pts,1024);

  // 6) three points functions
  sequential_source=appretto_malloc("Sequential source",loc_vol,colorspinspin);
  read_str_int("NSpec",&nspec);
  ith_spec=(int*)malloc(nspec*sizeof(int));
  imass_spec=(int*)malloc(nspec*sizeof(int));
  r_spec=(int*)malloc(nspec*sizeof(int));
  expect_str("iThetaMassR");
  for(int ispec=0;ispec<nspec;ispec++)
    {
      read_int(&(ith_spec[ispec]));
      read_int(&(imass_spec[ispec]));
      read_int(&(r_spec[ispec]));
      if(rank==0)printf(" spec %d: th=%g, m=%g, r=%d\n",ispec,theta[ith_spec[ispec]],mass[imass_spec[ispec]],r_spec[ispec]);
    }
  read_str_int("NContrThreePoints",&ncontr_3pts);
  contr_3pts=(complex*)malloc(sizeof(complex)*ncontr_3pts*glb_size[0]); 
  op1_3pts=(int*)malloc(sizeof(int)*ncontr_3pts);
  op2_3pts=(int*)malloc(sizeof(int)*ncontr_3pts);
  for(int icontr=0;icontr<ncontr_3pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1_3pts[icontr]));
      read_int(&(op2_3pts[icontr]));

      if(rank==0 && debug) printf(" contr.%d %d %d\n",icontr,op1_3pts[icontr],op2_3pts[icontr]);
    }

  read_str_int("NChromoContrThreePoints",&nch_contr_3pts);
  ch_contr_3pts=(complex*)malloc(sizeof(complex)*nch_contr_3pts*glb_size[0]);
  ch_op1_3pts=(int*)malloc(sizeof(int)*nch_contr_3pts);
  ch_op2_3pts=(int*)malloc(sizeof(int)*nch_contr_3pts);
  for(int icontr=0;icontr<nch_contr_3pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(ch_op1_3pts[icontr]));
      read_int(&(ch_op2_3pts[icontr]));

      if(rank==0 && debug) printf(" ch-contr.%d %d %d\n",icontr,ch_op1_3pts[icontr],ch_op2_3pts[icontr]);
    }

  read_str_str("OutfileThreePoints",outfile_3pts,1024);

  close_input();

  ////////////////////////////////////// end of input reading/////////////////////////////////

  //allocate gauge conf, Pmunu and all the needed spincolor and colorspinspin
  conf=appretto_malloc("conf",loc_vol+loc_bord+loc_edge,quad_su3);
  Pmunu=appretto_malloc("Pmunu",loc_vol,as2t_su3);

  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_local_gauge_conf(conf,conf_path);
  communicate_gauge_borders(conf);
  communicate_gauge_edges(conf);
  
  double gplaq=global_plaquette(conf);
  if(rank==0) printf("plaq: %.18g\n",gplaq);
  
  Pmunu_term(Pmunu,conf);
  
  //Put the anti-periodic condition on the temporal border
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,1);

  //Allocate all the S0 colorspinspin vectors
  npropS0=ntheta*nmass;
  S0[0]=(colorspinspin**)malloc(sizeof(colorspinspin*)*npropS0);
  S0[1]=(colorspinspin**)malloc(sizeof(colorspinspin*)*npropS0);
  for(int iprop=0;iprop<npropS0;iprop++)
    {
      S0[0][iprop]=appretto_malloc("S0[0]",loc_vol,colorspinspin);
      S0[1][iprop]=appretto_malloc("S0[1]",loc_vol,colorspinspin);
    }

  //Allocate nmass spincolors, for the cgmms solutions
  cgmms_solution=(spincolor**)malloc(sizeof(spincolor*)*nmass);
  for(int imass=0;imass<nmass;imass++) cgmms_solution[imass]=appretto_malloc("cgmms_solution",loc_vol+loc_bord,spincolor);
  reco_solution[0]=appretto_malloc("reco_solution[0]",loc_vol,spincolor);
  reco_solution[1]=appretto_malloc("reco_solution[1]",loc_vol,spincolor);
  
  //Allocate one spincolor for the source
  source=appretto_malloc("source",loc_vol+loc_bord,spincolor);
  original_source=appretto_malloc("original_source",loc_vol,colorspinspin);

  //Allocate one colorspinspin for the chromo-contractions
  ch_colorspinspin=appretto_malloc("chromo-colorspinspin",loc_vol,colorspinspin);

  //Allocate all the S1 colorspinspin vectors
  npropS1=ntheta*nmass;
  S1=(colorspinspin**)malloc(sizeof(colorspinspin*)*loc_vol);
  for(int iprop=0;iprop<npropS0;iprop++) S1[iprop]=appretto_malloc("S1[i]",loc_vol,colorspinspin);
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

  appretto_free(Pmunu);appretto_free(conf);
  for(int iprop=0;iprop<npropS0;iprop++){appretto_free(S0[0][iprop]);appretto_free(S0[1][iprop]);appretto_free(S1[iprop]);}
  appretto_free(reco_solution[0]);appretto_free(reco_solution[1]);
  appretto_free(ch_colorspinspin);appretto_free(sequential_source);
  for(int imass=0;imass<nmass;imass++) appretto_free(cgmms_solution[imass]);
  appretto_free(source);appretto_free(original_source);
  close_appretto();
}

//calculate the standard propagators
void calculate_S0()
{
  for(int id=0;id<4;id++)
    { //loop over the source dirac index
      for(int ivol=0;ivol<loc_vol;ivol++)
	{
	  get_spincolor_from_colorspinspin(source[ivol],original_source[ivol],id);
	  //put the g5
	  for(int id1=2;id1<4;id1++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic][ri]*=-1;
	}
      
      communicate_lx_spincolor_borders(source);
      
      for(int itheta=0;itheta<ntheta;itheta++)
	{ //adapt the border condition
	  put_theta[1]=put_theta[2]=put_theta[3]=theta[itheta];
	  adapt_theta(conf,old_theta,put_theta,1,1);
	  
	  double part_time=-take_time();
	  inv_Q2_cgmms(cgmms_solution,source,NULL,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	  part_time+=take_time();ninv_tot++;inv_time+=part_time;
	  if(rank==0) printf("Finished the inversion of S0 theta %d, dirac index %d in %g sec\n",itheta,id,part_time);

	  for(int imass=0;imass<nmass;imass++)
	    { //reconstruct the doublet
	      reconstruct_doublet(reco_solution[0],reco_solution[1],cgmms_solution[imass],conf,kappa,mass[imass]);
	      if(rank==0) printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	      for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
		for(int i=0;i<loc_vol;i++)
		  put_spincolor_into_colorspinspin(S0[r][iprop_of(itheta,imass)][i],reco_solution[r][i],id);
	    }
	}
    }
  
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
      rotate_vol_colorspinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
}

//calculate the sequential propagators
void calculate_S1(int ispec)
{
  for(int id=0;id<4;id++)
    { 
      for(int ivol=0;ivol<loc_vol;ivol++) //avoid the g5 insertion
	get_spincolor_from_colorspinspin(source[ivol],sequential_source[ivol],id);
      communicate_lx_spincolor_borders(source);

      for(int itheta=0;itheta<ntheta;itheta++)
	{ //adapt the border condition
	  put_theta[1]=put_theta[2]=put_theta[3]=theta[itheta];
	  adapt_theta(conf,old_theta,put_theta,1,1);

	  double part_time=-take_time();
	  inv_Q2_cgmms(cgmms_solution,source,NULL,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	  part_time+=take_time();ninv_tot++;inv_time+=part_time;
	  if(rank==0) printf("Finished the inversion of S1 theta %d, dirac index %d in %g sec\n",itheta,id,part_time);
	  
	  for(int imass=0;imass<nmass;imass++)
	    { //reconstruct the doublet: r(S1)=!r(spec), so we have to multiply by Q+ if r(spec)==1 and Q- if 0
	      double reco_mass=-mass[imass];
	      if(r_spec[ispec]==1) reco_mass=-reco_mass;
	      //use reco_solution[0] as temporary storage
	      apply_Q(reco_solution[0],cgmms_solution[imass],conf,kappa,reco_mass);
	      if(rank==0) printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	      for(int i=0;i<loc_vol;i++) put_spincolor_into_colorspinspin(S1[iprop_of(itheta,imass)][i],reco_solution[0][i],id);
	    }
	}
    }

  //put the (1+-ig5)/sqrt(2) factor. On the source rotate as r_spec, on the sink as !r_spec
  for(int ipropS1=0;ipropS1<npropS1;ipropS1++) //but, being D^-1, everything is swapped
    rotate_vol_colorspinspin_to_physical_basis(S1[ipropS1],!(r_spec[ispec]),(r_spec[ispec]));
}

//Calculate and print to file the 2pts
void calculate_all_2pts()
{
  FILE *fout=open_text_file_for_output(outfile_2pts);

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
		  
		  if(rank==0) fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",
				      mass[im1],theta[ith1],r1,mass[im2],theta[ith2],r2);
		  
		  meson_two_points(contr_2pts,op1_2pts,S0[r1][ip1],op2_2pts,S0[r2][ip2],ncontr_2pts);
		  ncontr_tot+=ncontr_2pts;
		  print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,twall,"");
		  
		  if(nch_contr_2pts>0)
		    {
		      unsafe_apply_chromo_operator_to_colorspinspin(ch_colorspinspin,Pmunu,S0[r2][ip2]);
		      meson_two_points(ch_contr_2pts,ch_op1_2pts,S0[r1][ip1],ch_op2_2pts,ch_colorspinspin,nch_contr_2pts);
		      ncontr_tot+=nch_contr_2pts;
		      print_contractions_to_file(fout,nch_contr_2pts,ch_op1_2pts,ch_op2_2pts,ch_contr_2pts,twall,"CHROMO-");
		    }
		  if(rank==0) fprintf(fout,"\n");
		}
    }
  contr_time+=take_time();

  if(rank==0) fclose(fout);
}

//Calculate and print to file the 3pts
void calculate_all_3pts(int ispec)
{ //generate a different filename for each spectator
  char outfile_3pts_spec[1024];
  if(nspec>1) sprintf(outfile_3pts_spec,"%s_%d",outfile_3pts,ispec);
  else memcpy(outfile_3pts_spec,outfile_3pts,1024);

  FILE *fout=open_text_file_for_output(outfile_3pts_spec);
 
  contr_time-=take_time();
  
  int r1=r_spec[ispec];
  int r2=!r1;
  for(int ith2=0;ith2<ntheta;ith2++)
    for(int im2=0;im2<nmass;im2++)
      for(int ith1=0;ith1<ntheta;ith1++)
	for(int im1=0;im1<nmass;im1++)
	  {
	    int ip1=iprop_of(ith1,im1),ip2=iprop_of(ith2,im2);
	    
	    if(rank==0) fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",
				mass[im1],theta[ith1],r1,mass[im2],theta[ith2],r2);

	    meson_two_points(contr_3pts,op1_3pts,S0[r1][ip1],op2_3pts,S1[ip2],ncontr_3pts);
	    ncontr_tot+=ncontr_3pts;
	    print_contractions_to_file(fout,ncontr_3pts,op1_3pts,op2_3pts,contr_3pts,twall,"");
	    
	    if(nch_contr_3pts>0)
	      {
		unsafe_apply_chromo_operator_to_colorspinspin(ch_colorspinspin,Pmunu,S1[ip2]);
		meson_two_points(ch_contr_3pts,ch_op1_3pts,S0[r1][ip1],ch_op2_3pts,ch_colorspinspin,nch_contr_3pts);
		ncontr_tot+=nch_contr_3pts;
		print_contractions_to_file(fout,nch_contr_3pts,ch_op1_3pts,ch_op2_3pts,ch_contr_3pts,twall,"CHROMO-");
	      }
	    if(rank==0) fprintf(fout,"\n");
	  }

  contr_time+=take_time();
  
  if(rank==0) fclose(fout);
}

//check all the two points
void check_two_points(int ispec)
{
  FILE *fout=open_text_file_for_output("2pts_check");
  int spat_vol=glb_size[1]*glb_size[2]*glb_size[3];

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
		fprintf(fout," # P5%s\t%+016.16g\t%+016.16g\n",gtag[op2_2pts[icontr]],(contr_2pts+icontr*glb_size[0])[twall][0]/spat_vol,(contr_2pts+icontr*glb_size[0])[twall][1]/spat_vol);
	    fprintf(fout,"\n");
	  }
      }
  
  if(rank==0) fclose(fout);
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

  tot_time-=take_time();
  initialize_semileptonic(arg[1]);
  
  generate_source();
  calculate_S0();
  calculate_all_2pts();

  if(rank==0) printf("Ok, the 2points part is finished\n");
  
  for(int ispec=0;ispec<nspec;ispec++)
    {
      generate_sequential_source(ispec);
      calculate_S1(ispec);
      check_two_points(ispec);
      calculate_all_3pts(ispec);
    }
  
  tot_time+=take_time();
  close_semileptonic();
  
  return 0;
}
