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
spincolor **cgm_solution,*temp_vec[2];

//cgm inverter parameters
double *stopping_residues;
int niter_max=100000;

// derivative variables
spincolor *der_source;
su3spinspin **S0_der[2][3][3]; //propagateur "derive" su3spinspin
spincolor *temp_der_vec[2];
spincolor **cgm_der_solution;

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

//Generate a sequential source for S1
void generate_sequential_source(int ispec)
{
  int r=r_spec[ispec];

  master_printf("Creating the sequential source\n");
  NISSA_LOC_VOL_LOOP(ivol)
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
  master_printf("Sequential source created\n");
  
  set_borders_invalid(sequential_source);
}  

//Parse all the input file
void initialize_semileptonic(char *input_path)
{
  open_input(input_path);

  // 1) Read information about the gauge conf
  int L,T;
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

      if(nissa_verbosity) master_printf(" contr.%d %d %d\n",icontr,op1_2pts[icontr],op2_2pts[icontr]);
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

      if(nissa_verbosity) master_printf(" ch-contr.%d %d %d\n",icontr,ch_op1_2pts[icontr],ch_op2_2pts[icontr]);
    }

  read_str_str("OutfileTwoPoints",outfile_2pts,1024);

  // 6) three points functions
  
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

      if(nissa_verbosity) master_printf(" contr.%d %d %d\n",icontr,op1_3pts[icontr],op2_3pts[icontr]);
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

      if(nissa_verbosity) master_printf(" ch-contr.%d %d %d\n",icontr,ch_op1_3pts[icontr],ch_op2_3pts[icontr]);
    }

  read_str_str("OutfileThreePoints",outfile_3pts,1024);

  close_input();

  ////////////////////////////////////// end of input reading/////////////////////////////////

  //allocate gauge conf, Pmunu and all the needed spincolor and su3spinspin
  conf=nissa_malloc("or_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  sme_conf=nissa_malloc("sm_conf",loc_vol+bord_vol,quad_su3);
  Pmunu=nissa_malloc("Pmunu",loc_vol,as2t_su3);

  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_ildg_gauge_conf(conf,conf_path);
  Pmunu_term(Pmunu,conf);
  //prepare the smeared version and calculate plaquette
  ape_spatial_smear_conf(sme_conf,conf,ape_alpha,ape_niter);

  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  master_printf("smeared plaq: %.18g\n",global_plaquette_lx_conf(sme_conf));
    
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
	
  //Allocate all the S0_der su3spinspin vectors
  for(int i=0;i<3;i++)
    {
      for(int j=0;j<3;j++) 
	{
	  S0_der[0][i][j]=nissa_malloc("S0_der[0]",npropS0,su3spinspin*);
	  S0_der[1][i][j]=nissa_malloc("S0_der[1]",npropS0,su3spinspin*);
	  
	  for(int iprop=0;iprop<npropS0;iprop++)
	    {
	      S0_der[0][i][j][iprop]=nissa_malloc("S0_der[0]",loc_vol,su3spinspin);
	      S0_der[1][i][j][iprop]=nissa_malloc("S0_der[1]",loc_vol,su3spinspin);							
	    }
	}
      
    }
  
  //Allocate nmass spincolors, for the cgm solutions
  cgm_solution=nissa_malloc("cgm_solution",nmass,spincolor*);
  for(int imass=0;imass<nmass;imass++) cgm_solution[imass]=nissa_malloc("cgm_solution",loc_vol+bord_vol,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol+bord_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol+bord_vol,spincolor);
  
  //Allocate nmass spincolors, for the cgm "derivative" solutions
  cgm_der_solution=nissa_malloc("cgm_solution",nmass,spincolor*);
  for(int imass=0;imass<nmass;imass++) cgm_der_solution[imass]=nissa_malloc("cgm_der_solution",loc_vol+bord_vol,spincolor);
  temp_der_vec[0]=nissa_malloc("temp_der_vec[0]",loc_vol+bord_vol,spincolor);
  temp_der_vec[1]=nissa_malloc("temp_der_vec[1]",loc_vol+bord_vol,spincolor);
  
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  original_source=nissa_malloc("original_source",loc_vol,su3spinspin);
  der_source=nissa_malloc("der_source",loc_vol+bord_vol,spincolor);

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
  master_printf("\n");
  master_printf("Total time: %g, of which:\n",tot_time);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d contr. (%2.2gs avg)\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);

  nissa_free(Pmunu);nissa_free(conf);nissa_free(sme_conf);
  for(int iprop=0;iprop<npropS0;iprop++){nissa_free(S0[0][iprop]);nissa_free(S0[1][iprop]);nissa_free(S1[iprop]);}
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++) 
      for(int iprop=0;iprop<npropS0;iprop++)
	{nissa_free(S0_der[0][i][j][iprop]);nissa_free(S0_der[1][i][j][iprop]);}
  nissa_free(S0[0]);nissa_free(S0[1]);nissa_free(S1);
  for(int i=0;i<3;i++)
	for(int j=0;j<3;j++) 
	  {nissa_free(S0_der[0][i][j]); nissa_free(S0_der[1][i][j]);}
  nissa_free(temp_der_vec[0]);nissa_free(temp_der_vec[1]);
  nissa_free(temp_vec[0]);nissa_free(temp_vec[1]);
  nissa_free(ch_su3spinspin);nissa_free(sequential_source);
  nissa_free(contr_2pts);nissa_free(ch_contr_2pts);
  nissa_free(contr_3pts);nissa_free(ch_contr_3pts);
  nissa_free(op1_2pts);nissa_free(op2_2pts);
  nissa_free(op1_3pts);nissa_free(op2_3pts);
  nissa_free(ch_op1_2pts);nissa_free(ch_op2_2pts);
  nissa_free(ch_op1_3pts);nissa_free(ch_op2_3pts);
  nissa_free(ith_spec);nissa_free(r_spec);nissa_free(imass_spec);
  for(int imass=0;imass<nmass;imass++) {nissa_free(cgm_solution[imass]); nissa_free(cgm_der_solution[imass]);}
  nissa_free(cgm_solution);nissa_free(cgm_der_solution);
  nissa_free(source);nissa_free(original_source);nissa_free(der_source);
  close_nissa();
}

//smear addditivily a su3spinspin
void smear_additive_su3spinspin(su3spinspin *out,su3spinspin *in,int ism_lev,int *jacobi_niter)
{
  spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  
  int nsme=jacobi_niter[ism_lev];
  if(ism_lev>0) nsme-=jacobi_niter[ism_lev-1];

  for(int ic=0;ic<3;ic++)
    for(int id=0;id<4;id++)
      {
	NISSA_LOC_VOL_LOOP(ivol)
	  get_spincolor_from_su3spinspin(temp[ivol],in[ivol],id,ic);
	
	jacobi_smearing(temp,temp,sme_conf,jacobi_kappa,nsme);
	
	NISSA_LOC_VOL_LOOP(ivol)
	  put_spincolor_into_su3spinspin(out[ivol],temp[ivol],id,ic);
      }
  
  nissa_free(temp);
}

//Apply the covariant derivative having direction mu on a spincolor
void apply_nabla_i_spincolor(spincolor *out,spincolor *in,quad_su3 *conf,int mu)
{
  memset(out,0,loc_vol*sizeof(spincolor));
  communicate_lx_spincolor_borders(in);
  
  NISSA_LOC_VOL_LOOP(ix)
    {
      int Xup,Xdw;
      Xup=loclx_neighup[ix][mu];
      Xdw=loclx_neighdw[ix][mu];
      
      unsafe_su3_prod_spincolor(out[ix],conf[ix][mu],in[Xup]);
      unsafe_su3_dag_subt_the_prod_spincolor(out[ix],conf[Xdw][mu],in[Xdw]);
    }
  
  set_borders_invalid(out);
}

/*
void apply_nabla_dag_i_spincolor(spincolor *out,spincolor *in,quad_su3 *conf,int mu)
{
  memset(out,0,loc_vol*sizeof(spincolor));
  communicate_lx_spincolor_borders(in);
  
  NISSA_LOC_VOL_LOOP(ix)
    {
      int Xup,Xdw;
      Xup=loclx_neighup[ix][mu];
      Xdw=loclx_neighdw[ix][mu];
      
      safe_su3_dag_prod_spincolor(out[ix], conf[ix][mu],in[Xup]);
      su3_subt_the_prod_spincolor(out[ix],conf[Xdw][mu],in[Xdw]); 
    }
  
  set_borders_invalid(out);
}
*/

//calculate the standard propagators
void calculate_S0(int ism_lev_so)
{
  smear_additive_su3spinspin(original_source,original_source,ism_lev_so,jacobi_niter_so);
  
  for(int ic=0;ic<3;ic++)
    for(int id=0;id<4;id++)
      { //loop over the source dirac index
	NISSA_LOC_VOL_LOOP(ivol)
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
            inv_tmQ2_cgm(cgm_solution,conf,kappa,mass,nmass,niter_max,stopping_residues,source);
	    part_time+=take_time();ninv_tot++;inv_time+=part_time;
	    master_printf("Finished the inversion of S0 theta %d, dirac index %d in %g sec\n",itheta,id,part_time);
	    
	    for(int imass=0;imass<nmass;imass++)
	      { //reconstruct the doublet
		reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass[imass],cgm_solution[imass]);
		master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
		for(int r=0;r<2;r++) //convert the id-th spincolor into the su3spinspin
		  for(int i=0;i<loc_vol;i++)
		    put_spincolor_into_su3spinspin(S0[r][iprop_of(itheta,imass)][i],temp_vec[r][i],id,ic);
	      }
	  }
      }
  
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
      rotate_vol_su3spinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
}

// standard propagators
void calculate_S0_derivative(int ism_lev_so)
{
  smear_additive_su3spinspin(original_source,original_source,ism_lev_so,jacobi_niter_so);
  
  for(int ic=0;ic<3;ic++)
    for(int id=0;id<4;id++)
      { //loop over the source dirac index
	NISSA_LOC_VOL_LOOP(ivol)
	  {
	    get_spincolor_from_su3spinspin(source[ivol],original_source[ivol],id,ic);
	    //put the g5
	    for(int id1=2;id1<4;id1++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic][ri]*=-1;
	  }
	
	//derivative operator applied D_i on the smeared source
	for (int idir=1;idir<4;idir++)
	  {
	    //master_printf("apply the derivative on the source \n");
	    apply_nabla_i_spincolor(der_source, source, conf, idir);
	    
	    for(int itheta=0;itheta<ntheta;itheta++)
	      {
		//adapt the border condition
		put_theta[1]=put_theta[2]=put_theta[3]=theta[itheta];
		adapt_theta(conf,old_theta,put_theta,1,1);
		
		double part_time=-take_time();
		//master_printf("calculate the solution with derivative source \n");
		inv_tmQ2_cgm(cgm_solution,conf,kappa,mass,nmass,niter_max,stopping_residues,der_source);
		part_time+=take_time();ninv_tot++;inv_time+=part_time;
		master_printf("Finished the inversion of S0_der theta %d, dirac index %d in %g sec\n",itheta,id,part_time);
		
		for(int imass=0;imass<nmass;imass++)
		  {
		    reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass[imass],cgm_solution[imass]);
		    master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
		    
		    //derivative_operator applied on the solution
		    for(int jdir=1;jdir<4;jdir++) 
		      { 				
			//reconstruct the doublet										
			apply_nabla_i_spincolor(temp_der_vec[0], temp_vec[0], conf, jdir);
			apply_nabla_i_spincolor(temp_der_vec[1], temp_vec[1], conf, jdir);
			
			for(int r=0;r<2;r++) //convert the id-th spincolor into the su3spinspin
			  for(int i=0;i<loc_vol;i++)
			    put_spincolor_into_su3spinspin(S0_der[r][idir-1][jdir-1][iprop_of(itheta,imass)][i],temp_der_vec[r][i],id,ic);
			
		      }// end of derivative_operator applied on the solution
		  }// loop imass
		
	      }//loop itheta
	    
	  } //end of iteration nabla_i (derivative_source)
      }
  
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
      for(int idir=1;idir<4;idir++)
	for(int jdir=1;jdir<4;jdir++)
	  rotate_vol_su3spinspin_to_physical_basis(S0_der[r][idir-1][jdir-1][ipropS0],!r,!r);	
}



//Calculate and print to file the 2pts
void calculate_all_2pts(int ism_lev_so,int ism_lev_si)
{
  for(int r=0;r<2;r++)
    for(int iprop=0;iprop<npropS0;iprop++)
      smear_additive_su3spinspin(S0[r][iprop],S0[r][iprop],ism_lev_si,jacobi_niter_si);
  
  for(int r=0;r<2;r++)
    for(int iprop=0;iprop<npropS0;iprop++)
      for(int idir=1;idir<4;idir++)
	for(int jdir=1;jdir<4;jdir++)
	  smear_additive_su3spinspin(S0_der[r][idir-1][jdir-1][iprop],S0_der[r][idir-1][jdir-1][iprop],ism_lev_si,jacobi_niter_si);
  
  char path[1024];
  sprintf(path,"%s_%02d_%02d",outfile_2pts,jacobi_niter_so[ism_lev_so],jacobi_niter_si[ism_lev_si]);
  FILE *fout=open_text_file_for_output(path);
  
  contr_time-=take_time();
  
  for(int ispec=0;ispec<nspec;ispec++)
    {
      int ith1=ith_spec[ispec];		
      for(int idir=1;idir<4;idir++)			
	for(int jdir=1;jdir<4;jdir++)
	  for(int ith2=0;ith2<ntheta;ith2++)
	    for(int im2=0;im2<nmass;im2++)
	      for(int r2=0;r2<2;r2++)
		for(int im1=0;im1<nmass;im1++)
		  for(int r1=0;r1<2;r1++)
		    {
		      int ip1=iprop_of(ith1,im1),ip2=iprop_of(ith2,im2);
		      
		      if(rank==0)
			{
			  fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d idir=%d jdir=%d \n",
				  mass[im1],theta[ith1],r1,mass[im2],theta[ith2],r2,idir,jdir);
			  fprintf(fout," smear_source=%d smear_sink=%d\n",jacobi_niter_so[ism_lev_so],jacobi_niter_si[ism_lev_si]);
			}
		      
		      meson_two_points(contr_2pts,op1_2pts,S0_der[r1][idir-1][jdir-1][ip1],op2_2pts,S0[r2][ip2],ncontr_2pts);
		      ncontr_tot+=ncontr_2pts;
		      print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,source_coord[0],"",1);
		      
		      if(rank==0) fprintf(fout,"\n");
		    }
    }
  contr_time+=take_time();
  
  if(rank==0) fclose(fout);
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();
  
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  tot_time-=take_time();
  initialize_semileptonic(arg[1]);
  
  generate_delta_source(original_source,source_coord);
  
  for(int sm_lev_so=0;sm_lev_so<nsm_lev_so;sm_lev_so++)
    {
      calculate_S0(sm_lev_so);
      calculate_S0_derivative(sm_lev_so);
      for(int sm_lev_si=0;sm_lev_si<nsm_lev_si;sm_lev_si++) calculate_all_2pts(sm_lev_so,sm_lev_si);
    }
  
  tot_time+=take_time();
  close_semileptonic();
  
  return 0;
}

