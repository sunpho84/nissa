#include "nissa.h"
#include "addendum.c"

//gauge info
int ngauge_conf;
char conf_path[1024];
quad_su3 *conf,*sme_conf;
double kappa;

int nmass;
double *mass;
double put_theta[4],old_theta[4]={0,0,0,0};

//source data
int source_coord[4]={0,0,0,0};
spincolor *source;
su3spinspin *original_source;

//smearing parameters
double jacobi_kappa,ape_alpha;
int *jacobi_niter,ape_niter;
int nsm_lev;

//vectors for the spinor data
int npropS0;
su3spinspin **S0[2];
spincolor **cgmms_solution,*temp_vec[2];

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
int wall_time;
int ninv_tot=0,ncontr_tot=0;
double tot_time=0,inv_time=0,contr_time=0;

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

//Generate the source for the dirac index
void generate_source()
{ //reset
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
  //Walltime
  read_str_int("Walltime",&wall_time);
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the masses
  
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
  
  read_str_int("NGaugeConf",&ngauge_conf);  
  
  ////////////////////////////////////// end of input reading/////////////////////////////////
  
  //allocate gauge conf and all the needed spincolor and su3spinspin
  conf=nissa_malloc("or_conf",loc_vol+loc_bord,quad_su3);
  sme_conf=nissa_malloc("sm_conf",loc_vol+loc_bord,quad_su3);
  
  //Allocate all the S0 su3spinspin vectors
  npropS0=nmass;
  S0[0]=nissa_malloc("S0[0]",npropS0,su3spinspin*);
  S0[1]=nissa_malloc("S0[1]",npropS0,su3spinspin*);
  for(int iprop=0;iprop<npropS0;iprop++)
    {
      S0[0][iprop]=nissa_malloc("S0[0][iprop]",loc_vol,su3spinspin);
      S0[1][iprop]=nissa_malloc("S0[1][iprop]",loc_vol,su3spinspin);
    }
  
  //Allocate nmass spincolors, for the cgmms solutions
  cgmms_solution=nissa_malloc("cgmms_solution",nmass,spincolor*);
  for(int imass=0;imass<nmass;imass++) cgmms_solution[imass]=nissa_malloc("cgmms_solution[imass]",loc_vol+loc_bord,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol+loc_bord,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol+loc_bord,spincolor);
  
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+loc_bord,spincolor);
  original_source=nissa_malloc("orig_source",loc_vol,su3spinspin);
}

//load the conf, smear it and put boundary cond
void load_gauge_conf()
{
  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_ildg_gauge_conf(conf,conf_path);
  //prepare the smerded version and calculate plaquette
  ape_smear_conf(sme_conf,conf,ape_alpha,ape_niter);
  communicate_lx_gauge_borders(conf);
  communicate_lx_gauge_borders(sme_conf);
  
  double gplaq=global_plaquette_lx_conf(conf);
  master_printf("plaq: %.18g\n",gplaq);
  gplaq=global_plaquette_lx_conf(sme_conf);
  master_printf("smerded plaq: %.18g\n",gplaq);

  //Put the anti-periodic condition on the temporal border
  old_theta[0]=0;
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,0);
}

//Finalization
void close_semileptonic()
{
  nissa_free(original_source);nissa_free(source);
  nissa_free(temp_vec[1]);nissa_free(temp_vec[0]);
  for(int imass=0;imass<nmass;imass++) nissa_free(cgmms_solution[imass]);
  nissa_free(cgmms_solution);
  for(int r=0;r<2;r++)
    {
      for(int iprop=0;iprop<npropS0;iprop++)
	nissa_free(S0[r][iprop]);
      nissa_free(S0[r]);
    }
  nissa_free(sme_conf);
  nissa_free(conf);
  nissa_free(contr_2pts);
  nissa_free(op1_2pts);
  nissa_free(op2_2pts);

  if(rank==0)
    {
      printf("\n");
      printf("Total time: %g, of which:\n",tot_time);
      printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
      printf(" - %02.2f%s to perform %d contr. (%2.2gs avg)\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
    }
  
  close_nissa();
}

//calculate the standard propagators
void calculate_S0(int sm_lev_sour)
{
  for(int id=0;id<4;id++)
    for(int ic=0;ic<3;ic++)
      { //loop over the source dirac index
	for(int ivol=0;ivol<loc_vol;ivol++)
	  {
	    get_spincolor_from_su3spinspin(source[ivol],original_source[ivol],id,ic);
	    //put the g5
	    for(int id1=2;id1<4;id1++) for(int ic1=0;ic1<3;ic1++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic1][ri]*=-1;
	  }
	
	//smerd the source if needed
	if(sm_lev_sour!=0) jacobi_smearing(source,source,sme_conf,jacobi_kappa,jacobi_niter[sm_lev_sour]);
	
	double part_time=-take_time();
	communicate_lx_spincolor_borders(source);
	inv_tmQ2_cgmms(cgmms_solution,source,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	part_time+=take_time();ninv_tot++;inv_time+=part_time;
	master_printf("Finished the inversion of S0, dirac index %d, color %d in %g sec\n",id,ic,part_time);
	
	for(int imass=0;imass<nmass;imass++)
	  { //reconstruct the doublet
	    reconstruct_tm_doublet(temp_vec[0],temp_vec[1],cgmms_solution[imass],conf,kappa,mass[imass]);
	    master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	    for(int r=0;r<2;r++) //convert the id-th spincolor into the su3spinspin
	      for(int i=0;i<loc_vol;i++)
		put_spincolor_into_su3spinspin(S0[r][imass][i],temp_vec[r][i],id,ic);
	  }
      }
  
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
      rotate_vol_su3spinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
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
	      for(int ic=0;ic<3;ic++)
		{
		  for(int ivol=0;ivol<loc_vol;ivol++) get_spincolor_from_su3spinspin(temp_vec[0][ivol],S0[r][imass][ivol],id,ic);
		  jacobi_smearing(temp_vec[1],temp_vec[0],sme_conf,jacobi_kappa,jacobi_niter[sm_lev_sink]-jacobi_niter[sm_lev_sink-1]);
		  for(int ivol=0;ivol<loc_vol;ivol++) put_spincolor_into_su3spinspin(S0[r][imass][ivol],temp_vec[1][ivol],id,ic);
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
		print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,source_coord[0],"",1);
		
		if(rank==0) fprintf(fout,"\n");
	      }
      contr_time+=take_time();
      
      if(rank==0) fclose(fout);
    }
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();
  
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  tot_time-=take_time();
  initialize_semileptonic(arg[1]);
  
  int iconf=0;
  do
    {
      //Gauge path
      int ex=0;
      do
	{
	  read_str(conf_path,1024);
	  read_int(&source_coord[0]);
	  read_str(outfile_2pts,1024);
	  
	  char check_path[1024];
	  sprintf(check_path,"%s_%d%d",outfile_2pts,0,0);
	  ex=file_exists(check_path);
	  if(ex) master_printf("%s already analized, skiping.\n",conf_path);
	  else
	    {
	      master_printf("%s not already analized.\n",conf_path);
	      FILE *fout=open_text_file_for_output(check_path);
	      if(rank==0) fclose(fout);
	    }
	}
      while(ex);
      load_gauge_conf();
      generate_source();
      
      for(int sm_lev_sour=0;sm_lev_sour<nsm_lev;sm_lev_sour++)
	{
	  calculate_S0(sm_lev_sour);
	  calculate_all_2pts(sm_lev_sour);
	}
      iconf++;
    }
  while(iconf<ngauge_conf && take_time()+tot_time<wall_time);
  
  tot_time+=take_time();
  close_semileptonic();
  
  return 0;
}
