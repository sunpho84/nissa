#include "nissa.h"

typedef int int2[2];
typedef int2 interv[2];

//gauge info
int ngauge_conf,nanalized_conf;
char conf_path[1024];
quad_su3 *conf,*unfix_conf;
double kappa;

//mass list
int nmass;
double *mass;
double put_theta[4],old_theta[4]={0,0,0,0};

//output parameters
int write_fixed_conf;
int work_in_physical_base;
int X_space_prop_prec,P_space_prop_prec;
int n_X_interv,n_P_interv;
interv *X_interv,*P_interv;

//source data
int source_coord[4]={0,0,0,0};
spincolor *source;
su3spinspin *original_source;

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
char outfolder[1024];

//timings
int ninv_tot=0,ncontr_tot=0,nsaved_prop_tot=0;
int wall_time;
double tot_time=0,inv_time=0,contr_time=0,fix_time=0;
double fft_time=0,save_prop_time=0,load_time=0,filter_prop_time=0;

//Read the information regarding the format in which to save prop
int read_prop_prec(const char *name)
{
  int prec;
  
  read_str_int(name,&prec);
  if(prec!=0 && prec!=32 && prec!=64)
    crash("Error, asking to save %s in %d precision (only 0, 32, 64 available)",name,prec);
  
  return prec;
}

//Read the list of subset of X or P prop to save
interv* read_subset_list(int *n_subset,const char *name,const char *tag)
{
  interv *inte;
  
  read_str_int(name,n_subset);
  inte=nissa_malloc(tag,*n_subset,interv);
  for(int isub=0;isub<(*n_subset);isub++)
    for(int mu=0;mu<2;mu++)
      for(int iext=0;iext<2;iext++)
	read_int(inte[isub][mu]+iext);
  
  return inte;
}

//write all ptopagators
void write_all_propagators(const char *name,int prec)
{
  for(int r=0;r<2;r++)
    for(int imass=0;imass<nmass;imass++)
      {
	char path[1024];
	sprintf(path,"%s/%sprop/r%1d_im%02d",outfolder,name,r,imass);
	save_prop_time-=take_time();
	write_su3spinspin(path,S0[r][imass],prec);
	save_prop_time+=take_time();
	nsaved_prop_tot++;
      }
}

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
      lx[idir]=source_coord[idir]-proc_coord[idir]*loc_size[idir];
      islocal&=(lx[idir]>=0);
      islocal&=(lx[idir]<loc_size[idir]);
    }
  
  if(islocal)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	original_source[loclx_of_coord(lx)][ic][ic][id][id][0]=1;
}

//Parse all the input file
void initialize_Zcomputation(char *input_path)
{
  open_input(input_path);

  // 1) Read information about the gauge conf
  
  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Wall_time
  read_str_int("WallTime",&wall_time);
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the masses
  
  //Read the masses
  read_list_of_doubles("NMass",&nmass,&mass);

  // 3) Info about the inverter

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
  
  // 4) contraction list for two points
  
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
  
  // 5) Information on output
  
  //save fixed conf?
  read_str_int("WriteFixedConf",&write_fixed_conf);
  
  //work in physical space?
  read_str_int("WorkInPhysicalBase",&work_in_physical_base);
  
  //precision in which to save the X and P prop (0, 32, 64)
  X_space_prop_prec=read_prop_prec("XSpacePropPrec");
  P_space_prop_prec=read_prop_prec("PSpacePropPrec");
    
  //read subsets of X and P space props to save
  X_interv=read_subset_list(&n_X_interv,"NXSpacePropInterv","XInterv");
  P_interv=read_subset_list(&n_P_interv,"NPSpacePropInterv","PInterv");
    
  read_str_int("NGaugeConf",&ngauge_conf);  
  
  ////////////////////////////////////// end of input reading/////////////////////////////////
  
  //allocate gauge conf and all the needed spincolor and su3spinspin
  conf=nissa_malloc("or_conf",loc_vol+loc_bord,quad_su3);
  unfix_conf=nissa_malloc("unfix_conf",loc_vol+loc_bord,quad_su3);
  
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
  reco_solution[0]=nissa_malloc("reco_sol[0]",loc_vol+loc_bord,spincolor);
  reco_solution[1]=nissa_malloc("reco_sol[1]",loc_vol+loc_bord,spincolor);
  
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+loc_bord,spincolor);
  original_source=nissa_malloc("orig_source",loc_vol,su3spinspin);
}

//load the conf, fix it and put boundary cond
void load_gauge_conf()
{
  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  load_time-=take_time();
  read_gauge_conf(unfix_conf,conf_path);
  load_time+=take_time();
  //prepare the fixed version and calculate plaquette
  double elaps_time=-take_time();
  landau_gauge_fix(conf,unfix_conf,1.e-20);
  elaps_time+=take_time();
  fix_time+=elaps_time;
  master_printf("Fixed conf in %lg sec\n",elaps_time);
  communicate_gauge_borders(conf);
  communicate_gauge_borders(unfix_conf);
  
  if(write_fixed_conf)
    {
      char temp[1024];
      sprintf(temp,"%s/fixed_conf",outfolder);
      write_gauge_conf(temp,conf);
    }    
  
  master_printf("plaq: %.18g\n",global_plaquette(conf));
  master_printf("unfix plaq: %.18g\n",global_plaquette(unfix_conf));

  //Put the anti-periodic condition on the temporal border
  old_theta[0]=0;
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,0);
}

//Finalization
void close_Zcomputation()
{
  nissa_free(original_source);nissa_free(source);
  nissa_free(reco_solution[1]);nissa_free(reco_solution[0]);
  for(int imass=0;imass<nmass;imass++) nissa_free(cgmms_solution[imass]);
  nissa_free(cgmms_solution);
  for(int r=0;r<2;r++)
    {
      for(int iprop=0;iprop<npropS0;iprop++)
	nissa_free(S0[r][iprop]);
      nissa_free(S0[r]);
    }
  nissa_free(unfix_conf);
  nissa_free(conf);
  nissa_free(contr_2pts);
  nissa_free(op1_2pts);
  nissa_free(op2_2pts);
  nissa_free(X_interv);
  nissa_free(P_interv);
  
  master_printf("\n");
  master_printf("Total time: %lg sec (%lg average per conf), of which:\n",tot_time,tot_time/nanalized_conf);
  master_printf(" - %02.2f%s to load %d conf. (%2.2gs avg)\n",load_time/tot_time*100,"%",nanalized_conf,load_time/nanalized_conf);
  master_printf(" - %02.2f%s to fix %d conf. (%2.2gs avg)\n",fix_time/tot_time*100,"%",nanalized_conf,fix_time/nanalized_conf);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d contr. (%2.2gs avg)\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
  master_printf(" - %02.2f%s to save %d su3spinspins. (%2.2gs avg)\n",save_prop_time/tot_time*100,"%",nsaved_prop_tot,save_prop_time/nsaved_prop_tot);
  master_printf(" - %02.2f%s to filter propagators. (%2.2gs avg per conf)\n",filter_prop_time/tot_time*100,"%",filter_prop_time/nanalized_conf);
  
  close_nissa();
}

//calculate the standard propagators
void calculate_S0()
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
	
	double part_time=-take_time();
	communicate_lx_spincolor_borders(source);
	inv_Q2_cgmms(cgmms_solution,source,NULL,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	part_time+=take_time();ninv_tot++;inv_time+=part_time;
	master_printf("Finished the inversion of S0, dirac index %d, color %d in %g sec\n",id,ic,part_time);
	
	//reconstruct the doublet
	for(int imass=0;imass<nmass;imass++)
	  {
	    reconstruct_doublet(reco_solution[0],reco_solution[1],cgmms_solution[imass],conf,kappa,mass[imass]);
	    
	    for(int r=0;r<2;r++) //convert the id-th spincolor into the su3spinspin
	      for(int i=0;i<loc_vol;i++)
		put_spincolor_into_su3spinspin(S0[r][imass][i],reco_solution[r][i],id,ic);
	    
	    master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	  }
      }
  
  //rotate only if working in the physical base asked
  if(work_in_physical_base)
    for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
      for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
      	rotate_vol_su3spinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
	
  //save propagators if asked
  if(X_space_prop_prec) write_all_propagators("X",X_space_prop_prec);
}

//compute the fft of all propagators
void compute_fft(double sign)
{
  fft_time-=take_time();
  for(int r=0;r<2;r++)
    for(int imass=0;imass<nmass;imass++)
      {
	fft4d((complex*)S0[r][imass],(complex*)S0[r][imass],144,sign,0);
	//multiply by the conjugate of the fft of the source
	for(int imom=0;imom<loc_vol;imom++)
	  {
	    double arg=0;
	    for(int mu=0;mu<4;mu++) arg+=((double)glb_coord_of_loclx[imom][mu]*source_coord[mu])/glb_size[mu];
	    arg*=-sign*2*M_PI;
	    complex f={cos(arg),sin(arg)};

	    for(int ic_si=0;ic_si<3;ic_si++)
	      for(int ic_so=0;ic_so<3;ic_so++)
		for(int id_si=0;id_si<4;id_si++)
		  for(int id_so=0;id_so<4;id_so++)
		    safe_complex_prod(S0[r][imass][imom][ic_si][ic_so][id_si][id_so],
				      S0[r][imass][imom][ic_si][ic_so][id_si][id_so],
				      f);
	  }
      }
  fft_time+=take_time();

  //save propagators if asked
  if(P_space_prop_prec) write_all_propagators("P",P_space_prop_prec);
}

//filter the propagators
void print_propagator_subset(const char *name,int nsubset,interv *inte)
{
  filter_prop_time-=take_time();
  
  colorspincolorspin *buf=nissa_malloc("buf",2*nmass,colorspincolorspin);
  
  //build oputput file name
  char outfile_fft[1024];
  sprintf(outfile_fft,"%s/%s",outfolder,name);

  //open oputput file for concurent access from different ranks
  MPI_File fout;
  int rc=MPI_File_open(cart_comm,outfile_fft,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fout);
  if(rc) decript_MPI_error(rc,"Unable to open file: %s",outfile_fft);
  
  //loop over momenta subsets
  int ip=0;
  for(int isub=0;isub<nsubset;isub++)
    {
      //loop over momenta in each se
      int glb_ip[4];
      for(glb_ip[0]=inte[isub][0][0];glb_ip[0]<=inte[isub][0][1];glb_ip[0]++)
	for(glb_ip[1]=inte[isub][1][0];glb_ip[1]<=inte[isub][1][1];glb_ip[1]++)
	  for(glb_ip[2]=inte[isub][1][0];glb_ip[2]<=inte[isub][1][1];glb_ip[2]++)
	    for(glb_ip[3]=inte[isub][1][0];glb_ip[3]<=inte[isub][1][1];glb_ip[3]++)
	      {
		//identify the rank hosting this element
		int hosting=rank_hosting_site_of_coord(glb_ip);
		
		//if the hosting is the current rank write the data at the correct location
		if(hosting==cart_rank)
		  {
		    //find local index
		    int loc_ip[4];
		    for(int mu=0;mu<4;mu++) loc_ip[mu]=glb_ip[mu]%loc_size[mu];
		    int ilp=loclx_of_coord(loc_ip);
		    
		    //bufferize data
		    for(int r=0;r<2;r++)
		      for(int imass=0;imass<nmass;imass++)
			for(int ic_so=0;ic_so<3;ic_so++)
			  for(int id_so=0;id_so<4;id_so++)
			    for(int ic_si=0;ic_si<3;ic_si++)
			      for(int id_si=0;id_si<4;id_si++)
				memcpy(buf[r*nmass+imass][ic_so][id_so][ic_si][id_si],S0[r][imass][ilp][ic_si][ic_so][id_si][id_so],
				       sizeof(complex));
		    
		    //find the position in the file where to write data
		    int offset=ip*2*nmass*sizeof(colorspincolorspin);
		    
		    MPI_File_write_at(fout,offset,buf,nmass*2*16,MPI_SU3,MPI_STATUS_IGNORE);
		  }
		ip++;
	      }
    }
  
  nissa_free(buf);
  
  MPI_File_close(&fout);

  filter_prop_time+=take_time();
}

//Calculate and print to file the 2pts
void calculate_all_2pts()
{
  char outfile_2pts[1024];
  sprintf(outfile_2pts,"%s/2pts",outfolder);
  FILE *fout=open_text_file_for_output(outfile_2pts);
  
  //perform the contractions
  contr_time-=take_time();      
  for(int im2=0;im2<nmass;im2++)
    for(int r2=0;r2<2;r2++)
      for(int im1=0;im1<nmass;im1++)
	for(int r1=0;r1<2;r1++)
	  {
	    master_fprintf(fout," # m1=%f r1=%d , m2=%f r2=%d \n",mass[im1],r1,mass[im2],r2);
	    
	    meson_two_points(contr_2pts,op1_2pts,S0[r1][im1],op2_2pts,S0[r2][im2],ncontr_2pts);
	    ncontr_tot+=ncontr_2pts;
	    print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,source_coord[0],"",1);
	    
	    if(rank==0) fprintf(fout,"\n");
	  }
  contr_time+=take_time();
}

//read the conf parameters
int read_conf_parameters(int *iconf)
{
  int ok_conf;
  
  do
    {
      //Gauge path
      read_str(conf_path,1024);
      
      //Source position
      read_int(&(source_coord[0]));
      read_int(&(source_coord[1]));
      read_int(&(source_coord[2]));
      read_int(&(source_coord[3]));
      
      //Folder
      read_str(outfolder,1024);
      master_printf("Considering configuration %s\n",conf_path);
      ok_conf=!(dir_exists(outfolder));
      if(ok_conf)
	{
	  if(rank==0)
	    {
	      char temp[1024];
	      mkdir(outfolder,S_IRWXU);
	      sprintf(temp,"%s/Pprop",outfolder);
	      mkdir(temp,S_IRWXU);
	      sprintf(temp,"%s/Xprop",outfolder);
	      mkdir(temp,S_IRWXU);
	    }
	  master_printf("Configuration not already analized, starting.\n");
	}
      else
	master_printf("Configuration already analized, skipping.\n");
      (*iconf)++;
    }
  while(!ok_conf && (*iconf)<ngauge_conf);
  
  return ok_conf;
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
  if(enough_time) master_printf("Continuing!\n");
  else master_printf("Not enough time, exiting!\n");
  
  return enough_time;
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();
  
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  tot_time-=take_time();
  initialize_Zcomputation(arg[1]);
  
  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf && enough_time && read_conf_parameters(&iconf))
    {    
      load_gauge_conf();
      generate_source();
      
      //X space
      calculate_S0();
      calculate_all_2pts();
      if(n_X_interv) print_propagator_subset("Xprop/subset",n_X_interv,X_interv);
      
      //P space
      compute_fft(-1);
      if(n_P_interv) print_propagator_subset("Pprop/subset",n_P_interv,P_interv);      
      
      nanalized_conf++;
      
      enough_time=check_remaining_time();
    }
  
  if(iconf==ngauge_conf) master_printf("Finished all the conf!\n");

  tot_time+=take_time();
  close_Zcomputation();
  
  return 0;
}
