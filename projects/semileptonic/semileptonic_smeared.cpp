#include <stdlib.h>
#include <string.h>
#include <math.h>

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
  
  If specified, the program computes also the 2 points correlations funcions with derivative
  operator on the source. It is possible to spcify to compute derivative-source propagator 
  only for a subset of masses: the restriction is imposed to Q1.
  
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

//Wilson clover or Tm?
int Wclov_tm;
double cSW;

//gauge info
int ngauge_conf,nanalyzed_conf=0;
char conf_path[1024],outfolder[1024];
quad_su3 *conf,*sme_conf;
as2t_su3 *Pmunu;
double kappa;
double put_theta[4],old_theta[4]={0,0,0,0};

//list of masses and theta
int which_r_S0;
int use_cgm_S0,use_cgm_S1;
int nmassS0,nthetaS0;
int nmassS1,nthetaS1;
double *massS0,*thetaS0;
double *massS1,*thetaS1;
int start_massS0der;
int nmassS0der;
double *massS0der;

//source data
int save_source;
int seed,noise_type;
coords source_coord;
spincolor *source;
prop_type *original_source;

//smearing parameters
double gaussian_kappa,ape_alpha;
int ape_niter;
int *gaussian_niter_so,nsm_lev_so;
int *gaussian_niter_si,nsm_lev_si;
int *gaussian_niter_se,nsm_lev_se;

//vectors for the S0 props
int compute_der,nmuS;
int npropS0;
int load_S0,save_S0;
prop_type **S0[2];
int ncgm_solution;
spincolor **cgm_solution,*temp_vec[2];
int ithetaS0_min,ithetaS0_max;

//cgm inverter parameters
double *stopping_residues_S0;
double *stopping_residues_S1;
int niter_max=100000;

//two points contractions
int ncontr_2pts;
int only_standing_2pts;
int only_charged_2pts;
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
int contr_3pts_up_to_S0_mass;
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
int ipropS0(int itheta,int imass,int mu_der)
{
  if(mu_der==0) return itheta*nmassS0+imass;
  else return nthetaS0*nmassS0+((mu_der-1)*nthetaS0+itheta)*nmassS0der+imass;
}
int ipropS1(int itheta,int imass)
{return itheta*nmassS1+imass;}

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
  generate_spindiluted_source(original_source,nissa_rnd_type_map[noise_type],source_coord[0]);
  
  //if asked, save the source
  if(save_source)
    {
      char outpath[1024];
      sprintf(outpath,"%s/source",outfolder);
      write_colorspinspin(outpath,original_source,64);
    }
#endif
}

//Generate a sequential source for S1
void generate_sequential_source(int ispec)
{
  int r=r_spec[ispec];
  
  select_propagator_timeslice(sequential_source,sequential_source,0);
  
  master_printf("\nCreating the sequential source for spectator %d\n",ispec);
  nissa_loc_vol_loop(ivol)
    {
      //put to zero everywhere but on the slice
      if(glb_coord_of_loclx[ivol][0]!=(source_coord[0]+tsep)%glb_size[0])
	memset(sequential_source[ivol],0,sizeof(prop_type));
      else
	{
	  //avoid to put g5, beacuse commutate with (i+-g5)/sqrt(2) and cancel with those of the QQ
	  memcpy(sequential_source[ivol],S0[r][ipropS0(ith_spec[ispec],imass_spec[ispec],0)][ivol],sizeof(prop_type));
	  if(Wclov_tm==1) //if doing tm
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

//Apply the covariant derivative having direction mu on a spincolor
void apply_nabla_i(spincolor *out,spincolor *in,quad_su3 *conf,int mu)
{
  spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
  
  memset(temp,0,loc_vol*sizeof(spincolor));
  communicate_lx_spincolor_borders(in);
  
  nissa_loc_vol_loop(ix)
    {
      int Xup,Xdw;
      Xup=loclx_neighup[ix][mu];
      Xdw=loclx_neighdw[ix][mu];
      
      unsafe_su3_prod_spincolor(             temp[ix],conf[ix][mu] ,in[Xup]);
      unsafe_su3_dag_subt_the_prod_spincolor(temp[ix],conf[Xdw][mu],in[Xdw]);
    }
  
  memcpy(out,temp,sizeof(spincolor)*loc_vol);
  nissa_free(temp);
  
  set_borders_invalid(out);
}

void apply_nabla_i(prop_type *out,prop_type *in,quad_su3 *conf,int mu)
{
  spincolor *temp_in=nissa_malloc("temp_in",loc_vol+bord_vol,spincolor);
  spincolor *temp_out=nissa_malloc("temp_out",loc_vol,spincolor);
  
#ifdef POINT_SOURCE_VERSION
  for(int ic=0;ic<3;ic++)
#endif
    for(int id=0;id<4;id++)
      {
#ifdef POINT_SOURCE_VERSION
	nissa_loc_vol_loop(ivol) get_spincolor_from_su3spinspin(temp_in[ivol],in[ivol],id,ic);
#else
	nissa_loc_vol_loop(ivol) get_spincolor_from_colorspinspin(temp_in[ivol],in[ivol],id);
#endif
	set_borders_invalid(temp_in);
	apply_nabla_i(temp_out,temp_in,conf,mu);
#ifdef POINT_SOURCE_VERSION
	nissa_loc_vol_loop(ivol) put_spincolor_into_su3spinspin(out[ivol],temp_out[ivol],id,ic);
#else
	nissa_loc_vol_loop(ivol) put_spincolor_into_colorspinspin(out[ivol],temp_out[ivol],id);
#endif
	set_borders_invalid(out);
      }
  
  nissa_free(temp_in);
  nissa_free(temp_out);
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
  //Decide if twisted(1) or clover run
  read_str_int("TwistedMassRun",&Wclov_tm);
  //Kappa is read really only for tm
  if(Wclov_tm==1) read_str_double("Kappa",&kappa);
  read_str_double("cSW",&cSW);
  
  // 2) Read information about the source
  
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  //Read the noise type
  read_str_int("NoiseType",&noise_type);
  //read whether we want to save the source
  read_str_int("SaveSource",&save_source);

  // 3) Smearing parameters
  
  //Smearing parameters
  read_str_double("ApeAlpha",&ape_alpha);
  read_str_int("ApeNiter",&ape_niter);
  read_str_double("GaussianKappa",&gaussian_kappa);
  read_list_of_ints("GaussianNiterSo",&nsm_lev_so,&gaussian_niter_so);
  for(int iter=1;iter<nsm_lev_so;iter++)
    if(gaussian_niter_so[iter]<gaussian_niter_so[iter-1])
      crash("Error, gaussian lev sou %d minor than %d (%d, %d)!\n",iter,iter-1,gaussian_niter_so[iter],gaussian_niter_so[iter-1]);
  read_list_of_ints("GaussianNiterSe",&nsm_lev_se,&gaussian_niter_se);
  for(int iter=1;iter<nsm_lev_se;iter++)
    if(gaussian_niter_se[iter]<gaussian_niter_se[iter-1])
      crash("Error, gaussian lev seq %d minor than %d (%d, %d)!\n",iter,iter-1,gaussian_niter_se[iter],gaussian_niter_se[iter-1]);
  read_list_of_ints("GaussianNiterSi",&nsm_lev_si,&gaussian_niter_si);
  for(int iter=1;iter<nsm_lev_si;iter++)
    if(gaussian_niter_si[iter]<gaussian_niter_si[iter-1])
      crash("Error, gaussian lev seq %d minor than %d (%d, %d)!\n",iter,iter-1,gaussian_niter_si[iter],gaussian_niter_si[iter-1]);
  
  // 4) Read list of masses and of thetas

  //the number or fr can be chosen only if tm, otherwise use 0
  if(Wclov_tm==1) read_str_int("WhichRS0",&which_r_S0);
  else which_r_S0=0;
  
  //the cgm can be (for the time being) used only for both r, therefore also only for tm
  if(Wclov_tm==0) use_cgm_S0=0;
  else
    if(which_r_S0!=2) read_str_int("UseCgmS0",&use_cgm_S0);
    else use_cgm_S0=1;
  
  read_list_of_double_pairs("MassResiduesS0",&nmassS0,&massS0,&stopping_residues_S0);
  read_list_of_doubles("NThetaS0",&nthetaS0,&thetaS0);
  read_str_int("SaveS0",&save_S0);
  if(save_S0==0) read_str_int("LoadS0",&load_S0);
  else load_S0=0;
  
  // 5) contraction list for two points
  
  read_str_int("ComputeDerivativeCorrelations",&compute_der);
  if(compute_der)
    {
      nmuS=4;
      read_str_int("StartDerivativeInversionFromIMass",&start_massS0der);
      nmassS0der=nmassS0-start_massS0der;
      massS0der=massS0+start_massS0der;
    }
  else nmuS=1;
  
  read_str_int("OnlyStandingTwoPoints",&only_standing_2pts);
  read_str_int("OnlyChargedTwoPoints",&only_charged_2pts);
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
  
  // 6) Read list of masses and of thetas for S1
  
  if(Wclov_tm==1) read_str_int("UseCgmS1",&use_cgm_S1);
  read_list_of_double_pairs("MassResiduesS1",&nmassS1,&massS1,&stopping_residues_S1);
  read_list_of_doubles("NThetaS1",&nthetaS1,&thetaS1);
  read_str_int("ContrThreePointsUpToS0Mass",&contr_3pts_up_to_S0_mass);
  
  // 7) three points functions
  
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
      
      if(ith_spec[ispec]<0||ith_spec[ispec]>=nthetaS0)    crash("theta for ispec %d out of bounds",ispec);
      if(imass_spec[ispec]<0||imass_spec[ispec]>=nmassS0) crash("mass for ispec %d out of bounds",ispec);
      if(r_spec[ispec]<0||r_spec[ispec]>=2)               crash("r for ispec %d out of bounds",ispec);
      if(which_r_S0!=2&&r_spec[ispec]!=which_r_S0)        crash("r for ispec %d uncomputed",ispec);
      
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
  conf=nissa_malloc("or_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  sme_conf=nissa_malloc("sm_conf",loc_vol+bord_vol,quad_su3);
  Pmunu=nissa_malloc("Pmunu",loc_vol,as2t_su3);
  
  //Allocate all the S0 prop_type vectors
  npropS0=nthetaS0*nmassS0;
  if(compute_der) npropS0+=3*nthetaS0*nmassS0der;
  S0[0]=nissa_malloc("S0[0]",npropS0,prop_type*);
  S0[1]=nissa_malloc("S0[1]",npropS0,prop_type*);
  for(int iprop=0;iprop<npropS0;iprop++)
    for(int r=0;r<2;r++)
      if(which_r_S0==2||which_r_S0==r) S0[r][iprop]=nissa_malloc("S0[r]",loc_vol,prop_type);
  
  //Allocate nmass spincolors, for the cgm solutions
  ncgm_solution=max_int(nmassS0,nmassS1);
  cgm_solution=nissa_malloc("cgm_solution",ncgm_solution,spincolor*);
  for(int imass=0;imass<ncgm_solution;imass++) cgm_solution[imass]=nissa_malloc("cgm_solution",loc_vol+bord_vol,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol,spincolor);
  
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  original_source=nissa_malloc("original_source",loc_vol,prop_type);
  
  //Allocate one prop_type for the chromo-contractions
  ch_prop=nissa_malloc("chromo-prop",loc_vol,prop_type);
  
  //Allocate all the S1 prop_type vectors
  npropS1=nthetaS1*nmassS1;
  S1=nissa_malloc("S1",loc_vol,prop_type*);
  for(int iprop=0;iprop<npropS1;iprop++) S1[iprop]=nissa_malloc("S1[i]",loc_vol,prop_type);
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
      
      //Check if the conf has been finished or is already running
      master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
      char fin_file[1024],run_file[1024];
      sprintf(fin_file,"%s/finished",outfolder);
      sprintf(run_file,"%s/running",outfolder);
      ok_conf=!(file_exists(fin_file)) && !(file_exists(run_file));
      
      //if not finished
      if(ok_conf)
	{
	  master_printf(" Configuration \"%s\" not yet analyzed, starting",conf_path);
	  if(!dir_exists(outfolder))
	    {
	      int ris=create_dir(outfolder);
	      if(ris==0) master_printf(" Output path \"%s\" not present, created.\n",outfolder);
	      else
		crash(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
	    }
	  file_touch(run_file);
	}
      else
	master_printf(" In output path \"%s\" terminating file already present: configuration \"%s\" already analyzed, skipping.\n",outfolder,conf_path);
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
  for(int iprop=0;iprop<npropS0;iprop++)
    for(int r=0;r<2;r++)
      if(which_r_S0==2||which_r_S0==r) nissa_free(S0[r][iprop]);
  for(int iprop=0;iprop<npropS1;iprop++) nissa_free(S1[iprop]);
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
  for(int imass=0;imass<ncgm_solution;imass++) nissa_free(cgm_solution[imass]);
  nissa_free(cgm_solution);
  nissa_free(source);nissa_free(original_source);
  
  close_nissa();
}

//smear addditivily a propagator
void smear_additive_propagator(prop_type *out,prop_type *in,int ism_lev,int *gaussian_niter)
{
  spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  
  int nsme=gaussian_niter[ism_lev];
  if(ism_lev>0) nsme-=gaussian_niter[ism_lev-1];
  
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
	
	gaussian_smearing(temp,temp,sme_conf,gaussian_kappa,nsme);
	
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
  smear_additive_propagator(original_source,original_source,ism_lev_so,gaussian_niter_so);
  master_printf("\n");
  
  //loop over derivative of the source
  for(int muS=0;muS<nmuS;muS++)
    {
      //decide parameters of inverter
      double *mass=(muS==0)?massS0:massS0der;
      int nmass=(muS==0)?nmassS0:nmassS0der;
      double *stopping_residues=(muS==0)?stopping_residues_S0:stopping_residues_S0+start_massS0der;
      
      //loop over the source dirac index
#ifdef POINT_SOURCE_VERSION
      for(int ic=0;ic<3;ic++)
#endif
	for(int id=0;id<4;id++)
	  { 
	    //put the g5
#ifdef POINT_SOURCE_VERSION
	    get_spincolor_from_su3spinspin(source,original_source,id,ic);
#else
	    get_spincolor_from_colorspinspin(source,original_source,id);
#endif
	    //add gamma5 apart if using cg or tm 
	    if(Wclov_tm==0||use_cgm_S0||(Wclov_tm==1&&cSW!=0)) safe_dirac_prod_spincolor(source,base_gamma[5],source);
	    
	    //if needed apply nabla
	    if(muS>0)
	      {
		//remove the border condition
		put_theta[1]=put_theta[2]=put_theta[3]=0;
		adapt_theta(conf,old_theta,put_theta,1,1);
		
		apply_nabla_i(source,source,conf,muS);
	      }
	    
	    for(int itheta=0;itheta<nthetaS0;itheta++)
	      {
		//adapt the border condition
		put_theta[1]=put_theta[2]=put_theta[3]=thetaS0[itheta];
		adapt_theta(conf,old_theta,put_theta,1,1);
		
		//invert
		if(!load_S0)
		  {
		    double part_time=-take_time();
		    
		    //decide if to use multimass or single mass
		    if(use_cgm_S0)
		      {
			if(cSW==0) inv_tmQ2_cgm(cgm_solution,conf,kappa,mass,nmass,niter_max,stopping_residues,source);
			else inv_tmclovQ2_cgm(cgm_solution,conf,kappa,cSW,Pmunu,mass,nmass,niter_max,stopping_residues,source);
		      }
		    else
		      for(int imass=0;imass<nmass;imass++)
			{
			  //the sign of mass depends on r
			  double m=mass[imass];
			  if(Wclov_tm==1)
			    {
			      master_printf("Ma no!\n");
			      if(which_r_S0==0) m*=-1;
			      
			      if(cSW==0) inv_tmD_cg_eoprec_eos(cgm_solution[imass],NULL,conf,kappa,m,niter_max,5,stopping_residues[imass],source);
			      else inv_tmclovQ_cg(cgm_solution[imass],NULL,conf,kappa,cSW,Pmunu,m,niter_max,5,stopping_residues[imass],source);
			    }
			  else //m=kappa
			    inv_WclovQ_cg(cgm_solution[imass],NULL,conf,m,cSW,Pmunu,niter_max,5,stopping_residues[imass],source);
			  
			  master_printf("Finished submass[%d]=%lg\n",imass,m);
			}
		    part_time+=take_time();ninv_tot++;inv_time+=part_time;
		    master_printf("Finished the inversion of S0 theta %d, ",itheta);
		    if(compute_der) master_printf("source derivative %d ",muS);

#ifdef POINT_SOURCE_VERSION
		    master_printf("color index %d ",ic);
#endif	    
		    master_printf("dirac index %d in %g sec\n",id,part_time);
		  }
		
		//read or writ, if needed
		if(save_S0||load_S0)
		  for(int imass=0;imass<nmass;imass++)
		    {
		      int ip=ipropS0(itheta,imass,muS);
		      
		      //path for S0
		      char path[1024];
#ifdef POINT_SOURCE_VERSION
		      sprintf(path,"%s/S0_QD_sosm%02d_iprop%d.id%02d.ic%02d",outfolder,ism_lev_so,ip,id,ic);
#else
		      sprintf(path,"%s/S0_QD_sosm%02d_iprop%d.id%02d",outfolder,ism_lev_so,ip,id);
#endif
		      if(save_S0) write_spincolor(path,cgm_solution[imass],64);
		      else        read_spincolor(cgm_solution[imass],path);
		    }

		//reconstruct the doublet
		for(int imass=0;imass<nmass;imass++)
		  {
		    if(Wclov_tm==1&&use_cgm_S0)
		      {
			reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass[imass],cgm_solution[imass]);
			master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
		      }
		    else memcpy(temp_vec[which_r_S0],cgm_solution[imass],sizeof(spincolor)*loc_vol);
		    for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
		      if(which_r_S0==r||which_r_S0==2)
			{
#ifdef POINT_SOURCE_VERSION
			  put_spincolor_into_su3spinspin(S0[r][ipropS0(itheta,imass,muS)],temp_vec[r],id,ic);
#else
			  put_spincolor_into_colorspinspin(S0[r][ipropS0(itheta,imass,muS)],temp_vec[r],id);
#endif
			}
		  }
	      }
	  }
    }

  //rotate to physical basis if doing tm
  if(Wclov_tm==1)
    for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
      if(which_r_S0==r||which_r_S0==2)
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
  smear_additive_propagator(sequential_source,sequential_source,ism_lev_se,gaussian_niter_se);
  smear_additive_propagator(sequential_source,sequential_source,ism_lev_se,gaussian_niter_se);
  master_printf("\n");
  
  //loop over seq
#ifdef POINT_SOURCE_VERSION
  for(int ic=0;ic<3;ic++)
#endif
    for(int id=0;id<4;id++)
      { 
#ifdef POINT_SOURCE_VERSION
	get_spincolor_from_su3spinspin(source,sequential_source,id,ic);
#else
	get_spincolor_from_colorspinspin(source,sequential_source,id);
#endif
	safe_dirac_prod_spincolor(source,base_gamma[5],source);
	
	//if inverting Q
	if(Wclov_tm==0||use_cgm_S1||(Wclov_tm==1&&cSW!=0)) safe_dirac_prod_spincolor(source,base_gamma[5],source); 
	
	for(int itheta=0;itheta<nthetaS1;itheta++)
	  {
	    //adapt the border condition
	    put_theta[1]=put_theta[2]=put_theta[3]=thetaS1[itheta];
	    adapt_theta(conf,old_theta,put_theta,1,1);
	    
	    double part_time=-take_time();
	    
	    //decide to use one or the other inverters
	    if(use_cgm_S1)
	      {
		if(cSW==0) inv_tmQ2_cgm(cgm_solution,conf,kappa,massS1,nmassS1,niter_max,stopping_residues_S1,source);
		else inv_tmclovQ2_cgm(cgm_solution,conf,kappa,cSW,Pmunu,massS1,nmassS1,niter_max,stopping_residues_S1,source);
	      }
	    else
	      for(int imass=0;imass<nmassS1;imass++)
		{
		  //since r0==0 implicates r1=1, mass=-mass when r0=1
		  double m=massS1[imass];
		  if(Wclov_tm==1)
		    {
		      if(r_spec[ispec]==1) m*=-1;
		      if(cSW==0) inv_tmD_cg_eoprec_eos(cgm_solution[imass],NULL,conf,kappa,m,niter_max,5,stopping_residues_S1[imass],source);
		      else inv_tmclovQ_cg(cgm_solution[imass],NULL,conf,kappa,cSW,Pmunu,m,niter_max,5,stopping_residues_S1[imass],source);
		    }
		  else inv_WclovQ_cg(cgm_solution[imass],NULL,conf,m,cSW,Pmunu,niter_max,5,stopping_residues_S1[imass],source);
		}
	    
	    part_time+=take_time();ninv_tot++;inv_time+=part_time;
	    master_printf("Finished the inversion of S1 theta %d, seq sme lev %d, dirac index %d in %g sec\n",itheta,ism_lev_se,id,part_time);
	    
	    for(int imass=0;imass<nmassS1;imass++)
	      {
		//reconstruct the doublet: r(S1)=!r(spec), so we have to multiply by Q+ if r(spec)==1 and Q- if 0
		double reco_mass=-massS1[imass];
		if(r_spec[ispec]==1) reco_mass=-reco_mass;
		//use temp_vec[0] as temporary storage
		if(Wclov_tm==1&&use_cgm_S1)
		  {
		    apply_tmQ(temp_vec[0],conf,kappa,reco_mass,cgm_solution[imass]);
		    master_printf("Mass %d (%g) reconstructed \n",imass,massS1[imass]);
		  }
		else memcpy(temp_vec[0],cgm_solution[imass],sizeof(spincolor)*loc_vol);

#ifdef POINT_SOURCE_VERSION
		put_spincolor_into_su3spinspin(S1[ipropS1(itheta,imass)],temp_vec[0],id,ic);
#else
		put_spincolor_into_colorspinspin(S1[ipropS1(itheta,imass)],temp_vec[0],id);
#endif
	      }
	  }
      }
  
  //put the (1+-ig5)/sqrt(2) factor if tm. On the source rotate as r_spec, on the sink as !r_spec
  if(Wclov_tm==1)
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
  prop_type *temp_der=(compute_der==1)?nissa_malloc("temp_der",loc_vol,prop_type):NULL;
  
  smear_time-=take_time();
    
  for(int r=0;r<2;r++)
    if(which_r_S0==2||which_r_S0==r)
      for(int iprop=0;iprop<npropS0;iprop++)
	smear_additive_propagator(S0[r][iprop],S0[r][iprop],ism_lev_si,gaussian_niter_si);
  
  double temp_time=take_time();
  smear_time+=temp_time;
  contr_2pts_time-=temp_time;
  
  char path[1024];
  sprintf(path,"%s/2pts_%02d_%02d",outfolder,gaussian_niter_so[ism_lev_so],gaussian_niter_si[ism_lev_si]);
  FILE *fout=open_text_file_for_output(path);
  
  for(int ispec=0;ispec<nspec;ispec++)
    {
      int ith1=ith_spec[ispec];
      for(int muS_source=0;muS_source<nmuS;muS_source++)
	{
	  //if source is derivative also sink is
	  int start_muS_sink=(muS_source==0)?0:1;
	  int end_muS_sink=  (muS_source==0)?1:4;
	  
	  for(int muS_sink=start_muS_sink;muS_sink<end_muS_sink;muS_sink++)
	    for(int ith2=0;ith2<nthetaS0;ith2++)
	      if(!only_standing_2pts||ith2==ith1)
		for(int im2=0;im2<nmassS0;im2++)
		  {
		    int ip2=ipropS0(ith2,im2,0); //ip2 is not shifted
		    for(int r2=0;r2<2;r2++)
		      if(which_r_S0==2||which_r_S0==r2)
			{
			  if(nch_contr_2pts>0)
			    {
#ifdef POINT_SOURCE_VERSION
			      unsafe_apply_chromo_operator_to_su3spinspin(ch_prop,Pmunu,S0[r2][ip2]);
#else
			      unsafe_apply_chromo_operator_to_colorspinspin(ch_prop,Pmunu,S0[r2][ip2]);
#endif
			    }
			  
			  //decide parameters of mass 1
			  double *mass1=(muS_source==0)?massS0:massS0der;
			  int nmass1=(muS_source==0)?nmassS0:nmassS0der;
			  double *stopping_residues1=(muS_source==0)?stopping_residues_S0:stopping_residues_S0+start_massS0der;
			  
			  for(int im1=0;im1<nmass1;im1++)
			    {
			      int ip1=ipropS0(ith1,im1,muS_source);
			      
			      for(int r1=0;r1<2;r1++)
				if((which_r_S0==2&&(!only_charged_2pts||r2==r1))||which_r_S0==r1)
				  {
				    prop_type *S0_1;
				    
				    //if no derivative on the sink use S0, else derive the sink
				    if(muS_sink==0) S0_1=S0[r1][ip1];
				    else
				      {
					apply_nabla_i(temp_der,S0[r1][ip1],conf,muS_sink);
					S0_1=temp_der;
				      }
				    
				    //header
				    master_fprintf(fout," # m1=%lg th1=%lg res1=%lg r1=%d, m2=%lg th2=%lg res2=%lg r2=%d der_source=%d der_sink=%d",
						   mass1[im1], thetaS0[ith1],stopping_residues1[im1],r1,
						   massS0[im2],thetaS0[ith2],stopping_residues_S0[im2],  r2,
						   muS_source,muS_sink);
				    master_fprintf(fout," smear_source=%d smear_sink=%d\n",gaussian_niter_so[ism_lev_so],gaussian_niter_si[ism_lev_si]);
				    
				    //compute contractions
				    meson_two_points_Wilson_prop(contr_2pts,op1_2pts,S0_1,op2_2pts,S0[r2][ip2],ncontr_2pts);
				    ncontr_tot+=ncontr_2pts;
				    
				    //write 
				    contr_save_time-=take_time();
				    print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,source_coord[0],"",1.0);
				    contr_save_time+=take_time();
				    
				    //if chromo contractions
				    if(nch_contr_2pts>0)
				      {
					//compute them
					meson_two_points_Wilson_prop(ch_contr_2pts,ch_op1_2pts,S0_1,ch_op2_2pts,ch_prop,nch_contr_2pts);
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
    }
  
  if(compute_der) nissa_free(temp_der);
  
  contr_2pts_time+=take_time();
  if(rank==0) fclose(fout);
}

//Calculate and print to file the 3pts
void calculate_all_3pts(int ispec,int ism_lev_so,int ism_lev_se)
{
  char path[1024];
  
  sprintf(path,"%s/3pts_sp%d_%02d_%02d",outfolder,ispec,gaussian_niter_so[ism_lev_so],gaussian_niter_se[ism_lev_se]);
  
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
	    for(int im1=0;im1<contr_3pts_up_to_S0_mass;im1++)
	      {
		int ip1=ipropS0(ith1,im1,0);
		
		//header
		master_fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d,",massS0[im1],thetaS0[ith1],r1,massS1[im2],thetaS1[ith2],r2);
		master_fprintf(fout," smear_source=%d smear_seq=%d\n",gaussian_niter_so[ism_lev_so],gaussian_niter_se[ism_lev_se]);
		
		//compute contractions
		meson_two_points_Wilson_prop(contr_3pts,op1_3pts,S0[r1][ip1],op2_3pts,S1[ip2],ncontr_3pts);
		ncontr_tot+=ncontr_3pts;
		
		//write them
		contr_save_time-=take_time();
		print_contractions_to_file(fout,ncontr_3pts,op1_3pts,op2_3pts,contr_3pts,source_coord[0],"",1.0);
		contr_save_time+=take_time();
		
		//if chromo contractions
		if(nch_contr_3pts>0)
		  {
		    //compute them
		    meson_two_points_Wilson_prop(ch_contr_3pts,ch_op1_3pts,S0[r1][ip1],ch_op2_3pts,ch_prop,nch_contr_3pts);
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
  sprintf(path,"%s/2pts_check_sp%d_%02d_%02d",outfolder,ispec,gaussian_niter_so[ism_lev_so],gaussian_niter_se[ism_lev_se]);
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
  init_nissa(narg,arg);
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
      char fin_file[1024],run_file[1024];
      sprintf(fin_file,"%s/finished",outfolder);
      sprintf(run_file,"%s/running",outfolder);
      file_touch(fin_file);
      rm(run_file);
      
      nanalyzed_conf++;
      enough_time=check_remaining_time();
    }
  
  if(iconf==ngauge_conf) master_printf("Finished all the conf!\n");
  
  tot_time+=take_time();
  close_semileptonic();
  
  return 0;
}
