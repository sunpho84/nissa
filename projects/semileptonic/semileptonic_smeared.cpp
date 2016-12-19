#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nissa.hpp"
#include "driver_corr.hpp"

using namespace nissa;

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
#define PROP_TYPE su3spinspin
#else
#define PROP_TYPE colorspinspin
#endif

//Wilson clover or Tm?
int Wclov_tm;
int rotate_to_phys_basis;
double cSW;

//gauge info
int ngauge_conf,nanalyzed_conf=0;
char conf_path[1024],outfolder[1024];
quad_su3 *conf,*sme_conf;
clover_term_t *Cl;
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
#ifndef POINT_SOURCE_VERSION
int seed;
rnd_t noise_type;
#endif
coords source_coord;
spincolor *source;
PROP_TYPE *original_source;

//smearing parameters
enum conf_smearing_t{no_conf_smearing,ape_conf_smearing,stout_conf_smearing};
conf_smearing_t conf_smearing;
ape_pars_t ape_smearing_pars;
stout_pars_t stout_smearing_pars;
double gaussian_kappa;
int *gaussian_niter_so,nsm_lev_so;
int *gaussian_niter_si,nsm_lev_si;
int *gaussian_niter_se,nsm_lev_se;

//vectors for the S0 props
int compute_der,nmuS;
int npropS0;
int load_S0,save_S0;
PROP_TYPE **S0[2];
int ncgm_solution;
spincolor **cgm_solution,*temp_vec[2];
int ithetaS0_min,ithetaS0_max;

//cgm inverter parameters
double *stop_res_S0;
double *stop_res_S1;
int niter_max=1000000;

//contraction method
int use_new_contraction_layout;
two_pts_comp_t two_pts_comp,three_pts_comp;

//two points contractions
int ncontr_2pts;
int only_standing_2pts;
int only_charged_2pts;
complex *contr_2pts;
double *new_contr_2pts;
int *op_sour_2pts,*op_sink_2pts;

//two points chromo contractions
int nch_contr_2pts;
complex *ch_contr_2pts;
int *ch_op_sour_2pts,*ch_op_sink_2pts;
PROP_TYPE *ch_prop;

//sequential props info
int nspec;
int tsep;
int *imass_spec,*r_spec,*ith_spec;
PROP_TYPE *sequential_source;

//sequential propagators
int npropS1;
PROP_TYPE **S1;

//three points contractions
int contr_3pts_up_to_S0_mass;
int ncontr_3pts;
complex *contr_3pts;
double *new_contr_3pts;
int *op_sour_3pts,*op_sink_3pts;

//two points chromo contractions
int nch_contr_3pts;
complex *ch_contr_3pts;
int *ch_op_sour_3pts,*ch_op_sink_3pts;

//timings
int wall_time;
double tot_prog_time=0;

int ninv_tot=0,ncontr_tot=0;
double inv_time=0,conf_smear_time=0;
double smear_time=0;
double load_save_S0_time=0;
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

//read the parameters to smear the conf entering covariant stuff
void read_conf_smearing_pars()
{
  //read the string
  char conf_smearing_str[100];
  read_str_str("ConfSmearingType",conf_smearing_str,100);
  if(strcasecmp(conf_smearing_str,"no")==0) conf_smearing=no_conf_smearing;
  else
    if(strcasecmp(conf_smearing_str,"ape")==0)
      {
	conf_smearing=ape_conf_smearing;
	read_ape_pars(ape_smearing_pars);
      }
    else
    if(strcasecmp(conf_smearing_str,"stout")==0)
      {
	conf_smearing=stout_conf_smearing;
	read_stout_pars(stout_smearing_pars);
      }
    else crash("Unkown conf smearing type: %s\n",conf_smearing_str);
}

//This function contract a source with a sequential spinor putting the passed list of operators
void contract_with_source(complex *corr,PROP_TYPE *S1,int *list_op,PROP_TYPE *source)
{
  //Temporary vector for the internal matrices
  dirac_matr t1[ncontr_2pts],t2[ncontr_2pts];
  complex *loc_corr=nissa_malloc("loc_corr",ncontr_2pts*glb_size[0],complex);
  
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Init the second
      t1[icontr]=base_gamma[0];
      t2[icontr]=base_gamma[list_op[icontr]];
    }
  
  //Call the routine which does the real contraction
#ifdef POINT_SOURCE_VERSION
  trace_g_ccss_dag_g_ccss(corr,loc_corr,t1,S1,t2,source,ncontr_2pts);
#else
  trace_g_css_dag_g_css(corr,loc_corr,t1,S1,t2,source,ncontr_2pts);
#endif
  
  nissa_free(loc_corr);
}

//generate the source
void generate_source()
{
#ifdef POINT_SOURCE_VERSION
  generate_delta_source(original_source,source_coord);
#else
  generate_spindiluted_source(original_source,noise_type,source_coord[0]);
  
  //if asked, save the source
  if(save_source)
    {
      char outpath[1024];
      sprintf(outpath,"%s/source",outfolder);
      write_double_vector(outpath,original_source,64,"source");
    }
#endif
}

//Generate a sequential source for S1
void generate_sequential_source(int ispec)
{
  int r=r_spec[ispec];
  
  master_printf("\nCreating the sequential source for spectator %d\n",ispec);
  NISSA_LOC_VOL_LOOP(ivol)
    {
      //put to zero everywhere but on the slice
      if(glb_coord_of_loclx[ivol][0]!=(source_coord[0]+tsep)%glb_size[0])
	memset(sequential_source[ivol],0,sizeof(PROP_TYPE));
      else
	{
	  //avoid to put g5, beacuse commute with (i+-g5)/sqrt(2) and cancel with those of the QQ
	  memcpy(sequential_source[ivol],S0[r][ipropS0(ith_spec[ispec],imass_spec[ispec],0)][ivol],sizeof(PROP_TYPE));
	  if(Wclov_tm and rotate_to_phys_basis) //if doing tm and want to rotate
	    for(int c=0;c<NCOL;c++) //rotate as r because it's D^-1
	      {
#ifdef POINT_SOURCE_VERSION
		for(int c1=0;c1<NCOL;c1++) rotate_spinspin_to_physical_basis(sequential_source[ivol][c][c1],r,r);
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
  //Decide if twisted(1) or clover run
  read_str_int("TwistedMassRun",&Wclov_tm);
  //Decide if to rotate to phsyical basis
  if(Wclov_tm)read_str_int("RotateToPhysBasis",&rotate_to_phys_basis);
  //Kappa is read really only for tm
  if(Wclov_tm) read_str_double("Kappa",&kappa);
  read_str_double("cSW",&cSW);
  
  // 2) Read information about the source
  
#ifndef POINT_SOURCE_VERSION
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  //Read the noise type
  char noise_type_str[20];
  read_str_str("NoiseType",noise_type_str,20);
  noise_type=convert_str_to_rnd_t(noise_type_str);
#endif
  //read whether we want to save the source
  read_str_int("SaveSource",&save_source);
  
  // 3) Smearing parameters
  
  //Smearing parameters
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
  read_conf_smearing_pars();
  
  // 4) Read list of masses and of thetas
  
  //the number or fr can be chosen only if tm, otherwise use 0
  if(Wclov_tm==1) read_str_int("WhichRS0",&which_r_S0);
  else which_r_S0=0;
  read_str_int("UseCgmS0",&use_cgm_S0);
  
  read_list_of_double_pairs("MassResiduesS0",&nmassS0,&massS0,&stop_res_S0);
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
  read_str_int("UseNewContractionLayout",&use_new_contraction_layout);
  if(use_new_contraction_layout)
    {
      char path[100];
      read_str_str("TwoPointsOpListFilePath",path,100);
      two_pts_comp=read_two_pts_sink_source_corr_from_file(path);
      ncontr_2pts=two_pts_comp.ncorr;
      verbosity_lv2_master_printf("Read %d corrs, corresponding to %d contr\n",ncontr_2pts,two_pts_comp.size());
      for(int icorr=0;icorr<two_pts_comp.ncorr;icorr++)
	verbosity_lv2_master_printf(" %s\n",two_pts_comp.corr_name[icorr].c_str());
      new_contr_2pts=nissa_malloc("new_contr_2pts",ncontr_2pts*glb_size[0],double);
    }
  else
    {
      read_str_int("NContrTwoPoints",&ncontr_2pts);
      contr_2pts=nissa_malloc("contr_2pts",ncontr_2pts*glb_size[0],complex);
      op_sour_2pts=nissa_malloc("op_sour_2pts",ncontr_2pts,int);
      op_sink_2pts=nissa_malloc("op_sink_2pts",ncontr_2pts,int);
      for(int icontr=0;icontr<ncontr_2pts;icontr++)
	{
	  //Read the operator pairs
	  read_int(&(op_sour_2pts[icontr]));
	  read_int(&(op_sink_2pts[icontr]));
	  
	  master_printf(" contr.%d %d %d\n",icontr,op_sour_2pts[icontr],op_sink_2pts[icontr]);
	}
      
      read_str_int("NChromoContrTwoPoints",&nch_contr_2pts);
      ch_contr_2pts=nissa_malloc("ch_contr_2pts",nch_contr_2pts*glb_size[0],complex);
      ch_op_sour_2pts=nissa_malloc("ch_op_sour_2pts",ncontr_2pts,int);
      ch_op_sink_2pts=nissa_malloc("ch_op_sink_2pts",ncontr_2pts,int);
      for(int icontr=0;icontr<nch_contr_2pts;icontr++)
	{
	  //Read the operator pairs
	  read_int(&(ch_op_sour_2pts[icontr]));
	  read_int(&(ch_op_sink_2pts[icontr]));
	  
	  master_printf(" ch-contr.%d %d %d\n",icontr,ch_op_sour_2pts[icontr],ch_op_sink_2pts[icontr]);
	}
    }
  
  // 6) Read list of masses and of thetas for S1
  
  if(Wclov_tm==1) read_str_int("UseCgmS1",&use_cgm_S1);
  read_list_of_double_pairs("MassResiduesS1",&nmassS1,&massS1,&stop_res_S1);
  read_list_of_doubles("NThetaS1",&nthetaS1,&thetaS1);
  read_str_int("ContrThreePointsUpToS0Mass",&contr_3pts_up_to_S0_mass);
  
  // 7) three points functions
  
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
      
      if(ith_spec[ispec]<0 or ith_spec[ispec]>=nthetaS0)    crash("theta for ispec %d out of bounds",ispec);
      if(imass_spec[ispec]<0 or imass_spec[ispec]>=nmassS0) crash("mass for ispec %d out of bounds",ispec);
      if(r_spec[ispec]<0 or r_spec[ispec]>=2)               crash("r for ispec %d out of bounds",ispec);
      if(which_r_S0!=2 and r_spec[ispec]!=which_r_S0)        crash("r for ispec %d uncomputed",ispec);
      
      master_printf(" spec %d: th=%g, m=%g, r=%d\n",ispec,thetaS0[ith_spec[ispec]],massS0[imass_spec[ispec]],r_spec[ispec]);
    }
  
  if(use_new_contraction_layout)
    {
      char path[100];
      read_str_str("ThreePointsOpListFilePath",path,100);
      three_pts_comp=read_two_pts_sink_source_corr_from_file(path);
      ncontr_3pts=three_pts_comp.ncorr;
      verbosity_lv2_master_printf("Read %d corrs, corresponding to %d contr\n",ncontr_3pts,three_pts_comp.size());      
      for(int icorr=0;icorr<three_pts_comp.ncorr;icorr++)
	verbosity_lv2_master_printf(" %s\n",three_pts_comp.corr_name[icorr].c_str());
      
      new_contr_3pts=nissa_malloc("new_contr_3pts",ncontr_3pts*glb_size[0],double);
    }
  else
    {
      read_str_int("NContrThreePoints",&ncontr_3pts);
      contr_3pts=nissa_malloc("contr_3pts",ncontr_3pts*glb_size[0],complex); 
      
      op_sour_3pts=nissa_malloc("op_sour_3pts",ncontr_3pts,int);
      op_sink_3pts=nissa_malloc("op_sink_3pts",ncontr_3pts,int);
      for(int icontr=0;icontr<ncontr_3pts;icontr++)
	{
	  //Read the operator pairs
	  read_int(&(op_sour_3pts[icontr]));
	  read_int(&(op_sink_3pts[icontr]));
	  
	  master_printf(" contr.%d %d %d\n",icontr,op_sour_3pts[icontr],op_sink_3pts[icontr]);
	}
      
      read_str_int("NChromoContrThreePoints",&nch_contr_3pts);
      ch_contr_3pts=nissa_malloc("ch_contr_3pts",nch_contr_3pts*glb_size[0],complex);
      ch_op_sour_3pts=nissa_malloc("ch_op_sour_3pts",nch_contr_3pts,int);
      ch_op_sink_3pts=nissa_malloc("ch_op_sink_3pts",nch_contr_3pts,int);
      for(int icontr=0;icontr<nch_contr_3pts;icontr++)
	{
	  //Read the operator pairs
	  read_int(&(ch_op_sour_3pts[icontr]));
	  read_int(&(ch_op_sink_3pts[icontr]));
	  
	  master_printf(" ch-contr.%d %d %d\n",icontr,ch_op_sour_3pts[icontr],ch_op_sink_3pts[icontr]);
	}
    }
  if(ncontr_3pts!=0 or nch_contr_3pts!=0) sequential_source=nissa_malloc("Sequential source",loc_vol,PROP_TYPE);
  
  read_str_int("NGaugeConf",&ngauge_conf);
  
  master_printf("\n");
  
  ////////////////////////////////////// end of input reading/////////////////////////////////
  
  //allocate gauge conf, Cl and all the needed spincolor and propagators
  conf=nissa_malloc("or_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  if(conf_smearing!=no_conf_smearing) sme_conf=nissa_malloc("sm_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  else sme_conf=conf;
  Cl=nissa_malloc("Cl",loc_vol,clover_term_t);
  
  //Allocate all the S0 PROP_TYPE vectors
  npropS0=nthetaS0*nmassS0;
  if(compute_der) npropS0+=3*nthetaS0*nmassS0der;
  S0[0]=nissa_malloc("S0[0]",npropS0,PROP_TYPE*);
  S0[1]=nissa_malloc("S0[1]",npropS0,PROP_TYPE*);
  for(int iprop=0;iprop<npropS0;iprop++)
    for(int r=0;r<2;r++)
      if(which_r_S0==2 or which_r_S0==r)
	S0[r][iprop]=nissa_malloc("S0[r]",loc_vol+bord_vol,PROP_TYPE);
  
  //Allocate nmass spincolors, for the cgm solutions
  ncgm_solution=std::max(nmassS0,nmassS1);
  cgm_solution=nissa_malloc("cgm_solution",ncgm_solution,spincolor*);
  for(int imass=0;imass<ncgm_solution;imass++) cgm_solution[imass]=nissa_malloc("cgm_solution",loc_vol+bord_vol,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol,spincolor);
  
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  original_source=nissa_malloc("original_source",loc_vol,PROP_TYPE);
  
  //Allocate one PROP_TYPE for the chromo-contractions
  if(nch_contr_2pts!=0 or nch_contr_3pts!=0) ch_prop=nissa_malloc("chromo-prop",loc_vol,PROP_TYPE);
  
  //Allocate all the S1 PROP_TYPE vectors
  npropS1=nthetaS1*nmassS1;
  S1=nissa_malloc("S1",npropS1,PROP_TYPE*);
  for(int iprop=0;iprop<npropS1;iprop++) S1[iprop]=nissa_malloc("S1[i]",loc_vol,PROP_TYPE);
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
      ok_conf=!(file_exists(fin_file)) and !(file_exists(run_file));
      
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
  while(!ok_conf and (*iconf)<ngauge_conf);
  
  master_printf("\n");
  
  return ok_conf;
}

//read the conf and setup it
void setup_conf()
{
  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_ildg_gauge_conf(conf,conf_path);
  master_printf("plaq: %16.16lg\n",global_plaquette_lx_conf(conf));
  
  clover_term(Cl,cSW,conf);
  
  conf_smear_time-=take_time();
  
  //prepare the smerded version and compute plaquette
  switch(conf_smearing)
    {
    case no_conf_smearing:
      break;
    case ape_conf_smearing:
      ape_spatial_smear_conf(sme_conf,conf,ape_smearing_pars.alpha,ape_smearing_pars.nlevels);
      break;
    case stout_conf_smearing:
      //workaround because stout is implemented only for eo conf
      quad_su3 *eo_conf[2];
      eo_conf[EVN]=nissa_malloc("new_conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
      eo_conf[ODD]=nissa_malloc("new_conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
      split_lx_vector_into_eo_parts(eo_conf,conf);
      stout_smear(eo_conf,eo_conf,&(stout_smearing_pars));
      paste_eo_parts_into_lx_vector(sme_conf,eo_conf);
      nissa_free(eo_conf[EVN]);
      nissa_free(eo_conf[ODD]);
      break;
    default:
      crash("unknown conf smearing type %d",(int)conf_smearing);
      break;
    }
  
  conf_smear_time+=take_time();
  
  if(conf_smearing!=no_conf_smearing) master_printf("smerded plaq: %.16g\n",global_plaquette_lx_conf(sme_conf));
  
  //put the anti-periodic condition on the temporal border
  old_theta[0]=old_theta[1]=old_theta[2]=old_theta[3]=0;
  put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
  adapt_theta(conf,old_theta,put_theta,1,1);
}

//Finalization
void close_semileptonic()
{
  close_input();
  
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g",tot_prog_time);
  
  double contr_time=contr_2pts_time+contr_3pts_time;
  master_printf(", of which:\n");
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf("  of which  %02.2f%s for %d cgm inversion overhead (%2.2gs avg)\n",cgm_inv_over_time/inv_time*100,"%",
		ninv_tot,cgm_inv_over_time/ninv_tot);
  master_printf(" - %02.2f%s to smear configuration\n",conf_smear_time*100.0/tot_prog_time,"%");
  master_printf(" - %02.2f%s to sink-smear propagators\n",smear_time*100.0/tot_prog_time,"%");
  master_printf(" - %02.2f%s to load or save propagators\n",load_save_S0_time*100.0/tot_prog_time,"%");
  master_printf(" - %02.2f%s to perform %d contr. (%2.2gs avg) of which:\n",contr_time/tot_prog_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
  master_printf("   * %02.2f%s to compute two points\n",contr_2pts_time*100.0/contr_time,"%");
  master_printf("   * %02.2f%s to compute three points\n",contr_3pts_time*100.0/contr_time,"%");
  master_printf(" - %02.2f%s to save correlations\n",contr_save_time*100.0/tot_prog_time,"%");
  nissa_free(Cl);nissa_free(conf);if(conf_smearing!=no_conf_smearing) nissa_free(sme_conf);
  for(int iprop=0;iprop<npropS0;iprop++)
    for(int r=0;r<2;r++)
      if(which_r_S0==2 or which_r_S0==r) nissa_free(S0[r][iprop]);
  for(int iprop=0;iprop<npropS1;iprop++) nissa_free(S1[iprop]);
  nissa_free(S0[0]);nissa_free(S0[1]);nissa_free(S1);
  nissa_free(temp_vec[0]);nissa_free(temp_vec[1]);
  if(nch_contr_2pts!=0 or nch_contr_3pts!=0) nissa_free(ch_prop);
  if(ncontr_3pts!=0 or nch_contr_3pts!=0) nissa_free(sequential_source);
  if(!use_new_contraction_layout)
    {
      nissa_free(contr_2pts);nissa_free(ch_contr_2pts);
      nissa_free(contr_3pts);nissa_free(ch_contr_3pts);
      nissa_free(op_sour_2pts);nissa_free(op_sink_2pts);
      nissa_free(op_sour_3pts);nissa_free(op_sink_3pts);
      nissa_free(ch_op_sour_2pts);nissa_free(ch_op_sink_2pts);
      nissa_free(ch_op_sour_3pts);nissa_free(ch_op_sink_3pts);
    }
  else
    {
      nissa_free(new_contr_2pts);
      nissa_free(new_contr_3pts);
    }
  nissa_free(ith_spec);nissa_free(r_spec);nissa_free(imass_spec);
  for(int imass=0;imass<ncgm_solution;imass++) nissa_free(cgm_solution[imass]);
  nissa_free(cgm_solution);
  nissa_free(source);nissa_free(original_source);
}

//smear addditivily a propagator
void smear_additive_propagator(PROP_TYPE *out,PROP_TYPE *in,int ism_lev,int *gaussian_niter)
{
  spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  
  int nsme=gaussian_niter[ism_lev];
  if(ism_lev>0) nsme-=gaussian_niter[ism_lev-1];
  
  //loop over dirac index
#ifdef POINT_SOURCE_VERSION
  for(int ic=0;ic<NCOL;ic++)
#endif
    for(int id=0;id<4;id++)
      {
#ifdef POINT_SOURCE_VERSION
	get_spincolor_from_su3spinspin(temp,in,id,ic);
#else
	get_spincolor_from_colorspinspin(temp,in,id);
#endif
	
	gaussian_smearing(temp,temp,sme_conf,gaussian_kappa,nsme);
	
#ifdef POINT_SOURCE_VERSION
	put_spincolor_into_su3spinspin(out,temp,id,ic);
#else
	put_spincolor_into_colorspinspin(out,temp,id);
#endif
      }
  
  nissa_free(temp);
}

//calculate the standard propagators
void calculate_all_S0(int ism_lev_so)
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
      double *stop_res=(muS==0)?stop_res_S0:stop_res_S0+start_massS0der;
      
      for(int itheta=0;itheta<nthetaS0;itheta++)
	{
	  //adapt the border condition
	  put_theta[1]=put_theta[2]=put_theta[3]=thetaS0[itheta];
	  adapt_theta(conf,old_theta,put_theta,1,1);
	  
	  //loop over the source dirac index
#ifdef POINT_SOURCE_VERSION
	  for(int ic=0;ic<NCOL;ic++)
#endif
	    for(int id=0;id<NDIRAC;id++)
	      { 
		//put the g5
#ifdef POINT_SOURCE_VERSION
		get_spincolor_from_su3spinspin(source,original_source,id,ic);
#else
		get_spincolor_from_colorspinspin(source,original_source,id);
#endif
		//add gamma5 apart if using cg or tm
		if(!Wclov_tm or use_cgm_S0 or (Wclov_tm and cSW!=0)) safe_dirac_prod_spincolor(source,base_gamma+5,source);
		
		//if needed apply nabla
		if(muS>0) apply_nabla_i(source,source,sme_conf,muS);
		
		//invert
		if(!load_S0)
		  {
		    double part_time=-take_time();
		    
		    //decide if to use multimass or single mass
		    if(use_cgm_S0)
		      {
			if(cSW==0) inv_tmQ2_cgm(cgm_solution,conf,kappa,mass,nmass,niter_max,stop_res,source);
			else       inv_tmclovQ2_cgm(cgm_solution,conf,kappa,Cl,mass,nmass,niter_max,stop_res,source);
		      }
		    else
		      for(int imass=0;imass<nmass;imass++)
			{
			  //the sign of mass depends on r
			  double m=mass[imass];
			  if(Wclov_tm)
			    {
			      if(which_r_S0==0) m*=-1;
			      if(cSW==0) inv_tmD_cg_eoprec(cgm_solution[imass],NULL,conf,kappa,m,niter_max,stop_res[imass],source);
			      else inv_tmclovQ_cg(cgm_solution[imass],NULL,conf,kappa,Cl,m,niter_max,stop_res[imass],source);
			    }
			  else //m=kappa
			    inv_WclovQ_cg(cgm_solution[imass],NULL,conf,m,Cl,niter_max,stop_res[imass],source);
			  
			  master_printf("Finished submass[%d]=%lg\n",imass,m);
			}
		    
		    part_time+=take_time();ninv_tot++;inv_time+=part_time;
		    
		    master_printf("Finished the inversion of S0 theta %d, ",itheta);
		    if(compute_der) master_printf("source derivative %d ",muS);

#ifdef POINT_SOURCE_VERSION
		    master_printf("color index %d ",ic);
#endif
		    master_printf("dirac index %d",id);
		    
		    master_printf(", in %lg s\n",part_time);
		  }
		
		//read or write, if needed
		if(save_S0 or load_S0)
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
		      
		      load_save_S0_time-=take_time();
		      
		      if(save_S0) write_double_vector(path,cgm_solution[imass],64,"S0");
		      else        read_real_vector(cgm_solution[imass],path,"S0");
		      
		      load_save_S0_time+=take_time();
		    }
		
		//reconstruct the doublet
		for(int imass=0;imass<nmass;imass++)
		  {
		    if(Wclov_tm and use_cgm_S0)
		      {
			if(cSW==0) reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass[imass],cgm_solution[imass]);
			else       reconstruct_tmclov_doublet(temp_vec[0],temp_vec[1],conf,kappa,Cl,mass[imass],cgm_solution[imass]);
			master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
		      }
		    else memcpy(temp_vec[which_r_S0],cgm_solution[imass],sizeof(spincolor)*loc_vol);
		    for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
		      if(which_r_S0==r or which_r_S0==2)
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
  if(Wclov_tm and rotate_to_phys_basis)
    {
      for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
	if(which_r_S0==r or which_r_S0==2)
	  for(int ipropS0=0;ipropS0<npropS0;ipropS0++) //put the (1+ig5)/sqrt(2) factor
	    {
#ifdef POINT_SOURCE_VERSION
	      rotate_vol_su3spinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
#else
	      rotate_vol_colorspinspin_to_physical_basis(S0[r][ipropS0],!r,!r);
#endif
	    }
      master_printf("Propagators rotated\n");
    }
  
  master_printf("\n");
}

//calculate the sequential propagators
void calculate_all_S1(int ispec,int ism_lev_se)
{
  //smear additively the seq
  master_printf("\nSeq Smearing level: %d (will be applied twice!)\n",ism_lev_se);
  smear_additive_propagator(sequential_source,sequential_source,ism_lev_se,gaussian_niter_se);
  smear_additive_propagator(sequential_source,sequential_source,ism_lev_se,gaussian_niter_se);
  master_printf("\n");
  
  for(int itheta=0;itheta<nthetaS1;itheta++)
    {
      //adapt the border condition
      put_theta[1]=put_theta[2]=put_theta[3]=thetaS1[itheta];
      adapt_theta(conf,old_theta,put_theta,1,1);
      
      //loop over seq
#ifdef POINT_SOURCE_VERSION
      for(int ic=0;ic<NCOL;ic++)
#endif
	for(int id=0;id<NDIRAC;id++)
	  { 
#ifdef POINT_SOURCE_VERSION
	    get_spincolor_from_su3spinspin(source,sequential_source,id,ic);
#else
	    get_spincolor_from_colorspinspin(source,sequential_source,id);
#endif
	    safe_dirac_prod_spincolor(source,base_gamma+5,source);
	    
	    //if inverting Q
	    if(!Wclov_tm or use_cgm_S1 or (Wclov_tm and cSW!=0)) safe_dirac_prod_spincolor(source,base_gamma+5,source); 
	    
	    double part_time=-take_time();
	    
	    //decide to use one or the other inverters
	    if(use_cgm_S1)
	      {
		if(cSW==0) inv_tmQ2_cgm(cgm_solution,conf,kappa,massS1,nmassS1,niter_max,stop_res_S1,source);
		else inv_tmclovQ2_cgm(cgm_solution,conf,kappa,Cl,massS1,nmassS1,niter_max,stop_res_S1,source);
	      }
	    else
	      for(int imass=0;imass<nmassS1;imass++)
		{
		  //since r0==0 implicates r1=1, mass=-mass when r0=1
		  double m=massS1[imass];
		  if(Wclov_tm)
		    {
		      if(r_spec[ispec]==1) m*=-1;
		      if(cSW==0) inv_tmD_cg_eoprec(cgm_solution[imass],NULL,conf,kappa,m,niter_max,stop_res_S1[imass],source);
		      else inv_tmclovQ_cg(cgm_solution[imass],NULL,conf,kappa,Cl,m,niter_max,stop_res_S1[imass],source);
		    }
		  else inv_WclovQ_cg(cgm_solution[imass],NULL,conf,m,Cl,niter_max,stop_res_S1[imass],source);
		}
	    
	    part_time+=take_time();ninv_tot++;inv_time+=part_time;
	    
	    master_printf("Finished the inversion of S1 theta %d, seq sme lev %d,",itheta,ism_lev_se);
	    
#ifdef POINT_SOURCE_VERSION
	    master_printf(" color index %d,",ic);
#endif
	    master_printf(" dirac index %d,",id);
	    master_printf(" in %g sec\n",part_time);
	    
	    for(int imass=0;imass<nmassS1;imass++)
	      {
		//reconstruct the doublet: r(S1)=!r(spec), so we have to multiply by Q+ if r(spec)==1 and Q- if 0
		double reco_mass=-massS1[imass];
		if(r_spec[ispec]==1) reco_mass=-reco_mass;
		//use temp_vec[0] as temporary storage
		if(Wclov_tm and use_cgm_S1)
		  {
		    if(cSW==0) apply_tmQ(temp_vec[0],conf,kappa,reco_mass,cgm_solution[imass]);
		    else       apply_tmclovQ(temp_vec[0],conf,kappa,Cl,reco_mass,cgm_solution[imass]);
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
  if(Wclov_tm and rotate_to_phys_basis)
    {
      for(int ipropS1=0;ipropS1<npropS1;ipropS1++) //but, being D^-1, everything is swapped
	{
#ifdef POINT_SOURCE_VERSION
	  rotate_vol_su3spinspin_to_physical_basis(S1[ipropS1],!(r_spec[ispec]),(r_spec[ispec]));
#else
	  rotate_vol_colorspinspin_to_physical_basis(S1[ipropS1],!(r_spec[ispec]),(r_spec[ispec]));
#endif
	}
      master_printf("Propagators rotated\n");
    }
  
  master_printf("\n");
}

//transpose dirac indices with spatial ones, so that destination
//is arranged as: t,id_si,id_so,ri,spat,color
THREADABLE_FUNCTION_2ARG(prepare_prop_for_new_contraction, PROP_TYPE*,prop, PROP_TYPE*,aux)
{
  GET_THREAD_ID();
  
  int spat_color_vol=loc_vol/loc_size[0]*sizeof(PROP_TYPE)/sizeof(spinspin);
  
  //fill aux
  for(int t=0;t<loc_size[0];t++)
    NISSA_PARALLEL_LOOP(ispat_color,0,spat_color_vol)
      for(int idirac_ri=0;idirac_ri<32;idirac_ri++)
	((double*)aux)[ispat_color+spat_color_vol*(idirac_ri+32*t)]=
	  ((double*)prop)[idirac_ri+32*(ispat_color+spat_color_vol*t)];
  THREAD_BARRIER();
  
  //copy aux to prop
  parallel_memcpy(prop,aux,sizeof(PROP_TYPE)*loc_vol);
}}
//do the opposite
THREADABLE_FUNCTION_2ARG(revert_prop_from_new_contraction, PROP_TYPE*,prop, PROP_TYPE*,aux)
{
  GET_THREAD_ID();
  
  int spat_color_vol=loc_vol/loc_size[0]*sizeof(PROP_TYPE)/sizeof(spinspin);
  
  //fill aux
  for(int t=0;t<loc_size[0];t++)
    for(int idirac_ri=0;idirac_ri<32;idirac_ri++)
      NISSA_PARALLEL_LOOP(ispat_color,0,spat_color_vol)
	((double*)aux)[idirac_ri+32*(ispat_color+spat_color_vol*t)]=
	((double*)prop)[ispat_color+spat_color_vol*(idirac_ri+32*t)];
  THREAD_BARRIER();
  
  //copy aux to prop
  parallel_memcpy(prop,aux,sizeof(PROP_TYPE)*loc_vol);
}}

//Uses the new layout
THREADABLE_FUNCTION_5ARG(new_meson_two_points, double*,glb_2pts, double*,loc_2pts, PROP_TYPE*,S_back, PROP_TYPE*,S_forw, two_pts_comp_t*,comp)
{
  GET_THREAD_ID();
  
  vector_reset(loc_2pts);
  
  //contract
  int slice_vol=loc_vol/loc_size[0]*sizeof(PROP_TYPE)/sizeof(spinspin);
  comp->summ_the_loc_forw_back_contractions(loc_2pts,(double*)S_forw,(double*)S_back,slice_vol,source_coord[0]);
  
  THREAD_BARRIER();
  
  //reduce
  if(IS_MASTER_THREAD) MPI_Reduce(loc_2pts,glb_2pts,glb_size[0]*comp->ncorr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  THREAD_BARRIER();
}}

//Calculate and print to file the 2pts
void calculate_all_2pts(int ism_lev_so,int ism_lev_si)
{
  PROP_TYPE *temp_der1=(compute_der>=1)?nissa_malloc("temp_der1",loc_vol+bord_vol,PROP_TYPE):NULL;
  PROP_TYPE *temp_der2=(compute_der>=2)?nissa_malloc("temp_der2",loc_vol+bord_vol,PROP_TYPE):NULL;
  PROP_TYPE *temp_transp=use_new_contraction_layout?nissa_malloc("temp_trans",loc_vol,PROP_TYPE):NULL;
  complex *loc_2pts;
  double *new_loc_2pts;
  if(!use_new_contraction_layout) loc_2pts=nissa_malloc("contr_2pts",std::max(ncontr_2pts,nch_contr_2pts)*glb_size[0],complex);
  else new_loc_2pts=nissa_malloc("contr_2pts",ncontr_2pts*glb_size[0],double);
  
  //smear additively the propagators
  smear_time-=take_time();
  
  for(int r=0;r<2;r++)
    if(which_r_S0==2 or which_r_S0==r)
      for(int iprop=0;iprop<npropS0;iprop++)
	smear_additive_propagator(S0[r][iprop],S0[r][iprop],ism_lev_si,gaussian_niter_si);
  
  //take intermediate timing
  double temp_time=take_time();
  smear_time+=temp_time;
  contr_2pts_time-=temp_time;
  
  //open output file
  char path[1024];
  sprintf(path,"%s/2pts_%02d_%02d",outfolder,gaussian_niter_so[ism_lev_so],gaussian_niter_si[ism_lev_si]);
  FILE *fout=open_text_file_for_output(path);
  
  //choose if to derive or not
  int nmuS1=(compute_der>=1)?4:1;
  int nmuS2=(compute_der>=2)?4:1;
  
  //change to new layout
  if(use_new_contraction_layout)
    for(int r=0;r<2;r++)
      if(which_r_S0==2 or which_r_S0==r)
	for(int iprop=0;iprop<npropS0;iprop++)
	  prepare_prop_for_new_contraction(S0[r][iprop],temp_transp);
  
  //loop on spectators (fixed theta)
  for(int ispec=0;ispec<nspec;ispec++)
    {
      int ith1=ith_spec[ispec];
      //loop on derivative at source and sink of second propagator
      for(int muS_source2=0;muS_source2<nmuS2;muS_source2++)
	for(int muS_sink2=0;muS_sink2<nmuS2;muS_sink2++)
	  //loop on derivative at source and sink of first propagator
	  for(int muS_source1=(muS_source2==0)?0:muS_source2+1;muS_source1<nmuS1;muS_source1++)
	    for(int muS_sink1=0;muS_sink1<nmuS1;muS_sink1++)
	      //loop on theta of second propagator
	      for(int ith2=0;ith2<nthetaS0;ith2++)
		if(!only_standing_2pts or ith2==ith1)
		  {
		    //decide parameters of mass 2
		    double *mass2=(muS_source2==0)?massS0:massS0der;
		    int nmass2=(muS_source2==0)?nmassS0:nmassS0der;
		    double *stop_res2=(muS_source2==0)?stop_res_S0:stop_res_S0+start_massS0der;
		    
		    //loop on mass of second propagator
		    for(int im2=0;im2<nmass2;im2++)
		      {
			int ip2=ipropS0(ith2,im2,muS_source2);
			//loop on r of second propagator
			for(int r2=0;r2<2;r2++)
			  if(which_r_S0==2 or which_r_S0==r2)
			    {
			      //if no derivative on the sink use S0, else derive the sink
			      PROP_TYPE *S0_2=(muS_sink2==0)?S0[r2][ip2]:temp_der2;
			      if(muS_sink2!=0)
				{
				  //in case revert from new layout
				  if(use_new_contraction_layout)
				    revert_prop_from_new_contraction(S0[r2][ip2],temp_transp);
				  
				  //derive the sink of second propagator
				  apply_nabla_i(temp_der2,S0[r2][ip2],sme_conf,muS_sink2);
				  
				  //and return back to new layout
				  if(use_new_contraction_layout)
                                    {
                                      prepare_prop_for_new_contraction(S0[r2][ip2],temp_transp);
                                      prepare_prop_for_new_contraction(temp_der2,temp_transp);
				    }
				}
			      
			      //apply chromo operator
			      if(nch_contr_2pts>0)
				{
				  //in case revert from new layout
				  if(use_new_contraction_layout)
				    revert_prop_from_new_contraction(S0_2,temp_transp);
				  chromo_operator_remove_cSW(Cl,cSW);
#ifdef POINT_SOURCE_VERSION
				  unsafe_apply_chromo_operator_to_su3spinspin(ch_prop,Cl,S0_2);
#else
				  unsafe_apply_chromo_operator_to_colorspinspin(ch_prop,Cl,S0_2);
#endif
				  chromo_operator_include_cSW(Cl,cSW);
				  //and return back to new layout
				  if(use_new_contraction_layout)
				    {
				      prepare_prop_for_new_contraction(S0_2,temp_transp);
				      prepare_prop_for_new_contraction(ch_prop,temp_transp);
				    }
				}
			      
			      //decide parameters of mass 1
			      double *mass1=(muS_source1==0)?massS0:massS0der;
			      int nmass1=(muS_source1==0)?nmassS0:nmassS0der;
			      double *stop_res1=(muS_source1==0)?stop_res_S0:stop_res_S0+start_massS0der;
			      
			      //loop on mass of first propagator
			      for(int im1=0;im1<nmass1;im1++)
				{
				  int ip1=ipropS0(ith1,im1,muS_source1);
				  //loop on r of first propagator
				  for(int r1=0;r1<2;r1++)
				    if((which_r_S0==2 and(!only_charged_2pts or r2==r1)) or which_r_S0==r1)
				      {
					//if no derivative on the sink use S0, else derive the sink
					PROP_TYPE *S0_1=(muS_sink1==0)?S0[r1][ip1]:temp_der1;
					if(muS_sink1!=0)
					  {
					    //in case revert from new layout
					    if(use_new_contraction_layout)
					      revert_prop_from_new_contraction(S0[r1][ip1],temp_transp);
					    
					    //derive the sink of first propagator
					    apply_nabla_i(temp_der1,S0[r1][ip1],sme_conf,muS_sink1);
					    
					    //and come back
					    if(use_new_contraction_layout)
					      {
						prepare_prop_for_new_contraction(S0[r1][ip1],temp_transp);
						prepare_prop_for_new_contraction(temp_der1,temp_transp);
					      }
					  }
					
					//header
					master_fprintf(fout," # m1=%lg th1=%lg res1=%lg r1=%d,"
						       " m2=%lg th2=%lg res2=%lg r2=%d,",
						       mass1[im1],thetaS0[ith1],stop_res1[im1],r1,
						       mass2[im2],thetaS0[ith2],stop_res2[im2],r2);
					master_fprintf(fout," dsrc1=%d dsrc2=%d, dsnk1=%d dsnk2=%d,",
						       muS_source1,muS_source2,muS_sink1,muS_sink2);
					master_fprintf(fout," sm_src=%d sm_snk=%d\n",
						       gaussian_niter_so[ism_lev_so],gaussian_niter_si[ism_lev_si]);
					
					//compute contractions
					if(use_new_contraction_layout)
					  new_meson_two_points(new_contr_2pts,new_loc_2pts,S0_1,S0_2,&two_pts_comp);
					 else meson_two_points_Wilson_prop(contr_2pts,loc_2pts,op_sour_2pts,S0_1,
									   op_sink_2pts,S0_2,ncontr_2pts);
					
					
					//write
					ncontr_tot+=ncontr_2pts;
					contr_save_time-=take_time();
					if(use_new_contraction_layout)
					  two_pts_comp.print_correlations_to_file(fout,new_contr_2pts);
					 else print_contractions_to_file(
				         fout,ncontr_2pts,op_sour_2pts,op_sink_2pts,contr_2pts,source_coord[0],"",1.0);
					contr_save_time+=take_time();
					//if chromo contractions
					/*
					if(nch_contr_2pts>0)
					  {
					    if(use_new_contraction_layout)
					      new_meson_two_points(ch_contr_2pts,loc_2pts,ch_op_sour_2pts,S0_1,
									 ch_op_sink_2pts,ch_prop,nch_contr_2pts);
					    else meson_two_points_Wilson_prop(ch_contr_2pts,loc_2pts,ch_op_sour_2pts,S0_1,
									      ch_op_sink_2pts,ch_prop,nch_contr_2pts);
					    
					    //write
#ifdef BENCH
					    ncontr_tot+=nch_contr_2pts;
					    contr_save_time-=take_time();
#endif
					    (use_new_contraction_layout?
					     print_optimized_contractions_to_file:print_contractions_to_file)
					      (fout,nch_contr_2pts,ch_op_sour_2pts,ch_op_sink_2pts,ch_contr_2pts,
					       source_coord[0],"CHROMO-",1.0);
#ifdef BENCH
					    contr_save_time+=take_time();
#endif
					  }
					*/
					master_fprintf(fout,"\n");
				      }
				}
			    }
		      }
		  }
    }
  
  //return to non new layout
  if(use_new_contraction_layout)
    for(int r=0;r<2;r++)
      if(which_r_S0==2 or which_r_S0==r)
	for(int iprop=0;iprop<npropS0;iprop++)
	  revert_prop_from_new_contraction(S0[r][iprop],temp_transp);
    
  //free memory
  if(compute_der>=1) nissa_free(temp_der1);
  if(compute_der>=2) nissa_free(temp_der2);
  if(use_new_contraction_layout)
    {
      nissa_free(temp_transp);
      nissa_free(new_loc_2pts);
    }
  else nissa_free(loc_2pts);
  
  contr_2pts_time+=take_time();
  
  close_file(fout);
}

//Calculate and print to file the 3pts
void calculate_all_3pts(int ispec,int ism_lev_so,int ism_lev_se)
{
  char path[1024];
  PROP_TYPE *temp_transp=use_new_contraction_layout?nissa_malloc("temp_trans",loc_vol,PROP_TYPE):NULL;
  complex *loc_3pts;
  double *new_loc_3pts;
  if(!use_new_contraction_layout) loc_3pts=nissa_malloc("contr_3pts",std::max(ncontr_3pts,nch_contr_3pts)*glb_size[0],complex);
  else new_loc_3pts=nissa_malloc("contr_3pts",ncontr_3pts*glb_size[0],double);
  
  //open output file and take time
  sprintf(path,"%s/3pts_sp%d_%02d_%02d",outfolder,ispec,gaussian_niter_so[ism_lev_so],gaussian_niter_se[ism_lev_se]);  
  FILE *fout=open_text_file_for_output(path);
  
  contr_3pts_time-=take_time();
  
  //select r
  int r1=r_spec[ispec];
  int r2=!r1;
  
  //pass to optimal layout for contractions
  if(use_new_contraction_layout)
    {
      for(int ith2=0;ith2<nthetaS1;ith2++)
	for(int im2=0;im2<nmassS1;im2++)
	  prepare_prop_for_new_contraction(S1[ipropS1(ith2,im2)],temp_transp);
      
      for(int ith1=0;ith1<nthetaS0;ith1++)
	for(int im1=0;im1<contr_3pts_up_to_S0_mass;im1++)
	  prepare_prop_for_new_contraction(S0[r1][ipropS0(ith1,im1,0)],temp_transp);
    }
  
  //loop on S1 prop
  for(int ith2=0;ith2<nthetaS1;ith2++)
    for(int im2=0;im2<nmassS1;im2++)
      {
	int ip2=ipropS1(ith2,im2);
	if(nch_contr_3pts>0)
	  {
	    //in case revert from new layout
	    chromo_operator_remove_cSW(Cl,cSW);
	    if(use_new_contraction_layout)
	      revert_prop_from_new_contraction(S1[ip2],temp_transp);
#ifdef POINT_SOURCE_VERSION
	    unsafe_apply_chromo_operator_to_su3spinspin(ch_prop,Cl,S1[ip2]);
#else
	    unsafe_apply_chromo_operator_to_colorspinspin(ch_prop,Cl,S1[ip2]);
#endif
	    chromo_operator_include_cSW(Cl,cSW);
	    //and return back to new layout
	    if(use_new_contraction_layout)
	      {
		prepare_prop_for_new_contraction(S1[ip2],temp_transp);
		prepare_prop_for_new_contraction(ch_prop,temp_transp);
	      }
	  }
	
	for(int ith1=0;ith1<nthetaS0;ith1++)
	    for(int im1=0;im1<contr_3pts_up_to_S0_mass;im1++)
	      {
		int ip1=ipropS0(ith1,im1,0);
		
		//header
		master_fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d,",massS0[im1],thetaS0[ith1],r1,massS1[im2],thetaS1[ith2],r2);
		master_fprintf(fout," smear_source=%d smear_seq=%d\n",gaussian_niter_so[ism_lev_so],gaussian_niter_se[ism_lev_se]);
		
		//compute contractions
		if(use_new_contraction_layout)
		  new_meson_two_points(new_contr_3pts,new_loc_3pts,S0[r1][ip1],S1[ip2],&three_pts_comp);
		else meson_two_points_Wilson_prop(contr_3pts,loc_3pts,op_sour_3pts,S0[r1][ip1],op_sink_3pts,
					     S1[ip2],ncontr_3pts);
		
		ncontr_tot+=ncontr_3pts;
		contr_save_time-=take_time();
		
		//write them
		if(use_new_contraction_layout) three_pts_comp.print_correlations_to_file(fout,new_contr_3pts);

		else print_contractions_to_file(fout,ncontr_3pts,op_sour_3pts,op_sink_3pts,contr_3pts,source_coord[0],"",1.0);
		
		contr_save_time+=take_time();
		
		//if chromo contractions
		/*
		if(nch_contr_3pts>0)
		  {
		    //compute them
		    if(use_new_contraction_layout)
		      new_meson_two_points(ch_contr_3pts,loc_3pts,ch_op_sour_3pts,S0[r1][ip1],
						 ch_op_sink_3pts,ch_prop,nch_contr_3pts);
		    else meson_two_points_Wilson_prop(ch_contr_3pts,loc_3pts,ch_op_sour_3pts,S0[r1][ip1],
						      ch_op_sink_3pts,ch_prop,nch_contr_3pts);
		    
		    //and write them
#ifdef BENCH
		    ncontr_tot+=nch_contr_3pts;
		    contr_save_time-=take_time();
#endif
		    (use_new_contraction_layout?
		     print_optimized_contractions_to_file:print_contractions_to_file)
		      (fout,nch_contr_3pts,ch_op_sour_3pts,ch_op_sink_3pts,ch_contr_3pts,source_coord[0],"CHROMO-",1.0);
#ifdef BENCH
		    contr_save_time+=take_time();
#endif
		  }
		*/
		master_fprintf(fout,"\n");
	    }
    }
  
  
  contr_3pts_time+=take_time();
  
  //convert back
  if(use_new_contraction_layout)
    {
      for(int ith2=0;ith2<nthetaS1;ith2++)
        for(int im2=0;im2<nmassS1;im2++)
	  revert_prop_from_new_contraction(S1[ipropS1(ith2,im2)],temp_transp);
      
      for(int ith1=0;ith1<nthetaS0;ith1++)
	for(int im1=0;im1<contr_3pts_up_to_S0_mass;im1++)
	  revert_prop_from_new_contraction(S0[r1][ipropS0(ith1,im1,0)],temp_transp);
    }
  
  //close and free
  close_file(fout);
  if(use_new_contraction_layout)
    {
      nissa_free(temp_transp);
      nissa_free(new_loc_3pts);
    }
  else nissa_free(loc_3pts);
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
	contract_with_source(contr_2pts,S1[ip2],op_sink_2pts,original_source);
	
	master_fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",
		       massS0[imass_spec[ispec]],thetaS0[ith_spec[ispec]],r_spec[ispec],massS1[im2],thetaS1[ith2],r_spec[ispec]);
	
	master_fprintf(fout,"\n");
	
	for(int icontr=0;icontr<ncontr_2pts;icontr++)
	  if(op_sour_2pts[icontr]==5)
	    master_fprintf(fout," # P5%s\t%+016.16g\t%+016.16g\n",gtag[op_sink_2pts[icontr]],(contr_2pts+icontr*glb_size[0])[source_coord[0]][0],(contr_2pts+icontr*glb_size[0])[source_coord[0]][1]);
	master_fprintf(fout,"\n");
      }
  
  close_file(fout);
}

//check if the time is enough
int check_remaining_time()
{
  int enough_time;
  
  //check remaining time
  double temp_time=take_time()+tot_prog_time;
  double ave_time=temp_time/nanalyzed_conf;
  double left_time=wall_time-temp_time;
  enough_time=left_time>(ave_time*1.1);
  
  master_printf("Remaining time: %lg sec\n",left_time);
  master_printf("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
  if(enough_time) master_printf("Continuing with next conf!\n");
  else master_printf("Not enough time, exiting!\n");
  
  return enough_time;
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //initialize the program
  if(narg<2) crash("Use: %s input_file",arg[0]);
  initialize_semileptonic(arg[1]);
  
  //loop over the configs
  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf and enough_time and !file_exists("stop") and read_conf_parameters(&iconf))
    {
      //smear the conf and generate the source
      setup_conf();
      generate_source();
      
      //loop on smearing of the source
      for(int sm_lev_so=0;sm_lev_so<nsm_lev_so;sm_lev_so++)
	{
	  //compute S0 propagator smearing the config the appropriate number of time the source
	  calculate_all_S0(sm_lev_so);
	  
	  /*
	  //hack to make twisted test
	  double loc[glb_size[0]],glb[glb_size[0]];
	  memset(loc,0,sizeof(double)*glb_size[0]);
	  NISSA_LOC_VOL_LOOP(ivol)
	  {
	    int t=glb_coord_of_loclx[ivol][0];
	    double a=0;
	    for(int id=0;id<4;id++)
	      for(int jd=0;jd<4;jd++)
		for(int ic=0;ic<3;ic++)
#ifdef POINT_SOURCE_VERSION
		  for(int jc=0;jc<3;jc++)
#endif
		    for(int ri=0;ri<2;ri++)
		      {
			double c=S0[0][0][ivol][ic]
#ifdef POINT_SOURCE_VERSION
			  [jc]
#endif
			  [jd][id][ri];
			a+=c*c;
		      }
	    double ph=0;
	    for(int mu=1;mu<4;mu++) ph+=glb_coord_of_loclx[ivol][mu];
	    ph*=M_PI*thetaS0[1]/glb_size[1];
	    loc[t]+=a*cos(ph);
	    //if(ivol==1) printf("ANNA2 %lg %lg\n",cos(ph),sin(ph));
	  }
	  MPI_Allreduce(loc,glb,glb_size[0],MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  for(int t=0;t<glb_size[0];t++)
	    master_printf("ANNA %d %16.16lg\n",t,glb[t]);
	  */
	  
	  //loop on spectator
	  if(ncontr_3pts!=0 or nch_contr_3pts!=0)
	    for(int ispec=0;ispec<nspec;ispec++)
	      {
		//select a timeslice and multiply by gamma5
		generate_sequential_source(ispec);
		
		//loop on smearing of the sequential prop
		for(int sm_lev_se=0;sm_lev_se<nsm_lev_se;sm_lev_se++)
		  {
		    calculate_all_S1(ispec,sm_lev_se);
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
      
      nanalyzed_conf++;
      enough_time=check_remaining_time();
    }
  
  if(iconf==ngauge_conf) master_printf("Finished all the conf!\n");
  
  tot_prog_time+=take_time();
  close_semileptonic();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
