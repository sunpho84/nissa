#include <nissa.hpp>
#include "./new_types.hpp"

int save_base_source;
int load_base_source;
int wall_time;

//conf
gauge_conf_t conf[2];

//source
int noise_type;
char base_source_path[1024];
int nsource;
int *obtain_source_from;
int *obtain_source_applying_sm_op;
in_source_t *source;

//list of masses ant thetas
int nmass_res_group,ntheta_group;
mass_res_group_t *mass_res_group;
theta_group_t *theta_group;

//smearing parameters
ape_smear_pars_t ape_smear_pars;
int ngauss_sm_op;
gauss_smear_pars_t *gauss_sm_op;

//list of propagators
int nprop_group;
prop_group_t *prop_group;
prop_group_command_t *prop_group_command;

//list of operators for 2pts correlation
int ntwo_pts_corr_group;
two_pts_corr_group_t *two_pts_corr_group;

//list of contraction commands
int ncorr_command;
corr_command_t *corr_command;

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
  double kappa;
  read_str_double("Kappa",&kappa);
  conf[0].set_kappa(kappa);

  // 2) Smearing parameters
  
  //Ape smearing parameters
  ape_smear_pars.read();
  //Read the gaussian operators
  read_str_int("NGaussSmOp",&ngauss_sm_op);
  gauss_sm_op=nissa_malloc("GaussSmOp",ngauss_sm_op,gauss_smear_pars_t);
  for(int iop=0;iop<ngauss_sm_op;iop++)
    {
      gauss_sm_op[iop].reset();
      expect_str(combine("GaussSmOp%d",iop).c_str());
      gauss_sm_op[iop].read();
    }
    
  // 3) Read information about the source
  
  //Read the seed and initialize the random generator
  int seed;
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  //Read the noise type
  read_str_int("NoiseType",&noise_type);
  //Read whether we want to load or save the source
  read_str_int("LoadBaseSource",&load_base_source);
  if(!load_base_source) read_str_int("SaveBaseSource",&save_base_source);
  if(load_base_source||save_base_source) read_str_str("BaseSourcePath",base_source_path,1024);
  //Read the number of additional sources
  read_str_int("NAddSource",&nsource);
  nsource++;
  //Allocate the sources and read how to obtain additional sources
  source=nissa_malloc("source*",nsource,in_source_t);
  obtain_source_from=nissa_malloc("obtain",nsource,int);
  obtain_source_applying_sm_op=nissa_malloc("applying",nsource,int);
  for(int iso=1;iso<nsource;iso++)
    {
      expect_str(combine("Source%d",iso).c_str());
      read_str_int("ObtainFromSource",&obtain_source_from[iso]);
      if(obtain_source_from[iso]>=iso||obtain_source_from[iso]<0) CRASH("source %d can be produced from sources in the range [0,%d)",iso,obtain_source_from[iso]);
      read_str_int("ObtainApplyingSmOp",&obtain_source_applying_sm_op[iso]);
      if(obtain_source_applying_sm_op[iso]>=ngauss_sm_op||obtain_source_applying_sm_op[iso]<0)
	crash("selected smearing operator %d not in the defined range [0,%d)",obtain_source_applying_sm_op[iso],ngauss_sm_op);
    }
  
  // 4) Read list of masses and of thetas

  read_str_int("NMassResGroup",&nmass_res_group);
  mass_res_group=nissa_malloc("MassResGroup",nmass_res_group,mass_res_group_t);
  for(int i=0;i<nmass_res_group;i++)
    {
      expect_str(combine("MassResGroup%d",i).c_str());
      mass_res_group[i].read();
    }

  read_str_int("NThetaGroup",&ntheta_group);
  theta_group=nissa_malloc("ThetaGroup",ntheta_group,theta_group_t);
  for(int i=0;i<ntheta_group;i++)
    {
      expect_str(combine("ThetaGroup%d",i).c_str());
      theta_group[i].read();
    }
  
  // 5) Read the list of propagators
  
  read_str_int("NPropGroup",&nprop_group);
  prop_group=nissa_malloc("PropGroup",nprop_group,prop_group_t);
  prop_group_command=nissa_malloc("PropGroupCommand",nprop_group,prop_group_command_t);
  
  for(int iprop_group=0;iprop_group<nprop_group;iprop_group++)
    {
      //read the prop group pars and command
      expect_str(combine("PropGroup%d",iprop_group).c_str());
      prop_group[iprop_group].read_pars(ntheta_group,theta_group,nmass_res_group,mass_res_group);
      prop_group_command[iprop_group].read_command(prop_group[iprop_group],source,prop_group,conf,gauss_sm_op);
    }
  
  // 6) Read the list of correlations
  
  read_str_int("NTwoPtsCorrGroup",&ntwo_pts_corr_group);
  two_pts_corr_group=nissa_malloc("TwoPtsCorrGroup",ntwo_pts_corr_group,two_pts_corr_group_t);
  
  for(int igroup=0;igroup<ntwo_pts_corr_group;igroup++)
    {
      two_pts_corr_group[igroup].reset();
      expect_str(combine("CorrGroup%d",igroup).c_str());
      two_pts_corr_group[igroup].read();
    }
  
  // 7) Read the list of correlations command
  
  read_str_int("NCorrCommand",&ncorr_command);
  corr_command=nissa_malloc("CorrCommand",ncorr_command,corr_command_t);
  
  for(int icomm=0;icomm<ncorr_command;icomm++)
    {
      corr_command[icomm].reset();
      expect_str(combine("CorrCommand%d",icomm).c_str());
      corr_command[icomm].read(ntwo_pts_corr_group,two_pts_corr_group,nprop_group,prop_group);
    }
}

void analysis(char *path,int tsource,char *out_path)
{
  strcpy(base_out_folder,out_path);
  
  //load the conf
  conf[0].read(path);
  
  //generate the ape smeared conf
  conf[1].copy(conf[0]);
  conf[1].ape_smear(ape_smear_pars);
  
  //if needed load the base source
  if(load_base_source) source[0].read(combine("%s/%s",base_out_folder,base_source_path).c_str());
  else
    {
      source[0].fill(nissa_rnd_type_map[RND_Z4],tsource);
      if(save_base_source) source[0].write(combine("%s/%s",base_out_folder,base_source_path).c_str());
    }
  
  //generate all the additional sources
  for(int iso=1;iso<nsource;iso++)
    {
      source[iso]=source[obtain_source_from[iso]];
      source[iso].smear(conf[1],gauss_sm_op[obtain_source_applying_sm_op[iso]]);
    }
  
  //generate all the propagators groups
  for(int iprop_group=0;iprop_group<nprop_group;iprop_group++)
    prop_group_command[iprop_group].exec();
  
  //make all the contractions
  for(int icomm=0;icomm<ncorr_command;icomm++)
    {
      corr_command[icomm].shift=tsource;
      corr_command[icomm].exec();
    }
}

void close_semileptonic()
{
  for(int imass_res_group=0;imass_res_group<nmass_res_group;imass_res_group++) mass_res_group[imass_res_group].destroy();
  nissa_free(mass_res_group);
  for(int iso=0;iso<nsource;iso++) source[iso].destroy();
  nissa_free(source);
  for(int iprop_group=0;iprop_group<nprop_group;iprop_group++) prop_group[iprop_group].destroy();
  nissa_free(prop_group);
  nissa_free(prop_group_command);
  for(int itheta_group=0;itheta_group<ntheta_group;itheta_group++) theta_group[itheta_group].destroy();
  nissa_free(theta_group);
  for(int iop=0;iop<ngauss_sm_op;iop++) gauss_sm_op[iop].destroy();
  nissa_free(gauss_sm_op);
  nissa_free(obtain_source_from);
  nissa_free(obtain_source_applying_sm_op);
  for(int igroup=0;igroup<ntwo_pts_corr_group;igroup++)
    two_pts_corr_group[igroup].destroy();
  nissa_free(two_pts_corr_group);
  conf[0].destroy();
  conf[1].destroy();
  for(int icomm=0;icomm<ncorr_command;icomm++)
    corr_command[icomm].destroy();
  nissa_free(corr_command);

  close_nissa();
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //initialize the program
  if(narg<2) CRASH("Use: %s input_file",arg[0]);
  initialize_semileptonic(arg[1]);
  
  int ngauge_conf;
  read_str_int("NGaugeConf",&ngauge_conf);
  
  for(int iconf=0;iconf<ngauge_conf;iconf++)
    {
      char conf_path[1024];
      char out_folder[1024];
      int tsource;
      
      read_str(conf_path,1024);
      read_int(&tsource);
      read_str(out_folder,1024);
      
      if(!file_exists(combine("%s/running",out_folder).c_str())&&!file_exists(combine("%s/finished",out_folder).c_str()))
	{
	  if(!dir_exists(out_folder)) create_dir(out_folder);
	  file_touch(combine("%s/running",out_folder).c_str());
	  
	  analysis(conf_path,tsource,out_folder);
	  
	  rm(combine("%s/running",out_folder).c_str());
	  file_touch(combine("%s/finished",out_folder).c_str());
	}
      else MASTER_PRINTF("Configuration %s already analized, skipping\n",out_folder);
    }

  close_input();
  
  close_semileptonic();
  
  return 0;
}
