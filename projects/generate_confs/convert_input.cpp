#include "nissa.hpp"

#include "driver.hpp"

using namespace nissa;

int parser_lex_destroy  (driver_t* yyscanner);

  //convert a string into smoothing method
  smooth_pars_t::method_t smooth_method_name_from_str(const char *name)
  {
    //database
    const int nmet_known=3;
    smooth_pars_t::method_t met_known[nmet_known]={smooth_pars_t::COOLING,smooth_pars_t::STOUT,smooth_pars_t::WFLOW};
    const char name_known[nmet_known][20]={"Cooling","Stouting","Wflowing"};
    
    //search
    int imet=0;
    while(imet<nmet_known && strcasecmp(name,name_known[imet])!=0) imet++;
    
    //check
    if(imet==nmet_known) crash("unknown smoothing method: %s",name);
    
    return met_known[imet];
  }

//read parameters to cool
void read_cool_pars(cool_pars_t &cool_pars)
{
  char gauge_action_name_str[1024];
  read_str_str("CoolAction",gauge_action_name_str,1024);
  cool_pars.gauge_action=gauge_action_name_from_str(gauge_action_name_str);
  read_str_int("CoolNSteps",&cool_pars.nsteps);
  read_str_int("CoolOverrelaxing",&cool_pars.overrelax_flag);
  if(cool_pars.overrelax_flag==1) read_str_double("CoolOverrelaxExp",&cool_pars.overrelax_exp);
}

//read parameters to flow
void legacy_read_Wflow_pars(Wflow_pars_t &pars)
{
  double T;
  read_str_double("T",&T);
  read_str_double("InteStep",&pars.dt);
  pars.nflows=T/pars.dt+0.5;
}

//read parameters to smooth
void read_smooth_pars(smooth_pars_t &smooth_pars,int flag=false)
{
  if(!flag==true) read_str_int("Smoothing",&flag);
  if(flag)
    {
      char smooth_method_name_str[1024];
      read_str_str("SmoothMethod",smooth_method_name_str,1024);
      smooth_pars.method=smooth_method_name_from_str(smooth_method_name_str);
      switch(smooth_pars.method)
	{
	case smooth_pars_t::COOLING: read_cool_pars(smooth_pars.cool);break;
	case smooth_pars_t::STOUT: read_stout_pars(smooth_pars.stout);break;
	case smooth_pars_t::WFLOW: legacy_read_Wflow_pars(smooth_pars.Wflow);break;
	default: crash("should not arrive here");break;
	}
      double each;
      read_str_double("MeasEach",&each);
      if(smooth_pars.method==smooth_pars_t::WFLOW) smooth_pars.meas_each_nsmooth=each/smooth_pars.Wflow.dt+0.5;
      else smooth_pars.meas_each_nsmooth=each;
    }
}

//topological potential
void read_topotential_pars(topotential_pars_t &pars,int flag=0)
{
  if(flag!=0) pars.flag=flag;
  else read_str_int("TopoPotential",&pars.flag);
  switch(pars.flag)
    {
    case 0: break;
    case 1: read_str_double("Potential",&pars.theta); break;
    case 2:
      pars.read_pars();
      break;
    default: crash("Not implemented yet"); break;
    }
  if(pars.flag) read_stout_pars(pars.stout_pars);
}

//degeneracy, mass, chpot and charge
void read_quark_content(quark_content_t &quark_content,bool flag=false)
{
  read_str_int("Degeneracy",&(quark_content.deg));
  read_str_double("Mass",&(quark_content.mass));
  read_str_double("RePotCh",&(quark_content.re_pot));
  read_str_double("ImPotCh",&(quark_content.im_pot));
  read_str_double("ElecCharge",&(quark_content.charge));
}

//read parameters of the background em field
void read_em_field_pars(em_field_pars_t &pars)
{
  int flag;
  read_str_int("PutBkgrdEMField",&flag);
  if(flag)
    {
      read_str_double("Ex",&(pars.E[0]));
      read_str_double("Ey",&(pars.E[1]));
      read_str_double("Ez",&(pars.E[2]));
      read_str_double("Bx",&(pars.B[0]));
      read_str_double("By",&(pars.B[1]));
      read_str_double("Bz",&(pars.B[2]));
    }
}

//read and return path
std::string read_path()
{
  char temp[1024];
  read_str_str("Path",temp,1024);
  return temp;
}

//read the theory_pars parameters of the theory
void read_theory_pars(theory_pars_t &theory_pars)
{
  //kind of action
  char gauge_action_name[1024];
  read_str_str("GaugeAction",gauge_action_name,1024);
  theory_pars.gauge_action_name=gauge_action_name_from_str(gauge_action_name);
  
  //beta for gauge action
  read_str_double("Beta",&theory_pars.beta);
  
  //topological potential
  read_topotential_pars(theory_pars.topotential_pars);
  
  //number of undegenerate flavs
  int nflavs;
  read_str_int("NDiffFlavs",&nflavs);
  theory_pars.quarks.resize(nflavs);
  
  //each flav parameters
  for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
    read_quark_content(theory_pars.quarks[iflav]);
  
  //additional parameters to read only if fermions are defined
  if(theory_pars.nflavs()!=0)
    {
      //stouting parameters
      read_stout_pars(theory_pars.stout_pars);
      
      //electric and magnetic field
      read_em_field_pars(theory_pars.em_field_pars);
    }
}

//read parameters to measure nucleon correlators
void read_nucleon_corr_meas_pars(std::vector<nucleon_corr_meas_pars_t> &pars,int itheory)
{
  int flag;
  read_str_int("MeasureNucleonCorr",&flag);
  if(flag)
    {
      pars.push_back(nucleon_corr_meas_pars_t());
      pars.back().itheory=itheory;
      pars.back().path=read_path();
      read_str_double("InvResidue",&pars.back().residue);
      read_str_int("NHits",&pars.back().nhits);
    }
}

//read how to measure staggered mesons
void read_stag_meson_corr_meas_pars(std::vector<meson_corr_meas_pars_t> &pars,int itheory)
{
  int flag;
  read_str_int("MeasureStagMesonCorr",&flag);
  if(flag)
    {
      pars.push_back(meson_corr_meas_pars_t());
      pars.back().itheory=itheory;
      pars.back().path=read_path();
      int nmesons;
      read_str_int("NMesons",&nmesons);
      for(int imeson=0;imeson<nmesons;imeson++)
	{
	  int spin,taste;
	  read_int(&spin);
	  read_int(&taste);
	  pars.back().mesons.push_back(std::make_pair(spin,taste));
	}
      read_str_double("InvResidue",&pars.back().residue);
      read_str_int("NCopies",&pars.back().ncopies);
      read_str_int("NHits",&pars.back().nhits);
    }
}

//read parameters to measure the fermionic gran mix
void read_fermionic_putpourri_meas_pars(std::vector<fermionic_putpourri_meas_pars_t> &pars,int itheory)
{
  int flag;
  read_str_int("MeasureFermionicPutpourri",&flag);
  if(flag)
    {
      pars.push_back(fermionic_putpourri_meas_pars_t());
      pars.back().itheory=itheory;
      pars.back().path=read_path();
      read_str_double("InvResidue",&pars.back().residue);
      read_str_int("ComputeSusceptivities",&pars.back().compute_susc);
      read_str_int("NCopies",&pars.back().ncopies);
      read_str_int("NHits",&pars.back().nhits);
    }
}

//read parameters to measure the quark density and its derivatives
void read_quark_rendens_meas_pars(std::vector<quark_rendens_meas_pars_t> &pars,int itheory)
{
  int flag;
  read_str_int("MeasureQuarkRendens",&flag);
  if(flag)
    {
      pars.push_back(quark_rendens_meas_pars_t());
      pars.back().itheory=itheory;
      pars.back().path=read_path();
      read_str_int("After",&pars.back().after);
      read_str_int("MaxOrder",&pars.back().max_order);
      read_str_double("InvResidue",&pars.back().residue);
      read_str_int("NCopies",&pars.back().ncopies);
      read_str_int("NHits",&pars.back().nhits);
    }
}

//read parameters to measure magnetization
void read_magnetization_meas_pars(std::vector<magnetization_meas_pars_t> &pars,int itheory)
{
  int flag;
  read_str_int("MeasureMagnetization",&flag);
  if(flag)
    {
      pars.push_back(magnetization_meas_pars_t());
      pars.back().itheory=itheory;
      pars.back().path=read_path();
      read_str_double("InvResidue",&pars.back().residue);
      read_str_int("NCopies",&pars.back().ncopies);
      read_str_int("NHits",&pars.back().nhits);
    }
}

//read parameters to measure the spin polarization
void read_spinpol_meas_pars(std::vector<spinpol_meas_pars_t> &pars,int itheory)
{
  int flag;
  read_str_int("MeasureSpinpol",&flag);
  if(flag)
    {
      pars.push_back(spinpol_meas_pars_t());
      pars.back().itheory=itheory;
      pars.back().path=read_path();
      read_str_double("InvResidue",&pars.back().residue);
      //read_str_int("Dir",&pars.back().dir);
      read_str_int("NHits",&pars.back().nhits);
      read_str_int("UseFermConfForGluons",&pars.back().use_ferm_conf_for_gluons);
      read_smooth_pars(pars.back().smooth_pars);
    }
}

//read parameters to smear a conf in time when computing gauge observables
void read_gauge_obs_temp_smear_pars(smooth_pars_t &pars)
{
  int use_hyp_or_ape_temp;
  read_str_int("UseHYPorAPE",&use_hyp_or_ape_temp);
  if(use_hyp_or_ape_temp==0)
    {
      pars.method=smooth_pars_t::HYP;
      read_str_double("HYPTempAlpha0",&pars.hyp.alpha0);
      read_str_double("HYPTempAlpha1",&pars.hyp.alpha1);
      read_str_double("HYPTempAlpha2",&pars.hyp.alpha2);
    }
  else
    {
      pars.method=smooth_pars_t::APE;
      read_str_double("APETempAlpha",&pars.ape.alpha);
      read_str_int("APETempNiters",&pars.ape.nlevels);
    }
}

//read the parameters to measure gauge observables
void read_gauge_obs_meas_pars(std::vector<gauge_obs_meas_pars_t> &pars,int flag=false)
{
  if(!flag) read_str_int("MeasureGaugeObs",&flag);
  if(flag)
    {
      pars.push_back(gauge_obs_meas_pars_t());
      pars.back().path=read_path();
    }
}

//read the parameters to compute polyakov correlators
void read_poly_corr_meas_pars(std::vector<poly_corr_meas_pars_t> &pars,int flag=false)
{
  if(!flag) read_str_int("MeasurePolyCorrs",&flag);
  if(flag)
    {
      pars.push_back(poly_corr_meas_pars_t());
      pars.back().path=read_path();
      read_gauge_obs_temp_smear_pars(pars.back().smear_pars);
      read_str_int("Dir",&pars.back().dir);
    }
}

//read parameters to study topology
void read_top_meas_pars(top_meas_pars_t &pars,int flag=false)
{
  if(!flag) read_str_int("MeasureTopology",&flag);
  if(flag)
    {
      pars.path=read_path();
      read_smooth_pars(pars.smooth_pars,true);
    }
}

//read parameters to measure all rectangles
void read_all_rect_meas_pars(std::vector<all_rects_meas_pars_t> &pars,int flag=false)
{
  if(!flag) read_str_int("MeasureAllRect",&flag);
  if(flag)
    {
      pars.push_back(all_rects_meas_pars_t());
      pars.back().path=read_path();
      read_gauge_obs_temp_smear_pars(pars.back().temp_smear_pars);
      read_str_int("Tint",&pars.back().Tmin);
      read_int(&pars.back().Tmax);
      read_str_int("Dint",&pars.back().Dmin);
      read_int(&pars.back().Dmax);
    }
}

//read parameters to measure flux tube
void read_watusso_meas_pars(std::vector<watusso_meas_pars_t> &pars,int flag=false)
{
  if(!flag) read_str_int("MeasureWatusso",&flag);
  if(flag)
    {
      pars.push_back(watusso_meas_pars_t());
      pars.back().path=read_path();
      read_gauge_obs_temp_smear_pars(pars.back().temp_smear_pars);
      pars.back().spat_smear_pars.method=smooth_pars_t::APE;
      read_ape_pars(pars.back().spat_smear_pars.ape);
      int size_min,size_step,size_max;
      read_str_int("SizeMin",&size_min);
      read_str_int("SizeStep",&size_step);
      read_str_int("SizeMax",&size_max);
      for(int size=size_min;size<=size_max;size+=size_step) pars.back().sizes.push_back(size);
      read_str_int("dmax",&pars.back().dmax);
    }
}

//the parameters relevant for hmc evolution
void read_hmc_evol_pars(hmc_evol_pars_t &pars,theory_pars_t &th)
{
  read_str_int("SkipMTestNTraj",&pars.skip_mtest_ntraj);
  read_str_double("HmcTrajLength",&pars.traj_length);
  read_str_int("NmdSteps",&pars.nmd_steps);
  read_str_int("NGaugeSubSteps",&pars.ngauge_substeps);
  read_str_double("MdResidue",&pars.md_residue);
  read_str_double("PfActionResidue",&pars.pf_action_residue);
  pars.npseudo_fs.resize(th.nflavs());
  expect_str("NPseudoFermions");
  for(int iflav=0;iflav<th.nflavs();iflav++) read_int(&pars.npseudo_fs[iflav]);
}

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<3) crash("Use: %s input_file output_file",arg[0]);
  
  //open input file
  open_input(arg[1]);
  
  driver_t *driver=new driver_t(input_global);
  
  //init the grid
  int L;
  read_str_int("L",&L);
  if(L>0)
    {
      read_str_int("T",&driver->T);
      driver->LX=driver->LY=driver->LZ=L;
    }
  else
    {
      read_str_int("LT",&driver->T);
      read_str_int("LX",&driver->LX);
      read_str_int("LY",&driver->LY);
      read_str_int("LZ",&driver->LZ);
   }
  
  //read number of additional theory to read
  int nvalence_theories;
  read_str_int("NValenceTheories",&nvalence_theories);
  int ntheories=nvalence_theories+1;
  driver->theories.resize(ntheories);
  std::vector<spinpol_meas_pars_t> spinpol_meas(ntheories);
  
  //read physical theory: theory 0 is the sea (simulated one)
  for(int itheory=0;itheory<ntheories;itheory++)
    {
      if(itheory==0) master_printf("Reading info on sea theory\n");
      else           master_printf("Reading info on additional (valence) theory %d/%d\n",itheory,nvalence_theories);
      read_theory_pars(driver->theories[itheory]);
      read_nucleon_corr_meas_pars(driver->nucleon_corr_meas,itheory);
      read_stag_meson_corr_meas_pars(driver->meson_corr_meas,itheory);
      read_fermionic_putpourri_meas_pars(driver->fermionic_putpourri_meas,itheory);
      read_quark_rendens_meas_pars(driver->quark_rendens_meas,itheory);
      read_spinpol_meas_pars(spinpol_meas,itheory);
      read_magnetization_meas_pars(driver->magnetization_meas,itheory);
    }
  
  //read if we want to measure gauge obs
  read_gauge_obs_meas_pars(driver->plaq_pol_meas);
  
  //read if we want to measure polyakov correlators
  read_poly_corr_meas_pars(driver->luppoli_meas);
  
  //read if we want to measure topological charge
  int ntop_meas;
  read_str_int("NTopMeas",&ntop_meas);
  driver->top_meas.resize(ntop_meas);
  for(int itop=0;itop<ntop_meas;itop++) read_top_meas_pars(driver->top_meas[itop]);
  
  //read if we want to measure all rectangles
  read_all_rect_meas_pars(driver->all_rects_meas);
  
  //read if we want to measure flux tube
  driver->watusso_meas.resize(1);
  read_watusso_meas_pars(driver->watusso_meas);
  
  //read the number of trajectory to evolve and the wall_time
  read_str_int("NTrajTot",&driver->evol_pars.ntraj_tot);
  read_str_int("WallTime",&driver->walltime);
  
  //read the seed
  read_str_int("Seed",&driver->seed);
  
  //if we want to produce something, let's do it, otherwise load the list of configurations to analyze
  if(driver->evol_pars.ntraj_tot>0)
    {
      driver->run_mode=driver_t::EVOLUTION_MODE;
      
      //load evolution info depending if is a quenched simulation or unquenched
      if(driver->theories[0].nflavs()!=0||driver->theories[0].topotential_pars.flag!=0)
        read_hmc_evol_pars(driver->evol_pars,driver->theories[0]);
      
      //read in and out conf path
      char temp[1024];
      read_str_str("ConfPath",temp,1024);
      driver->conf_pars.path=temp;
      read_str_str("StoreConfPath",temp,1024);
      driver->conf_pars.store_path=temp;
      driver->conf_pars.store_path+=".%04d";
      read_str_int("StoreConfEach",&driver->conf_pars.store_each);
      read_str_int("StoreRunningTempConf",&driver->conf_pars.store_running);
      
      //read if configuration must be generated cold or hot
      char start_conf_cond_str[1024];
      read_str_str("StartConfCond",start_conf_cond_str,1024);
      if(strcasecmp(start_conf_cond_str,"HOT")==0) driver->conf_pars.start_cond=HOT_START_COND;
      if(strcasecmp(start_conf_cond_str,"COLD")==0) driver->conf_pars.start_cond=COLD_START_COND;
    }
  else
    {
      driver->run_mode=driver_t::ANALYSIS_MODE;
      
      //load the number of configurations to analyze
      int nconf_to_analyze;
      read_str_int("NConfToAnalyze",&nconf_to_analyze);
      for(int iconf_to=0;iconf_to<nconf_to_analyze;iconf_to++)
        {
          char temp_path[1024];
          read_str(temp_path,1024);
	  driver->an_conf_list.push_back(temp_path);
        }
    }
  close_input();
  
  FILE *fout=open_file(arg[2],"w");
  driver->master_fprintf(fout);
  close_file(fout);
  
  delete driver;
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}

