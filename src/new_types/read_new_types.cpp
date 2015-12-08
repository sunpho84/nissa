#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "io/input.hpp"
#include "new_types_definitions.hpp"
#include "operations/stag/nucleon.hpp"

namespace nissa
{
  //convert a string into gauge action name
  gauge_action_name_t gauge_action_name_from_str(const char *name)
  {
    //database
    const int nact_known=3;
    gauge_action_name_t act_known[nact_known]={WILSON_GAUGE_ACTION,TLSYM_GAUGE_ACTION,IWASAKI_GAUGE_ACTION};
    const char name_known[nact_known][20]={"Wilson","tlSym","Iwasaki"};
    
    //search
    int iact=0;
    while(iact<nact_known && strcasecmp(name,name_known[iact])!=0) iact++;
    
    //check
    if(iact==nact_known) crash("unknown gauge action: %s",name);
    
    return act_known[iact];
  }
  
  //convert a string into smoothing method
  smooth_pars_t::method_t smooth_method_name_from_str(const char *name)
  {
    //database
    const int nmet_known=4;
    smooth_pars_t::method_t met_known[nmet_known]={smooth_pars_t::COOLING,smooth_pars_t::STOUTING,smooth_pars_t::ADAPTATIVE_STOUTING,smooth_pars_t::WFLOWING};
    const char name_known[nmet_known][20]={"Cooling","Stouting","AdaptativeStouting","Wflowing"};
    
    //search
    int imet=0;
    while(imet<nmet_known && strcasecmp(name,name_known[imet])!=0) imet++;
    
    //check
    if(imet==nmet_known) crash("unknown smoothing method: %s",name);
    
    return met_known[imet];
  }
  
  //read parameters to stout smear gauge action
  void read_stout_pars(stout_pars_t &stout_pars)
  {
    read_str_int("StoutingNLevel",&stout_pars.nlevels);
    if(stout_pars.nlevels!=0)
      {
	//isotropic or not?
	int iso;
	read_str_int("IsotropicStouting",&iso);
	
	//only iso implemented
	if(iso) read_str_double("StoutRho",&stout_pars.rho);
	else crash("Anisotropic stouting not yet implemented");
      }
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
  void read_Wflow_pars(Wflow_pars_t &pars)
  {
    read_str_double("FlowTime",&pars.T);
    read_str_double("InteStep",&pars.dt);
  }
  
  //read parameters to adaptative stout
  void read_adaptative_stout_pars(adaptative_stout_pars_t &pars)
  {
    read_str_int("Nlevels",&pars.nlevels);
    for(int ilev=0;ilev<pars.nlevels;ilev++)
      {
	double r;
	read_double(&r);
	pars.rho.push_back(r);
      }
  }
 
  //read parameters to smooth
  void read_smooth_pars(smooth_pars_t &smooth_pars,bool flag=false)
  {
    if(flag==true) smooth_pars.flag=true;
    else read_str_int("Smoothing",&smooth_pars.flag);
    if(smooth_pars.flag)
      {
	char smooth_method_name_str[1024];
	read_str_str("SmoothMethod",smooth_method_name_str,1024);
	smooth_pars.method=smooth_method_name_from_str(smooth_method_name_str);
	switch(smooth_pars.method)
	  {
	  case smooth_pars_t::COOLING: read_cool_pars(smooth_pars.cool_pars);break;
	  case smooth_pars_t::STOUTING: read_stout_pars(smooth_pars.stout_pars);break;
	  case smooth_pars_t::ADAPTATIVE_STOUTING: read_adaptative_stout_pars(smooth_pars.adaptative_stout_pars);break;
	  case smooth_pars_t::WFLOWING: read_Wflow_pars(smooth_pars.Wflow_pars);break;
	  case smooth_pars_t::UNSPEC_SMOOTH_METHOD: crash("should not arrive here");break;
	  }
	read_str_double("MeasEach",&smooth_pars.meas_each);
	if((smooth_pars.method==smooth_pars_t::COOLING||smooth_pars.method==smooth_pars_t::STOUTING)&&fabs(smooth_pars.meas_each-int(smooth_pars.meas_each))>=1.e-14)
	  crash("MeasEach must be integer if Cooling or Stouting method selected");
      }
  }
  
  //read and return path
  std::string read_path()
  {
    char temp[1024];
    read_str_str("Path",temp,1024);
    return temp;
  }
  
  //read parameters to study topology
  void read_top_meas_pars(top_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureTopology",&pars.flag);
    if(pars.flag)
      {
	pars.path=read_path();
	read_smooth_pars(pars.smooth_pars,true);
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
  
  //the parameters relevant for hmc evolution
  void read_hmc_evol_pars(hmc_evol_pars_t &pars,theory_pars_t &th)
  {
    read_str_int("SkipMTestNTraj",&pars.skip_mtest_ntraj);
    read_str_double("HmcTrajLength",&pars.traj_length);
    read_str_int("NmdSteps",&pars.nmd_steps);
    read_str_int("NGaugeSubSteps",&pars.ngauge_substeps);
    read_str_double("MdResidue",&pars.md_residue);
    read_str_double("PfActionResidue",&pars.pf_action_residue);
    pars.npseudo_fs=new int[th.nflavs()];
    expect_str("NPseudoFermions");
    for(int iflav=0;iflav<th.nflavs();iflav++) read_int(&pars.npseudo_fs[iflav]);
    pars.rat_appr=new rat_approx_t[3*th.nflavs()];
    for(int i=0;i<th.nflavs()*3;i++) pars.rat_appr[i].reset();
  }
  
  //read the parameters relevant for pure gauge evolution
  void read_pure_gauge_evol_pars(pure_gauge_evol_pars_t &pars)
  {
    //use or not hybrid Monte Carlo
    read_str_int("UseHMC",&pars.use_hmc);
    if(pars.use_hmc)
      {
	read_str_double("HmcTrajLength",&pars.traj_length);
	read_str_int("NmdSteps",&pars.nmd_steps);
	read_str_int("UseFacc",&pars.use_Facc);
	if(pars.use_Facc)
	  {
	    read_str_double("Kappa",&pars.kappa);
	    read_str_double("Residue",&pars.residue);
	  }
      }
    else
      {
	//heat bath parameters
	read_str_int("NHbSweeps",&pars.nhb_sweeps);
	read_str_int("NHbHits",&pars.nhb_hits);
	//overrelax parameters
	read_str_int("NOvSweeps",&pars.nov_sweeps);
	read_str_int("NOvHits",&pars.nov_hits);
      }
  }
  
  //read parameters to ape smear gauge action
  void read_ape_pars(ape_pars_t &ape_pars)
  {
    read_str_int("ApeNLevel",&ape_pars.nlev);
    read_str_double("ApeAlpha",&ape_pars.alpha);
  }
  
  //read parameters of the background em field
  void read_em_field_pars(em_field_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("PutBkgrdEMField",&pars.flag);
    if(pars.flag)
      {
	read_str_double("Ex",&(pars.E[0]));
	read_str_double("Ey",&(pars.E[1]));
	read_str_double("Ez",&(pars.E[2]));
	read_str_double("Bx",&(pars.B[0]));
	read_str_double("By",&(pars.B[1]));
	read_str_double("Bz",&(pars.B[2]));
      }
  }
  
  //read parameters to measure the fermionic gran mix
  void read_fermionic_putpourri_meas_pars(fermionic_putpourri_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureFermionicPutpourri",&pars.flag);
    if(pars.flag)
      {
	pars.path=read_path();
	read_str_double("InvResidue",&pars.residue);
	read_str_int("ComputeSusceptivities",&pars.compute_susc);
	read_str_int("NCopies",&pars.ncopies);
	read_str_int("NHits",&pars.nhits);
      }
  }
  
  //read parameters to measure the quark density and its derivatives
  void read_quark_rendens_meas_pars(quark_rendens_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureQuarkRendens",&pars.flag);
    if(pars.flag)
      {
	pars.path=read_path();
        read_str_double("InvResidue",&pars.residue);
        read_str_int("NCopies",&pars.ncopies);
        read_str_int("NHits",&pars.nhits);
      }
  }
  
  //read parameters to measure the spin polarization
  void read_spinpol_meas_pars(spinpol_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureSpinpol",&pars.flag);
    if(pars.flag)
      {
	pars.path=read_path();
        read_str_double("InvResidue",&pars.residue);
        read_str_int("Dir",&pars.dir);
        read_str_int("NHits",&pars.nhits);
	read_str_int("UseFermConfForGluons",&pars.use_ferm_conf_for_gluons);
	read_smooth_pars(pars.smooth_pars);
      }
  }
  
  //read parameters to measure magnetization
  void read_magnetization_meas_pars(magnetization_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureMagnetization",&pars.flag);
    if(pars.flag)
      {
	pars.path=read_path();
	read_str_double("InvResidue",&pars.residue);
	read_str_int("NCopies",&pars.ncopies);
	read_str_int("NHits",&pars.nhits);
      }
  }
  
  //read parameters to measure pseudoscalar correlators
  void read_pseudo_corr_meas_pars(pseudo_corr_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasurePseudoCorr",&pars.flag);
    if(pars.flag)
      {
	pars.path=read_path();
	read_str_double("InvResidue",&pars.residue);
	read_str_int("NHits",&pars.nhits);
      }
  }
  
  //read parameters to measure nucleon correlators
  void read_nucleon_corr_meas_pars(nucleon_corr_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureNucleonCorr",&pars.flag);
    if(pars.flag)
      {
	pars.path=read_path();
	read_str_double("InvResidue",&pars.residue);
	read_str_int("NHits",&pars.nhits);
      }
  }
  
  //read parameters to smear a conf in time when computing gauge observables
  void read_gauge_obs_temp_smear_pars(gauge_obs_temp_smear_pars_t &pars)
  {
    read_str_int("UseHYPorAPE",&pars.use_hyp_or_ape_temp);
    if(pars.use_hyp_or_ape_temp==0)
      {
	read_str_double("HYPTempAlpha0",&pars.hyp_temp_alpha0);
	read_str_double("HYPTempAlpha1",&pars.hyp_temp_alpha1);
	read_str_double("HYPTempAlpha2",&pars.hyp_temp_alpha2);
      }
    else
      {
	read_str_double("APETempAlpha",&pars.ape_temp_alpha);
	read_str_int("APETempNiters",&pars.nape_temp_iters);
      }
  }
  
  //read parameters to smear a conf in space and time when computing gauge observables
  void read_gauge_obs_temp_spat_smear_pars(gauge_obs_temp_spat_smear_pars_t &pars)
  {
    read_gauge_obs_temp_smear_pars(pars.gauge_temp_smear_pars);
    read_str_double("APESpatAlpha",&pars.ape_spat_alpha);
    read_list_of_ints("APESpatNlevels",&pars.nape_spat_levls,&pars.nape_spat_iters);
  }
  
  //read parameters to measure all rectangles
  void read_all_rect_meas_pars(all_rect_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureAllRect",&pars.flag);
    if(pars.flag)
      {
	pars.path=read_path();
	read_gauge_obs_temp_spat_smear_pars(pars.smear_pars);
	read_str_int("Tint",&pars.Tmin);
	read_int(&pars.Tmax);
	if(pars.Tmin<1||pars.Tmin>glb_size[0]) crash("Tint must start between [1,%d]",glb_size[0]);
	if(pars.Tmax<1||pars.Tmax>glb_size[0]) crash("Tint must end between [1,%d]",glb_size[0]);
	read_str_int("Dint",&pars.Dmin);
	read_int(&pars.Dmax);
	if(pars.Dmin<1||pars.Dmin>glb_size[1]) crash("Dint must start between [1,%d]",glb_size[1]);
	if(pars.Dmax<1||pars.Dmax>glb_size[1]) crash("Dint must end between [1,%d]",glb_size[1]);
      }
  }

  //read parameters to measure flux tube
  void read_watusso_meas_pars(watusso_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureWatusso",&pars.flag);
    if(pars.flag)
      {
	pars.path=read_path();
	read_gauge_obs_temp_spat_smear_pars(pars.smear_pars);
	read_str_int("SizeMin",&pars.size_min);
	read_str_int("SizeStep",&pars.size_step);
	read_str_int("SizeMax",&pars.size_max);
	read_str_int("dmax",&pars.dmax);
	for(int mu=0;mu<4;mu++)
	  {
	    if(pars.size_min>glb_size[mu]) crash("SizeMin must be between [1,%d]",glb_size[mu]);
	    if(pars.size_max>glb_size[mu]) crash("SizeMax must be between [1,%d]",glb_size[mu]);
	    if(pars.dmax>glb_size[mu]) crash("DMax must be between [1,%d]",glb_size[mu]);
	  }
      }
  }

  //read the parameters to compute polyakov correlators
  void read_poly_corr_meas_pars(poly_corr_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasurePolyCorrs",&pars.flag);
    if(pars.flag)
      {
	pars.path=read_path();
        read_gauge_obs_temp_smear_pars(pars.gauge_smear_pars);
	read_str_int("Dir",&pars.dir);
      }
  }
  
  //read the parameters to measure gauge observables
  void read_gauge_obs_meas_pars(gauge_obs_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureGaugeObs",&pars.flag);
    if(pars.flag) pars.path=read_path();
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
    read_str_int("NDiffFlavs",&(theory_pars.nstag_flavs));
    theory_pars.quark_content=nissa_malloc("quark_content",theory_pars.nflavs(),quark_content_t);
    
    //each flav parameters
    for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
      read_quark_content(theory_pars.quark_content[iflav]);
    
    //additional parameters to read only if fermions are defined
    if(theory_pars.nflavs()!=0)
      {
	//stouting parameters
	read_stout_pars(theory_pars.stout_pars);
	
	//electric and magnetic field
	read_em_field_pars(theory_pars.em_field_pars);
      }
  }
}
