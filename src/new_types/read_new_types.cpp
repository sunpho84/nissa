#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "new_types_definitions.hpp"
#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "io/input.hpp"

namespace nissa
{
  //convert a string into gauge action name
  gauge_action_name_t gauge_action_name_from_str(const char *name)
  {
    gauge_action_name_t ret=UNSPEC_GAUGE_ACTION;
    
    if(strcmp(name,"Wilson")==0) ret=WILSON_GAUGE_ACTION;
    else
      if(strcmp(name,"tlSym")==0) ret=TLSYM_GAUGE_ACTION;
      else crash("unknown gauge action: %s",name);
    
    return ret;
  }

  //read parameters to study topology
  void read_top_meas_pars(top_meas_pars_t &top_meas_pars,bool flag=false)
  {
    if(flag==true) top_meas_pars.flag=true;
    else read_str_int("MeasureTopology",&top_meas_pars.flag);
    if(top_meas_pars.flag)
      {
	char gauge_action_name_str[1024];
	read_str_str("TopPath",top_meas_pars.path,1024);
	read_str_str("TopCoolAction",gauge_action_name_str,1024);
	top_meas_pars.gauge_cooling_action=gauge_action_name_from_str(gauge_action_name_str);
	read_str_int("TopCoolNSteps",&top_meas_pars.cool_nsteps);
	read_str_int("TopCoolOverrelaxing",&top_meas_pars.cool_overrelax_flag);
	if(top_meas_pars.cool_overrelax_flag==1) read_str_double("TopCoolOverrelaxExp",&top_meas_pars.cool_overrelax_exp);
	read_str_int("TopCoolMeasEachNSteps",&top_meas_pars.meas_each);
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
      default: crash("Not implemented yet"); break;
      }
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
  void read_hmc_evol_pars(hmc_evol_pars_t &pars)
  {
    read_str_int("SkipMTestNTraj",&pars.skip_mtest_ntraj);
    read_str_double("HmcTrajLength",&pars.traj_length);
    read_str_int("NmdSteps",&pars.nmd_steps);
    read_str_int("NGaugeSubSteps",&pars.ngauge_substeps);
    read_str_double("MdResidue",&pars.md_residue);
    read_str_double("PfActionResidue",&pars.pf_action_residue);
  }
  
  //read the parameters relevant for pure gauge evolution
  void read_pure_gauge_evol_pars(pure_gauge_evol_pars_t &pars)
  {
    //heat bath parameters
    read_str_int("NHbSweeps",&pars.nhb_sweeps);
    read_str_int("NHbHits",&pars.nhb_hits);
    //overrelax parameters
    read_str_int("NOvSweeps",&pars.nov_sweeps);
    read_str_int("NOvHits",&pars.nov_hits);
  }
  
  //read parameters to stout smear gauge action
  void read_stout_pars(stout_pars_t &stout_pars)
  {
    read_str_int("StoutingNLevel",&stout_pars.nlev);
    if(stout_pars.nlev!=0)
      {
	//isotropic or not?
	int iso;
	read_str_int("IsotropicStouting",&iso);
	
	//only iso implemented
	if(iso)
	  {
	    double rho;
	    read_str_double("StoutRho",&rho);
	    for(int i=0;i<4;i++)
	      for(int j=0;j<4;j++)
		stout_pars.rho[i][j]=rho;
	  }
	else crash("Anisotropic stouting not yet implemented");
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
    else
      for(int i=0;i<3;i++)
	pars.E[i]=pars.B[i]=0;
  }
  
  //read parameters to measure chiral condensate
  void read_chiral_cond_pars(chiral_cond_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureChiralCond",&pars.flag);
    if(pars.flag)
      {
	read_str_str("ChiralCondPath",pars.path,1024);
	read_str_double("ChiralCondInvResidue",&pars.residue);
	read_str_int("ChiralCondNHits",&pars.nhits);
      }
  }
  
  //read parameters to measure magnetization
  void read_magnetization_pars(magnetization_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureMagnetization",&pars.flag);
    if(pars.flag)
      {
	read_str_str("MagnetizationPath",pars.path,1024);
	read_str_double("MagnetizationInvResidue",&pars.residue);
	read_str_int("MagnetizationNHits",&pars.nhits);
      }
  }
  
  //read parameters to measure pseudoscalar correlators
  void read_pseudo_corr_pars(pseudo_corr_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasurePseudoCorr",&pars.flag);
    if(pars.flag)
      {
	read_str_str("PseudoCorrPath",pars.path,1024);
	read_str_double("PseudoCorrInvResidue",&pars.residue);
	read_str_int("PseudoCorrNHits",&pars.nhits);
      }
  }
  
  //read parameters to measure all rectangles
  void read_all_rect_meas_pars(all_rect_meas_pars_t &pars,bool flag=false)
  {
    if(flag==true) pars.flag=true;
    else read_str_int("MeasureAllRect",&pars.flag);
    if(pars.flag)
      {
	read_str_str("AllRectPath",pars.path,1024);
	
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
	
	read_str_double("APESpatAlpha",&pars.ape_spat_alpha);
	read_list_of_ints("APESpatNlevels",&pars.nape_spat_levls,&pars.nape_spat_iters);
	
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
    read_str_int("NDiffFlavs",&(theory_pars.nflavs));
    theory_pars.quark_content=nissa_malloc("quark_content",theory_pars.nflavs,quark_content_t);
    
    //each flav parameters
    for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
      read_quark_content(theory_pars.quark_content[iflav]);
    
    //additional parameters to read only if fermions are defined
    if(theory_pars.nflavs!=0)
      {
	//stouting parameters
	read_stout_pars(theory_pars.stout_pars);
	
	//electric and magnetic field
	read_em_field_pars(theory_pars.em_field_pars);
	
	//info on pseudoscalar meson correlators, condensate and magnetization measure
	read_pseudo_corr_pars(theory_pars.pseudo_corr_pars);
	read_chiral_cond_pars(theory_pars.chiral_cond_pars);
	read_magnetization_pars(theory_pars.magnetization_pars);
      }
  }
}
