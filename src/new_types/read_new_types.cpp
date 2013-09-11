#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <string.h>

#include "new_types_definitions.h"
#include "base/debug.h"
#include "base/vectors.h"
#include "io/input.h"

//read parameters to study topology
void read_top_meas_pars(top_meas_pars_t &top_meas_pars)
{
  read_str_int("MeasureTopology",&top_meas_pars.flag);
  if(top_meas_pars.flag)
    {
      read_str_str("TopPath",top_meas_pars.path,1024);
      read_str_int("TopCoolNSteps",&top_meas_pars.cool_nsteps);
      read_str_int("TopCoolOverrelaxing",&top_meas_pars.cool_overrelax_flag);
      if(top_meas_pars.cool_overrelax_flag==1) read_str_double("TopCoolOverrelaxExp",&top_meas_pars.cool_overrelax_exp);
      read_str_int("TopCoolMeasEachNSteps",&top_meas_pars.meas_each);
    }
}

//read degeneracy, mass, chpot and charge
void read_quark_content(quark_content_t &quark_content)
{
  read_str_int("Degeneracy",&(quark_content.deg));
  read_str_double("Mass",&(quark_content.mass));
  read_str_double("RePotCh",&(quark_content.re_pot));
  read_str_double("ImPotCh",&(quark_content.im_pot));
  read_str_double("ElecCharge",&(quark_content.charge));
}

//read the parameters relevant for hmc evolution
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
void read_em_field_pars(em_field_pars_t &em_field_pars)
{
  read_str_int("PutBkgrdEMField",&em_field_pars.flag);
  if(em_field_pars.flag)
    {
      read_str_double("Ex",&(em_field_pars.E[0]));
      read_str_double("Ey",&(em_field_pars.E[1]));
      read_str_double("Ez",&(em_field_pars.E[2]));
      read_str_double("Bx",&(em_field_pars.B[0]));
      read_str_double("By",&(em_field_pars.B[1]));
      read_str_double("Bz",&(em_field_pars.B[2]));
    }
  else
    for(int i=0;i<3;i++)
      em_field_pars.E[i]=em_field_pars.B[i]=0;
}

//read parameters to measure chiral condensate
void read_chiral_cond_pars(chiral_cond_pars_t &pars)
{
  read_str_int("MeasureChiralCond",&pars.flag);
  if(pars.flag)
    {
      read_str_str("ChiralCondPath",pars.path,1024);
      read_str_double("ChiralCondInvResidue",&pars.residue);
      read_str_int("ChiralCondNHits",&pars.nhits);
    }
}

//read parameters to measure magnetization
void read_magnetization_pars(magnetization_pars_t &pars)
{
  read_str_int("MeasureMagnetization",&pars.flag);
  if(pars.flag)
    {
      read_str_str("MagnetizationPath",pars.path,1024);
      read_str_double("MagnetizationInvResidue",&pars.residue);
      read_str_int("MagnetizationNHits",&pars.nhits);
    }
}

//read parameters to measure pseudoscalar correlators
void read_pseudo_corr_pars(pseudo_corr_pars_t &pars)
{
  read_str_int("MeasurePseudoCorr",&pars.flag);
  if(pars.flag)
    {
      read_str_str("PseudoCorrPath",pars.path,1024);
      read_str_double("PseudoCorrInvResidue",&pars.residue);
      read_str_int("PseudoCorrNHits",&pars.nhits);
    }
}

//read parameters to measure all rectangles
void read_all_rect_meas_pars(all_rect_meas_pars_t &pars)
{
  read_str_int("MeasureAllRect",&pars.flag);
  if(pars.flag)
    {
      read_str_str("AllRectPath",pars.path,1024);
      
      read_str_double("HYPAlpha0",&pars.hyp_alpha0);
      read_str_double("HYPAlpha1",&pars.hyp_alpha1);
      read_str_double("HYPAlpha2",&pars.hyp_alpha2);
      
      read_str_double("APEAlpha",&pars.ape_alpha);
      read_list_of_ints("APENlevels",&pars.nape_levls,&pars.nape_iters);

      read_str_int("Tint",&pars.Tmin);
      read_int(&pars.Tmax);
      read_str_int("Dint",&pars.Dmin);
      read_int(&pars.Dmax);
    }
}

//read the theory_pars parameters of the theory
void read_theory_pars(theory_pars_t &theory_pars)
{
  //kind of action
  char gauge_action_name[1024];
  read_str_str("GaugeAction",gauge_action_name,1024);
  if(strcmp(gauge_action_name,"Wilson")==0) theory_pars.gauge_action_name=Wilson_action;
  else
    if(strcmp(gauge_action_name,"tlSym")==0) theory_pars.gauge_action_name=tlSym_action;
    else crash("unknown gauge action: %s",gauge_action_name);
  
  //beta for gauge action
  read_str_double("Beta",&theory_pars.beta);
  
  //read the number of undegenerate flavs
  read_str_int("NDiffFlavs",&(theory_pars.nflavs));
  theory_pars.quark_content=nissa_malloc("quark_content",theory_pars.nflavs,quark_content_t);
  
  //read each flav parameters
  for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
    read_quark_content(theory_pars.quark_content[iflav]);

  //additional parameters to read only if fermions are defined
  if(theory_pars.nflavs!=0)
    {
      //read stouting parameters
      read_stout_pars(theory_pars.stout_pars);

      //read electric and magnetic field
      read_em_field_pars(theory_pars.em_field_pars);
      
      //read info on pseudoscalar meson correlators, condensate and magnetization measure
      read_pseudo_corr_pars(theory_pars.pseudo_corr_pars);
      read_chiral_cond_pars(theory_pars.chiral_cond_pars);
      read_magnetization_pars(theory_pars.magnetization_pars);
    }
}
