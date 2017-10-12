#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "io/input.hpp"

#include "operations/stag/mesons.hpp"
#include "operations/stag/magnetization.hpp"
#include "operations/stag/nucleon.hpp"
#include "operations/stag/putpourri.hpp"
#include "operations/stag/rendens.hpp"
#include "operations/stag/spinpol.hpp"
#include "operations/su3_paths/topological_charge.hpp"

namespace nissa
{
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
    read_str_int("NFlows",&pars.nflows);
    read_str_double("FlowStep",&pars.dt);
  }
  
  //read and return path
  std::string read_path()
  {
    char temp[1024];
    read_str_str("Path",temp,1024);
    return temp;
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
  
  //read an ape smearing parameters
  void read_ape_pars(ape_pars_t &ape_pars)
  {
    read_str_int("ApeNLevel",&ape_pars.nlevels);
    read_str_double("ApeAlpha",&ape_pars.alpha);
  }
}
