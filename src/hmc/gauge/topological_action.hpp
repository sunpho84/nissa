#ifndef _TOPOLOGICAL_ACTION_HPP
#define _TOPOLOGICAL_ACTION_HPP

#include "new_types/metadynamics.hpp"
#include "operations/smearing/stout.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //parameters to add topological potential
  struct topotential_pars_t : meta_pars_t
  {
    int flag;
    double theta;
    stout_pars_t stout_pars;
    int def_flag(){return 0;}
    double def_theta(){return 0.0;}
    
    //methods inside measures/gauge/topological_charge.cpp
    void store_if_needed(eo_ptr<quad_su3> conf,int iconf);
    
    //methods inside operations/su3_paths/spectral_projectors.cpp
    void store_if_needed_sp(eo_ptr<quad_su3> conf,int iconf);
    
    int master_fprintf(FILE *fout,bool full=false);
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	meta_pars_t::is_nonstandard() or
	flag!=def_flag() or
	theta!=def_theta() or
	stout_pars.is_nonstandard();
    }
    
    topotential_pars_t() :
      meta_pars_t(),
      flag(def_flag()),
      theta(def_theta()){}
  };
  
  double topodynamical_potential(double Q,topotential_pars_t &pars);
  void save_topodynamical_potential(topotential_pars_t &pars);
  void load_topodynamical_potential(topotential_pars_t &pars,bool mandatory);
  double topotential_action(eo_ptr<quad_su3> ext_conf,topotential_pars_t &pars);
  double topotential_action(quad_su3 *lx_conf,topotential_pars_t &pars);
}


#endif
