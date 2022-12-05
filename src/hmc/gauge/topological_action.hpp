#ifndef _TOPOLOGICAL_ACTION_HPP
#define _TOPOLOGICAL_ACTION_HPP

#include "new_types/metadynamics.hpp"
#include "operations/smearing/stout.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //parameters to add topological potential
  struct topotential_pars_t :
    meta_pars_t
  {
    int flag;
    
    int def_flag() const
    {
      return 0;
    }
    
    double theta;
    
    double def_theta() const
    {
      return 0.0;
    }
    
    stout_pars_t stout_pars;
    
    //methods inside measures/gauge/topological_charge.cpp
    void store_if_needed(const EoField<quad_su3>& ext_conf,
			 const int& iconf) const;
    
    //methods inside operations/su3_paths/spectral_projectors.cpp
    void store_if_needed_sp(const EoField<quad_su3>& conf,
			    const int& iconf) const;
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const;
    
    std::string get_str(const bool& full=false) const;
    
    int is_nonstandard() const
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
  
  double topodynamical_potential(const double& Q,
				 const topotential_pars_t &pars);
  
  void save_topodynamical_potential(const topotential_pars_t &pars);
  
  void load_topodynamical_potential(topotential_pars_t& pars,
				    const bool& mandatory);
  
  double topotential_action(const EoField<quad_su3>& conf,
			    const topotential_pars_t &pars);
  
  double topotential_action(const LxField<quad_su3>& lx_conf,
			    const topotential_pars_t &pars);
}


#endif
