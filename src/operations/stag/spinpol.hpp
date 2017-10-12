#ifndef _SPINPOL_HPP
#define _SPINPOL_HPP

#include "stag.hpp"
#include "operations/smearing/recursive_Wflower.hpp"
#include "operations/smearing/smooth.hpp"

namespace nissa
{
  struct spinpol_meas_pars_t : base_fermionic_meas_t
  {
    std::vector<std::pair<int,int> > operators;
    int nops(){return operators.size();}
    
    int use_ferm_conf_for_gluons;
    int use_adjoint_flow;
    smooth_pars_t smooth_pars;
    
    std::string def_path(){return "pollo";}
    int def_use_adjoint_flow(){return 1;}
    int def_use_ferm_conf_for_gluons(){return 0;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	nops() or
	use_ferm_conf_for_gluons!=def_use_ferm_conf_for_gluons() or
	use_adjoint_flow!=def_use_adjoint_flow() or
	path!=def_path() or
	smooth_pars.is_nonstandard();
    }
    
    spinpol_meas_pars_t() :
      base_fermionic_meas_t(),
      use_ferm_conf_for_gluons(def_use_ferm_conf_for_gluons()),
      use_adjoint_flow(def_use_adjoint_flow())
    {path=def_path();}
    
    virtual ~spinpol_meas_pars_t(){}
  };
  
  void measure_spinpol(theory_pars_t *tp,spinpol_meas_pars_t *mp,int iconf,int conf_created,quad_su3 **glu_conf);
  
  //interface
  inline void measure_spinpol(quad_su3 **ferm__ignored_conf,theory_pars_t &tp,spinpol_meas_pars_t &mp,int iconf,int conf_created,stout_pars_t &stout_pars,quad_su3 **glu_conf)
  {measure_spinpol(&tp,&mp,iconf,conf_created,glu_conf);}
}

#endif
