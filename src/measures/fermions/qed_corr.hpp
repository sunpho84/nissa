#ifndef _QED_CORR_HPP
#define _QED_CORR_HPP

#include "hmc/theory_pars.hpp"
#include "fermionic_meas.hpp"

namespace nissa
{
  struct qed_corr_meas_pars_t : base_fermionic_meas_t
  {
    std::string def_path(){return "qed_corrs";}
    
    int master_fprintf(FILE *fout,bool full)
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard()||
	path!=def_path();
    }
    
    qed_corr_meas_pars_t() :
      base_fermionic_meas_t()
    {path=def_path();}
    virtual ~qed_corr_meas_pars_t(){}
  };
  
  void measure_qed_corr(eo_ptr<quad_su3> conf,theory_pars_t theory_pars,qed_corr_meas_pars_t meas_pars,int iconf,int conf_created);
}

#endif
