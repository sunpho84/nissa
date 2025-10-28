#ifndef _SPINPOL_HPP
#define _SPINPOL_HPP

#include "fermionic_meas.hpp"
#include "hmc/theory_pars.hpp"
#include "operations/smearing/recursive_Wflower.hpp"
#include "operations/smearing/smooth.hpp"

namespace nissa
{
  struct spinpol_meas_pars_t :
    base_fermionic_meas_t
  {
    std::vector<std::pair<int,int> > operators;

    int nops() const
    {
      return operators.size();
    }
    
    int use_ferm_conf_for_gluons;
    
    int use_adjoint_flow;
    
    smooth_pars_t smooth_pars;
    
    std::string def_path() const
    {
      return "pollo";
    }
    
    int def_use_adjoint_flow() const
    {
      return 1;
    }
    
    int def_use_ferm_conf_for_gluons() const
    {
      return 0;
    }
    
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"MeasSpinPol\n";
      os<<base_fermionic_meas_t::get_str(full);
      if(use_adjoint_flow!=def_use_adjoint_flow() or full) os<<" UseAdjointFlow\t=\t"<<use_adjoint_flow<<"\n";
      if(use_ferm_conf_for_gluons!=def_use_ferm_conf_for_gluons() or full) os<<" UseFermConfForGluons\t=\t"<<use_ferm_conf_for_gluons<<"\n";
      if(operators.size())
	{
	  os<<" Operators\t=\t{";
	  for(size_t i=0;i<operators.size();i++)
	    {
	      os<<"("<<operators[i].first<<","<<operators[i].second<<")";
	      if(i!=operators.size()-1) os<<",";
	    }
	  os<<"}\n";
	}
      os<<smooth_pars.get_str(full);
    
      return os.str();
    }
    
    bool is_nonstandard() const
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
    {
      path=def_path();
    }
    
    virtual ~spinpol_meas_pars_t()
    {
    }
  };
  
  void measure_spinpol(theory_pars_t *tp,spinpol_meas_pars_t *mp,int iconf,int conf_created,eo_ptr<quad_su3> glu_conf);
  
  //interface
  inline void measure_spinpol(eo_ptr<quad_su3> ferm__ignored_conf,theory_pars_t &tp,spinpol_meas_pars_t &mp,int iconf,int conf_created,stout_pars_t &stout_pars,eo_ptr<quad_su3> glu_conf)
  {measure_spinpol(&tp,&mp,iconf,conf_created,glu_conf);}
}

#endif
