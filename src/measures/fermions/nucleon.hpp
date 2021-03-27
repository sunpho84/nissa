#ifndef _NUCLEON_HPP
#define _NUCLEON_HPP

#include "fermionic_meas.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  //parameters to compute time nucleon correlator
  struct nucleon_corr_meas_pars_t : base_fermionic_meas_t
  {
    /// Smearing intensity
    double gaussSmeKappa;
    
    /// Default value for kappa
    double def_gaussSmeKappa()
    {
      return 0.0;
    }
    
    /// Number of smearing steps
    double gaussSmeNSteps;
    
    /// Default value for nsteps
    int def_gaussSmeNSteps()
    {
      return 0;
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Smearing intensity
    double apeSmeAlpha;
    
    /// Default value for alpha
    double def_apeSmeAlpha()
    {
      return 0.0;
    }
    
    /// Number of smearing steps
    double apeSmeNSteps;
    
    /// Default value for nsteps
    int def_apeSmeNSteps()
    {
      return 0;
    }
    
    /////////////////////////////////////////////////////////////////
    
    std::string def_path(){return "nucleon_corrs";}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	path!=def_path() or
	gaussSmeKappa!=def_gaussSmeKappa() or
	gaussSmeNSteps!=def_gaussSmeNSteps() or
	apeSmeAlpha!=def_apeSmeAlpha() or
	apeSmeNSteps!=def_apeSmeNSteps();
    }
    
    nucleon_corr_meas_pars_t() :
      base_fermionic_meas_t(),
      gaussSmeKappa(def_gaussSmeKappa()),
      gaussSmeNSteps(def_gaussSmeNSteps()),
      apeSmeAlpha(def_apeSmeAlpha()),
      apeSmeNSteps(def_apeSmeNSteps())
    {
      path=def_path();
    }
    virtual ~nucleon_corr_meas_pars_t(){}
  };
  
  void measure_nucleon_corr(eo_ptr<quad_su3> conf,theory_pars_t theory_pars,nucleon_corr_meas_pars_t meas_pars,int iconf,int conf_created);
}

#endif
