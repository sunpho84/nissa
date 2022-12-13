#ifndef _SPECTRAL_PROJECTORS_HPP
#define _SPECTRAL_PROJECTORS_HPP

#include "eigenvalues/eigenvalues.hpp"
#include "fermionic_meas.hpp"
#include "hmc/gauge/topological_action.hpp"
#include "hmc/theory_pars.hpp"
#include "operations/smearing/smooth.hpp"

namespace nissa
{
  //parameters to measure topology properties
  struct spectr_proj_meas_pars_t : base_fermionic_meas_t
  {
    int neigs;            //number of eigenvalues required
    
    double eig_precision; //relative precision between eigenvalues
    
    int wspace_size;      //size of Krylov space for the Arnoldi algorithm (it would be clipped in [2*neigs,mat_size])
    
    smooth_pars_t smooth_pars;
    
    std::string def_path() const
    {
      return "pettirosso";
    }
    
    int def_neigs() const
    {
      return 5;
    }
    
    double def_eig_precision() const
    {
      return 1e-5;
    }
    
    int def_wspace_size() const
    {
      return DEFAULT_EIGPROB_WSPACE_SIZE;
    }
    
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"MeasSpectrProj\n";
      os<<base_fermionic_meas_t::get_str(full);
      if(neigs!=def_neigs() or full) os<<" Neigs\t\t=\t"<<neigs<<"\n";
      if(eig_precision!=def_eig_precision() or full) os<<" EigPrecision\t\t=\t"<<eig_precision<<"\n";
      if(wspace_size!=def_wspace_size() or full) os<<" WSpaceSize\t\t=\t"<<wspace_size<<"\n";
      os<<smooth_pars.get_str(full);
      
      return os.str();
    }
    
    bool is_nonstandard() const
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	path!=def_path() or
	neigs!=def_neigs() or
	eig_precision!=def_eig_precision() or
	wspace_size!=def_wspace_size() or
	smooth_pars.is_nonstandard();
    }
    
    spectr_proj_meas_pars_t() :
      base_fermionic_meas_t(),
      neigs(def_neigs()),
      eig_precision(def_eig_precision()),
      wspace_size(def_wspace_size())
    {
      path=def_path();
    }
    
    virtual ~spectr_proj_meas_pars_t()
    {
    }
  };
  
  void measure_spectral_proj(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,spectr_proj_meas_pars_t &pars, int iconf,bool conf_created);
}

#endif
