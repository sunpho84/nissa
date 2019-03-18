#ifndef _SPECTRAL_PROPS_HPP
#define _SPECTRAL_PROPS_HPP

#include "eigenvalues/eigenvalues.hpp"
#include "fermionic_meas.hpp"
#include "hmc/gauge/topological_action.hpp"
#include "hmc/theory_pars.hpp"
#include "operations/smearing/smooth.hpp"

namespace nissa
{
  //parameters to measure topology properties
  struct spectr_props_meas_pars_t : base_fermionic_meas_t
  {
    std::string opname;
    int neigs;            //number of eigenvalues required
    bool minmax;
    double eig_precision; //relative precision between eigenvalues
    int wspace_size;      //size of Krylov space for the Arnoldi algorithm (it would be clipped in [2*neigs,mat_size])
    smooth_pars_t smooth_pars;

    std::string def_path(){return "spmeas_eigs_prs";}
    std::string def_opname(){return "iDst";}
    int def_neigs(){return 5;}
    bool def_minmax(){return 0;}
    double def_eig_precision(){return 1e-5;}
    int def_wspace_size(){return DEFAULT_EIGPROB_WSPACE_SIZE;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	path!=def_path() or
	opname!=def_opname() or
	neigs!=def_neigs() or
	minmax!=def_minmax() or
	eig_precision!=def_eig_precision() or
	wspace_size!=def_wspace_size() or
	smooth_pars.is_nonstandard();
    }
    
    spectr_props_meas_pars_t() :
      base_fermionic_meas_t(),
      opname(def_opname()),
      neigs(def_neigs()),
      minmax(def_minmax()),
      eig_precision(def_eig_precision()),
      wspace_size(def_wspace_size())
    {
      path=def_path();
    }
    virtual ~spectr_props_meas_pars_t(){}
  };
  
  void measure_spectral_props(quad_su3 **conf,theory_pars_t &theory_pars,spectr_props_meas_pars_t &pars, int iconf,bool conf_created);
}

#endif
