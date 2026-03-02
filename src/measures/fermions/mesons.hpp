#ifndef _MESONS_HPP
#define _MESONS_HPP

#include "fermionic_meas.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct meson_corr_meas_pars_t :
    base_fermionic_meas_t
  {
    std::vector<std::pair<int,int>> mesons;
    
    int dir;
    double gauss_kappa;
    int gauss_niter_src;
    int gauss_niter_snk;
    double ape_alpha;
    int ape_niter;

    double def_ape_alpha() const { return 0.5; }
    int def_ape_niter() const { return 20; }
    
    /// Time direction by defaul
    int def_dir() const
    {
      return 0;
    }

    // Smearing defaults: disabled
    double def_gauss_kappa() const { return 0.0; }
    int def_gauss_niter_src() const { return 0; }
    int def_gauss_niter_snk() const { return 0; }
    
    std::string def_path() const
    {
      return "meson_corrs";
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full=false)
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"MeasMesonCorrs\n";
      os<<base_fermionic_meas_t::get_str(full);

      if(dir!=def_dir() or full)
	os<<" Dir\t\t=\t"<<dir<<"\n";
      
      if(mesons.size() or full)
	{
	  os<<" Operators\t=\t{";
	  for(size_t i=0;i<mesons.size();i++)
	    {
	      os<<"("<<mesons[i].first<<","<<mesons[i].second<<")";
	      if(i!=mesons.size()-1) os<<",";
	    }
	  os<<"}\n";
	}

    if(ape_alpha!=def_ape_alpha() or full)
      os<<" ApeAlpha\t=\t"<<ape_alpha<<"\n";
    if(ape_niter!=def_ape_niter() or full)
      os<<" ApeNiter\t=\t"<<ape_niter<<"\n";
    if(gauss_kappa!=def_gauss_kappa() or full)
      os<<" GaussKappa\t=\t"<<gauss_kappa<<"\n";
    if(gauss_niter_src!=def_gauss_niter_src() or full)
      os<<" GaussNiterSrc\t=\t"<<gauss_niter_src<<"\n";
    if(gauss_niter_snk!=def_gauss_niter_snk() or full)
      os<<" GaussNiterSnk\t=\t"<<gauss_niter_snk<<"\n";
      
      return os.str();
    }
    
    int is_nonstandard() const
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	mesons.size() or
	dir!=def_dir() or
  gauss_kappa!=def_gauss_kappa() or
  gauss_niter_src!=def_gauss_niter_src() or
  gauss_niter_snk!=def_gauss_niter_snk() or
  ape_alpha!=def_ape_alpha() or
  ape_niter!=def_ape_niter() or
	path!=def_path();
    }
    
    meson_corr_meas_pars_t() :
      base_fermionic_meas_t(),
      dir(def_dir()),
      gauss_kappa(def_gauss_kappa()),
      gauss_niter_src(def_gauss_niter_src()),
      gauss_niter_snk(def_gauss_niter_snk()),
      ape_alpha(def_ape_alpha()),
      ape_niter(def_ape_niter())
    {
      path=def_path();
    }
    
    virtual ~meson_corr_meas_pars_t()
    {
    }
  };
  
  void measure_meson_corr(const EoField<quad_su3>& conf,
			  const theory_pars_t& tp,
			  const meson_corr_meas_pars_t& pars,
			  const int& iconf,
			  const int& conf_created);
}

#endif
