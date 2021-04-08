#ifndef _GAUGE_FIXING_HPP
#define _GAUGE_FIXING_HPP

#include "geometry/geometry_eo.hpp"
#include "io/input.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"

#include <sstream>

namespace nissa
{
  struct LC_gauge_fixing_pars_t
  {
    //direction from which to start
    enum gauge_t{Landau=0,Coulomb=1};
    static inline std::string gauge_tag(gauge_t gauge)
    {
      std::string res;
      
      switch(gauge)
	{
	case Landau: res="Landau";break;
	case Coulomb: res="Coulomb";break;
	default:
	  crash("unknown gauge %d",gauge);
	}
      
      return res;
    }
    static inline gauge_t gauge_from_tag(const char *tag)
    {
      if(strcasecmp(tag,"Landau")==0) return Landau;
      if(strcasecmp(tag,"Coulomb")==0) return Coulomb;
      crash("Unknown gauge %s, use \"Landau\" or \"Coulomb\"",tag);
      return Landau;
    }
    gauge_t def_gauge() const {return Landau;}
    gauge_t gauge;
    
    //precision for minimization
    double def_target_precision() const {return 1e-14;}
    double target_precision;
    
    //max number of iterations
    int def_nmax_iterations() const {return 100000;}
    int nmax_iterations;
    
    //iterations between reunitarization
    int def_unitarize_each() const {return 100;}
    int unitarize_each;
    
    //minimization method
    enum method_t{overrelax,exponentiate};
    static inline std::string method_tag(method_t method)
    {
      std::string res;
      
      switch(method)
	{
	case overrelax: res="Overrelax";break;
	case exponentiate: res="Exponentiate";break;
	default:
	  crash("unknwon method %d",method);
	}
      
      return res;
    }
    static inline method_t method_from_tag(const char *tag)
    {
      if(strcasecmp(tag,"Overrelax")==0) return overrelax;
      if(strcasecmp(tag,"Exponentiate")==0) return exponentiate;
      crash("Unknown method %s, use \"Overrelax\" or \"Exponentiate\"",tag);
      return overrelax;
    }
    method_t def_method() const {return exponentiate;}
    method_t method;
    
    //probability to overrelax
    double def_overrelax_prob() const {return 0.9;}
    double overrelax_prob;
    
    //parameter for alpha in exp(-i alpha der /2)
    double def_alpha_exp() const {return 0.16;}
    double alpha_exp;
    
    //use or not the adaptative search of 1405.5812
    int def_use_adaptative_search() const {return 1;}
    int use_adaptative_search;
    
    //use or not the generalized cg
    int def_use_generalized_cg() const {return 1;}
    int use_generalized_cg;
    
    //use or not fft to accelerate
    int def_use_fft_acc() const {return 1;}
    int use_fft_acc;
    
    int master_fprintf(FILE *fout,int full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(int full=false)
    {
      std::ostringstream os;
      if(full or is_nonstandard())
	{
	  os<<"GaugeFix\n";
	  if(full or gauge!=def_gauge()) os<<" Gauge\t=\t"<<gauge_tag(gauge)<<"\n";
	  if(full or target_precision!=def_target_precision()) os<<" TargetPrecision\t=\t"<<target_precision<<"\n";
	  if(full or nmax_iterations!=def_nmax_iterations()) os<<" TargetPrecision\t=\t"<<nmax_iterations<<"\n";
	  if(full or unitarize_each!=def_unitarize_each()) os<<" TargetPrecision\t=\t"<<unitarize_each<<"\n";
	  if(full or method!=def_method()) os<<" Method\t=\t"<<method_tag(method)<<"\n";
	  if(full or overrelax_prob!=def_overrelax_prob()) os<<" OverrelaxProb\t=\t"<<overrelax_prob<<"\n";
	  if(full or alpha_exp!=def_alpha_exp()) os<<" AlphaExp\t=\t"<<alpha_exp<<"\n";
	  if(full or use_adaptative_search!=def_use_adaptative_search()) os<<" UseAdaptativeSearch\t=\t"<<use_adaptative_search<<"\n";
	  if(full or use_generalized_cg!=def_use_generalized_cg()) os<<" UseGeneralizedCg\t=\t"<<use_generalized_cg<<"\n";
	  if(full or use_fft_acc!=def_use_fft_acc()) os<<" UseFFTacc\t=\t"<<use_fft_acc<<"\n";
	}
      return os.str();
    }
    
    bool is_nonstandard()
    {
      return
	gauge!=def_gauge() or
	target_precision!=def_target_precision() or
	nmax_iterations!=def_nmax_iterations() or
	unitarize_each!=def_unitarize_each() or
	method!=def_method() or
	overrelax_prob!=def_overrelax_prob() or
	alpha_exp!=def_alpha_exp() or
	use_adaptative_search!=def_use_adaptative_search() or
	use_generalized_cg!=def_use_generalized_cg() or
	use_fft_acc!=def_use_fft_acc();
    }
    
    LC_gauge_fixing_pars_t() :
      gauge(def_gauge()),
      target_precision(def_target_precision()),
      nmax_iterations(def_nmax_iterations()),
      unitarize_each(def_unitarize_each()),
      method(def_method()),
      overrelax_prob(def_overrelax_prob()),
      alpha_exp(def_alpha_exp()),
      use_adaptative_search(def_use_adaptative_search()),
      use_generalized_cg(def_use_generalized_cg()),
      use_fft_acc(def_use_fft_acc())
    {}
  };
  
  //read all Landau or Coulomb gauge fixing pars
  inline void read_LC_gauge_fixing_pars(LC_gauge_fixing_pars_t &pars)
  {
    char gauge_tag[200];
    read_str_str("Gauge",gauge_tag,200);
    pars.gauge=LC_gauge_fixing_pars_t::gauge_from_tag(gauge_tag);
    read_str_double("TargetPrecision",&pars.target_precision);
    read_str_int("NMaxIterations",&pars.nmax_iterations);
    read_str_int("UnitarizeEach",&pars.unitarize_each);
    char method_tag[200];
    read_str_str("Method",method_tag,200);
    pars.method=LC_gauge_fixing_pars_t::method_from_tag(method_tag);
    switch(pars.method)
      {
      case LC_gauge_fixing_pars_t::overrelax:
	read_str_double("OverrelaxProb",&pars.overrelax_prob);
	break;
      case LC_gauge_fixing_pars_t::exponentiate:
	read_str_double("AlphaExp",&pars.alpha_exp);
	read_str_int("UseAdaptativeSearch",&pars.use_adaptative_search);
	read_str_int("UseGeneralizedCg",&pars.use_generalized_cg);
	read_str_int("UseFFTacc",&pars.use_fft_acc);
	break;
      default:
	crash("unkown method %d",pars.method);
      }
  }
  
  void gauge_transform_conf(quad_su3 *uout,su3 *g,const quad_su3 *uin);
  void gauge_transform_conf(eo_ptr<quad_su3> uout,eo_ptr<su3> g,eo_ptr<quad_su3> uin);
  
  void gauge_transform_color(eo_ptr<color> out,eo_ptr<su3> g,eo_ptr<color> in);
  
  void Landau_or_Coulomb_gauge_fix(quad_su3 *conf_out,LC_gauge_fixing_pars_t *pars,quad_su3 *conf_in);
  
  void perform_random_gauge_transform(quad_su3 *conf_out,quad_su3 *conf_in);
  void perform_random_gauge_transform(eo_ptr<quad_su3> conf_out,eo_ptr<quad_su3> conf_in);
}

#endif
