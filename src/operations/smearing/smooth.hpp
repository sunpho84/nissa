#ifndef _SMOOTH_HPP
#define _SMOOTH_HPP

#include "APE.hpp"
#include "cooling.hpp"
#include "HYP.hpp"
#include "stout.hpp"
#include "Wflow.hpp"

namespace nissa
{
  //parameters to smooth a configuration
  struct smooth_pars_t
  {
    enum method_t{COOLING,STOUT,WFLOW,APE,HYP};
    
    //basic
    method_t method;
    int meas_each_nsmooth;
    method_t def_method(){return COOLING;}
    int def_meas_each_nsmooth(){return 1;}
    
    //pars
    cool_pars_t cool;
    stout_pars_t stout;
    Wflow_pars_t Wflow;
    ape_pars_t ape;
    hyp_pars_t hyp;
    
    std::string get_method_name()
    {
      switch(method)
	{
	case COOLING: return "Cooling";break;
	case STOUT: return "Stout";break;
	case WFLOW: return "WFlow";break;
	case APE: return "Ape";break;
	case HYP: return "Hyp";break;
	default: crash("not meant to be reached");return "";
	}
    }
    
    //return the next measure
    int next_nsmooth_meas(int nsmooth)
    {return (nsmooth+meas_each_nsmooth)/meas_each_nsmooth*meas_each_nsmooth;}
    
    //returns the number of smooth
    int nsmooth()
    {
      switch(method)
	{
	case COOLING:return cool.nsteps;break;
	case STOUT:return stout.nlevels;break;
	case WFLOW:return Wflow.nflows;break;
	case APE:return ape.nlevels;break;
	case HYP:return hyp.nlevels;break;
	default:crash("not meant to be reached");return 0;
	}
    }
    
    //returns the number of measurement, without 0
    int nmeas_nonzero()
    {return nsmooth()/meas_each_nsmooth;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	method!=def_method() or
	meas_each_nsmooth!=def_meas_each_nsmooth();
    }
    
    smooth_pars_t() :
      method(def_method()),
      meas_each_nsmooth(def_meas_each_nsmooth()) {}
  };
  
  void smooth_lx_conf_one_step(quad_su3 *smoothed_conf,smooth_pars_t &sp,int *dirs=all_dirs,int staple_min_dir=0);
  bool smooth_lx_conf_until_next_meas(quad_su3 *smoothed_conf,smooth_pars_t &sp,int &nsmooth,int *dirs=all_dirs,int staple_min_dir=0);
  void smooth_lx_conf(quad_su3 *smoothed_conf,smooth_pars_t &sp,int *dirs=all_dirs,int staple_min_dir=0);
}

#endif
