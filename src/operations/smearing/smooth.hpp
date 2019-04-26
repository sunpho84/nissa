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
    
    //space, time or full spacetime
    enum space_or_time_t{SPACE,TIME,SPACETIME};
    space_or_time_t space_or_time;
    space_or_time_t def_space_or_time(){return SPACETIME;}
    
    //returns the directions to smooth according to parameter
    static bool* get_dirs(space_or_time_t space_or_time)
    {
      switch(space_or_time)
      {
      case SPACE:
	return all_other_dirs[0];
	break;
      case TIME:
	return only_dir[0];
	break;
      case SPACETIME:
	return all_dirs;
	break;
      default:
	return NULL;
	crash("Unknown type");
      }
    }
    
    //returns the minimal staple direction according to parameter
    static int get_staple_min_dir(space_or_time_t space_or_time)
    {
      switch(space_or_time)
      {
      case SPACE:
	return 1;
	break;
      case TIME:
	return 1;
	break;
      case SPACETIME:
	return 0;
	break;
      default:
	return 0;
	crash("Unknown type");
      }
    }
    
    //convert a space_or_time_t into a str
    inline std::string space_or_time_str_from_name(space_or_time_t space_or_time)
    {
      switch(space_or_time)
      {
      case SPACE:
	return "Space";
	break;
      case TIME:
	return "Time";
	break;
      case SPACETIME:
	return "SpaceTime";
	break;
      default:
	return "Boh";
	crash("Unknown type");
      }
    }
    
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
    
    //return the next measure strictly after nsmooth
    int next_nsmooth_meas(int nsmooth)
    {
      return (nsmooth+meas_each_nsmooth)/meas_each_nsmooth*meas_each_nsmooth;
    }
    
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
	space_or_time!=def_space_or_time() or
	meas_each_nsmooth!=def_meas_each_nsmooth();
    }
    
    smooth_pars_t() :
      method(def_method()),
      meas_each_nsmooth(def_meas_each_nsmooth()),
      space_or_time(def_space_or_time())
    {}
  };
  
  void smooth_lx_conf_one_step(quad_su3 *smoothed_conf,smooth_pars_t &sp,bool *dirs=all_dirs,int staple_min_dir=0);
  bool smooth_lx_conf_until_next_meas(quad_su3 *smoothed_conf,smooth_pars_t &sp,int &nsmooth,bool *dirs=all_dirs,int staple_min_dir=0);
  void smooth_lx_conf(quad_su3 *smoothed_conf,smooth_pars_t &sp,bool *dirs=all_dirs,int staple_min_dir=0);
}

#endif
