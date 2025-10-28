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
    
    method_t def_method() const
    {
      return COOLING;
    }
    
    int def_meas_each_nsmooth() const
    {
      return 1;
    }
    
    //space, time or full spacetime
    enum space_or_time_t{SPACE,TIME,SPACETIME};
    
    space_or_time_t space_or_time;
    
    space_or_time_t def_space_or_time() const
    {
      return SPACETIME;
    }
    
    //returns the directions to smooth according to parameter
    static WhichDirs get_dirs(const space_or_time_t& space_or_time)
    {
      WhichDirs res={};
      
      switch(space_or_time)
	{
	case SPACE:
	  res=allOtherDirs[0];
	  break;
	case TIME:
	  res=onlyDir[0];
	  break;
	case SPACETIME:
	  res=allDirs;
	  break;
	default:
	  CRASH("Unknown type");
	}
      
      return res;
    }
    
    //returns the minimal staple direction according to parameter
    static int get_staple_min_dir(const space_or_time_t& space_or_time)
    {
      int res=0;
      
      switch(space_or_time)
	{
	case SPACE:
	  res=1;
	  break;
	case TIME:
	  res=1;
	  break;
	case SPACETIME:
	  res=0;
	  break;
	default:
	  res=0;
	  CRASH("Unknown type");
	}
      
      return res;
    }
    
    //convert a space_or_time_t into a str
    inline std::string space_or_time_str_from_name(const space_or_time_t& space_or_time) const
    {
      std::string res;
      
      switch(space_or_time)
      {
      case SPACE:
	res="Space";
	break;
      case TIME:
	res="Time";
	break;
      case SPACETIME:
	res="SpaceTime";
	break;
      default:
	res="Boh";
	CRASH("Unknown type");
      }
      
      return res;
    }
    
    //pars
    cool_pars_t cool;
    
    stout_pars_t stout;
    
    Wflow_pars_t Wflow;
    
    ape_pars_t ape;
    
    hyp_pars_t hyp;
    
    std::string get_method_name() const
    {
      std::string res;
      
      switch(method)
	{
	case COOLING:
	  res="Cooling";
	  break;
	case STOUT:
	  res="Stout";
	  break;
	case WFLOW:
	  res="WFlow";
	  break;
	case APE:
	  res="Ape";
	  break;
	case HYP:
	  res="Hyp";
	  break;
	default:
	  CRASH("not meant to be reached");
	  res="";
	}
      
      return res;
    }
    
    /// return the next measure strictly after nsmooth
    int next_nsmooth_meas(int nsmooth) const
    {
      return (nsmooth+meas_each_nsmooth)/meas_each_nsmooth*meas_each_nsmooth;
    }
    
    //returns the number of smooth
    int nsmooth() const
    {
      int res=0;
      
      switch(method)
	{
	case COOLING:
	  res=cool.nsteps;
	  break;
	case STOUT:
	  res=stout.nlevels;
	  break;
	case WFLOW:
	  res=Wflow.nflows;
	  break;
	case APE:
	  res=ape.nlevels;
	  break;
	case HYP:
	  res=hyp.nlevels;
	  break;
	default:
	  CRASH("not meant to be reached");
	  res=0;
	}
      
      return res;
    }
    
    //returns the number of measurement, without 0
    int nmeas_nonzero() const
    {
      return nsmooth()/meas_each_nsmooth;
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      if(full or is_nonstandard())
	{
	  if(full or method!=def_method() or
	     (method==COOLING and cool.is_nonstandard()) or
	     (method==STOUT and stout.is_nonstandard()) or
	     (method==WFLOW and Wflow.is_nonstandard()) or
	     (method==APE and ape.is_nonstandard()) or
	     (method==HYP and hyp.is_nonstandard()))
	    {
	      os<<" SmoothMethod\t=\t";
	      switch(method)
		{
		case COOLING: os<<cool.get_str(full);break;
		case STOUT: os<<stout.get_str(full);break;
		case WFLOW: os<<Wflow.get_str(full);break;
		case APE: os<<ape.get_str(full);break;
		case HYP: os<<hyp.get_str(full);break;
		}
	      //os<<" /* alternatives: Cooling, Stout, WFlow, Ape, Hyp */\n";
	    }
	  if(full or space_or_time!=def_space_or_time()) os<<" SpaceOrTime\t=\t"<<space_or_time_str_from_name(space_or_time)<<"\n";
	  if(full or meas_each_nsmooth!=def_meas_each_nsmooth()) os<<" MeasEachNSmooth\t=\t"<<meas_each_nsmooth<<"\n";
	}
      
      return os.str();
    }
    
    int is_nonstandard() const
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
    {
    }
  };
  
  void smooth_lx_conf_one_step(LxField<quad_su3>& smoothed_conf,
			       const smooth_pars_t &sp,
			       const WhichDirs& dirs=allDirs,
			       const int& staple_min_dir=0);
  
  bool smooth_lx_conf_until_next_meas(LxField<quad_su3>& smoothed_conf,
				      const smooth_pars_t &sp,
				      int &nsmooth,
				      const WhichDirs& dirs,
				      const int& staple_min_dir=0);
  
  void smooth_lx_conf(quad_su3 *smoothed_conf,smooth_pars_t &sp,const WhichDirs& dirs=allDirs,int staple_min_dir=0);
}

#endif
