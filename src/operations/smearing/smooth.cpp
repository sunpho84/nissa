#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "smooth.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"

namespace nissa
{
  //smooth a configuration until measure is due
  bool smooth_lx_conf_until_next_meas(quad_su3 *smoothed_conf,smooth_pars_t &sp,double &t,double &tnext_meas,int *dirs,int staple_min_dir)
  {
    if(sp.method==smooth_pars_t::COOLING and dirs!=all_dirs) crash("not implemented");
    
    bool finished=1;
    while(t+1e-10<tnext_meas)
      switch(sp.method)
	{
	case smooth_pars_t::COOLING: cool_lx_conf(smoothed_conf,get_sweeper(sp.cool.gauge_action));t++;finished=(t+1e-10>sp.cool.nsteps);break;
	case smooth_pars_t::STOUT: stout_smear_single_level(smoothed_conf,smoothed_conf,sp.stout.rho);t++;finished=(t+1e-10>sp.stout.nlevels);break;
	case smooth_pars_t::WFLOW: Wflow_lx_conf(smoothed_conf,sp.Wflow.dt,dirs);t+=sp.Wflow.dt;finished=(t+1e-10>sp.Wflow.T);break;
	case smooth_pars_t::HYP: hyp_smear_conf(smoothed_conf,smoothed_conf,sp.hyp.alpha0,sp.hyp.alpha1,sp.hyp.alpha2,dirs);t+=1;finished=(t+1e-10>sp.hyp.nlevels);break;
	case smooth_pars_t::APE: ape_smear_conf(smoothed_conf,smoothed_conf,sp.ape.alpha,1,dirs,staple_min_dir);t+=1;finished=(t+1e-10>sp.ape.nlevels);break;
	}
    if(not finished) tnext_meas+=sp.meas_each;
    
    return finished;
  }
  
  //smooth a configuration as imposed
  void smooth_lx_conf(quad_su3 *smoothed_conf,smooth_pars_t &sp,int *dirs,int staple_min_dir)
  {
    if(sp.method==smooth_pars_t::COOLING and dirs!=all_dirs) crash("not implemented");
    
    switch(sp.method)
      {
      case smooth_pars_t::COOLING: for(int icool=0;icool<sp.cool.nsteps;icool++) cool_lx_conf(smoothed_conf,get_sweeper(sp.cool.gauge_action));break;
      case smooth_pars_t::STOUT: for(int istout=0;istout<sp.stout.nlevels;istout++) stout_smear(smoothed_conf,smoothed_conf,&sp.stout,dirs);break;
      case smooth_pars_t::WFLOW: for(double t=0;t<sp.Wflow.T+1e-10;t+=sp.Wflow.dt) Wflow_lx_conf(smoothed_conf,sp.Wflow.dt,dirs);break;
      case smooth_pars_t::HYP: for(int ihyp=0;ihyp<sp.hyp.nlevels;ihyp++) hyp_smear_conf(smoothed_conf,smoothed_conf,sp.hyp.alpha0,sp.hyp.alpha1,sp.hyp.alpha2,dirs);break;
      case smooth_pars_t::APE: for(int iape=0;iape<sp.ape.nlevels;iape++) ape_smear_conf(smoothed_conf,smoothed_conf,sp.ape.alpha,1,dirs,staple_min_dir);break;
      }
  }
  
  std::string smooth_pars_t::get_str(bool full)
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
	if(full or meas_each!=def_meas_each()) os<<" MeasEach\t=\t"<<meas_each<<"\n";
      }
    
    return os.str();
  }
}
