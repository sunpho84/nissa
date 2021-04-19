#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "smooth.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"

namespace nissa
{
  //smooth a configuration for one step
  void smooth_lx_conf_one_step(quad_su3 *smoothed_conf,smooth_pars_t &sp,const Coords<bool>& dirs,const Direction& staple_min_dir)
  {
    verbosity_lv3_master_printf("smoothing one step\n");
    switch(sp.method)
      {
      case smooth_pars_t::COOLING: cool_lx_conf(smoothed_conf,get_sweeper(sp.cool.gauge_action));break;
      case smooth_pars_t::STOUT: stout_smear_single_level(smoothed_conf,smoothed_conf,sp.stout.rho,dirs);;break;
      case smooth_pars_t::WFLOW: Wflow_lx_conf(smoothed_conf,sp.Wflow.dt,dirs);break;
      case smooth_pars_t::HYP: hyp_smear_conf(smoothed_conf,smoothed_conf,sp.hyp.alpha0,sp.hyp.alpha1,sp.hyp.alpha2,dirs);break;
      case smooth_pars_t::APE: ape_smear_conf(smoothed_conf,smoothed_conf,sp.ape.alpha,1,dirs,staple_min_dir);break;
      }
  }
  
  //smooth a configuration until measure is due
  bool smooth_lx_conf_until_next_meas(quad_su3 *smoothed_conf,smooth_pars_t &sp,int &nsmooth,const Coords<bool>& dirs,const Direction& staple_min_dir)
  {
    crash("fix the check below");
    //if(sp.method==smooth_pars_t::COOLING and dirs!=all_dirs) crash("not implemented");
    
    const int next_nsmooth_meas=sp.next_nsmooth_meas(nsmooth);
    const bool finished=(next_nsmooth_meas>sp.nsmooth());
    
    if(not finished)
      while(nsmooth<next_nsmooth_meas)
	{
	  verbosity_lv3_master_printf("smoothing %d to %d\n",nsmooth,next_nsmooth_meas);
	  smooth_lx_conf_one_step(smoothed_conf,sp,dirs,staple_min_dir);
	  nsmooth++;
	}
    
    return finished;
  }
  
  //smooth a configuration as imposed
  void smooth_lx_conf(quad_su3 *smoothed_conf,smooth_pars_t &sp,const Coords<bool>& dirs,const Direction& staple_min_dir)
  {
    for(int ismooth=0;ismooth<sp.nsmooth();ismooth++)
      smooth_lx_conf_one_step(smoothed_conf,sp,dirs,staple_min_dir);
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
	if(full or space_or_time!=def_space_or_time()) os<<" SpaceOrTime\t=\t"<<space_or_time_str_from_name(space_or_time)<<"\n";
	if(full or meas_each_nsmooth!=def_meas_each_nsmooth()) os<<" MeasEachNSmooth\t=\t"<<meas_each_nsmooth<<"\n";
      }
    
    return os.str();
  }
}
