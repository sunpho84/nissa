#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "smooth.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"

namespace nissa
{
  //smooth a configuration for one step
  void smooth_lx_conf_one_step(LxField<quad_su3>& smoothed_conf,
			       const smooth_pars_t &sp,
			       const WhichDirs& dirs,
			       const int& staple_min_dir)
  {
    crash("reimplement"); //reimplement
    // verbosity_lv3_master_printf("smoothing one step\n");
    // switch(sp.method)
    //   {
    //   case smooth_pars_t::COOLING: cool_lx_conf(smoothed_conf,get_sweeper(sp.cool.gauge_action));break;
    //   case smooth_pars_t::STOUT: stout_smear_single_level(smoothed_conf,smoothed_conf,sp.stout.rho,dirs);;break;
    //   case smooth_pars_t::WFLOW: Wflow_lx_conf(smoothed_conf,sp.Wflow.dt,dirs);break;
    //   case smooth_pars_t::HYP: hyp_smear_conf(smoothed_conf,smoothed_conf,sp.hyp.alpha0,sp.hyp.alpha1,sp.hyp.alpha2,dirs);break;
    //   case smooth_pars_t::APE: ape_smear_conf(smoothed_conf,smoothed_conf,sp.ape.alpha,1,dirs,staple_min_dir);break;
    //   }
  }
  
  /// Smooth a configuration until measure is due
  bool smooth_lx_conf_until_next_meas(LxField<quad_su3>& smoothed_conf,
				      const smooth_pars_t &sp,
				      int &nsmooth,
				      const WhichDirs& dirs,
				      const int& staple_min_dir)
  {
    bool diff=false;
    for(int mu=0;mu<NDIM;mu++)
      diff|=dirs[mu]!=allDirs[mu];
    
    if(sp.method==smooth_pars_t::COOLING and diff) crash("not implemented");
    
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
  void smooth_lx_conf(quad_su3 *smoothed_conf,smooth_pars_t &sp,const WhichDirs& dirs,int staple_min_dir)
  {
    crash("reimplement");
    
    // for(int ismooth=0;ismooth<sp.nsmooth();ismooth++)
    //   smooth_lx_conf_one_step(smoothed_conf,sp,dirs,staple_min_dir);
  }
}
