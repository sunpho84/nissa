#ifndef _GAUGECONF_HPP
#define _GAUGECONF_HPP

#include <sstream>

#include "operations/su3_paths/gauge_sweeper.hpp"
#include "operations/smearing/smooth.hpp"

namespace nissa
{
  //Starting condition for a gauge conf
  enum start_conf_cond_t{UNSPEC_START_COND,HOT_START_COND,COLD_START_COND};
  
  //Boundary conditions
  enum boundary_cond_t{UNSPEC_BOUNDARY_COND,PERIODIC_BOUNDARY_COND,OPEN_BOUNDARY_COND};
  
  //results of a unitarity check
  struct unitarity_check_result_t
  {
    std::int64_t nbroken_links;
    double average_diff;
    double max_diff;
    
    unitarity_check_result_t () : nbroken_links(0),average_diff(0.0),max_diff(0.0) {}
  };
  
  //parameters to compute gauge observabls
  struct gauge_obs_meas_pars_t
  {
    int each;
    int after;
    std::string path;
    int meas_plaq;
    int meas_energy;
    int meas_poly;
    int use_smooth;
    smooth_pars_t smooth_pars;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    std::string def_path(){return "gauge_obs";}
    int def_meas_plaq(){return 1;}
    int def_meas_energy(){return 0;}
    int def_meas_poly(){return 1;}
    int def_use_smooth(){return 0;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
      
    int is_nonstandard()
    {
      return
	each!=def_each() or
	after!=def_after() or
	path!=def_path() or
	meas_plaq!=def_meas_plaq() or
	meas_energy!=def_meas_energy() or
	meas_poly!=def_meas_poly() or
	use_smooth!=def_use_smooth() or
	smooth_pars.is_nonstandard();
    }
    
    gauge_obs_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()),
      meas_plaq(def_meas_plaq()),
      meas_energy(def_meas_energy()),
      meas_poly(def_meas_poly()),
      use_smooth(def_use_smooth())
    {}
  };
  
  /////////////////////////////////////////////////////////////
  //to be improved
  void average_gauge_energy(double *energy,
			    const LxField<quad_su3>& conf);
  
  inline double average_gauge_energy(const LxField<quad_su3>& conf)
  {
    double energy;
    average_gauge_energy(&energy,conf);
    
    return energy;
  }
  
  void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis);
  void ac_rotate_vector(void *out,void *in,int axis,size_t bps);
  void adapt_theta(quad_su3 *conf,momentum_t& old_theta,const momentum_t& put_theta,int putonbords,int putonedges);
  void cool_lx_conf(quad_su3 *conf,gauge_sweeper_t *sweeper);
  void generate_cold_eo_conf(EoField<quad_su3>& conf);
  void generate_hot_eo_conf(EoField<quad_su3>& conf);
  
  /// Generate an identical conf
  template <typename C>
  void generate_cold_lx_conf(C& conf)
  {
    NISSA_LOC_VOL_LOOP(ivol)
      for(int mu=0;mu<NDIM;mu++)
	su3_put_to_id(conf[ivol][mu]);
    
    set_borders_invalid(conf);
  }
  
  /// Generate a random conf
  template <typename C>
  void generate_hot_lx_conf(C& conf)
  {
    if(loc_rnd_gen_inited==0)
      crash("random number generator not inited");
    
    NISSA_LOC_VOL_LOOP(ivol)
      for(int mu=0;mu<NDIM;mu++)
	su3_put_to_rnd(conf[ivol][mu],loc_rnd_gen[ivol]);
    
    set_borders_invalid(conf);
  }
  
  void heatbath_lx_conf(LxField<quad_su3>& conf,gauge_sweeper_t* sweeper,const double& beta,const int& nhits);
  
  void overrelax_lx_conf(LxField<quad_su3>& conf,gauge_sweeper_t* sweeper,int nhits);
  
  void put_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
  void rem_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
  
  /// Perform a unitarity check on a lx conf
  template <typename C>
  void unitarity_check_lx_conf(unitarity_check_result_t &result,const C& conf)
  {
    //results
    LxField<double> locAvg("locAvg");
    LxField<double> locMax("locMax");
    LxField<int64_t> locNbroken("locNbroken");
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	double a=0;
	double m=0;
	int n=0;
	
	for(int idir=0;idir<NDIM;idir++)
	{
	  const double err=
	    su3_get_non_unitariness(conf[ivol][idir]);
	  
	  //compute average and max deviation
	  a=err;
	  m=err;
	  n+=(err>1e-13);
	}
	
	locAvg[ivol]=a;
	locMax[ivol]=m;
	locNbroken[ivol]=n;
      }
    NISSA_PARALLEL_LOOP_END;
    
    glb_reduce(&result.average_diff,locAvg,locVol);
    result.average_diff/=glbVol*NDIM;
    
    glbReduce(&result.max_diff,locMax,locVol,
	      [] CUDA_DEVICE (double& res,const double& acc) INLINE_ATTRIBUTE
	      {
		if(acc>res)
		  res=acc;
	      });
    glb_reduce(&result.nbroken_links,locNbroken,locVol);
  }
  
  void unitarize_lx_conf_orthonormalizing(quad_su3 *conf);
  void unitarize_lx_conf_maximal_trace_projecting(quad_su3 *conf);
  void unitarize_eo_conf_maximal_trace_projecting(eo_ptr<quad_su3> conf);
}

#endif
