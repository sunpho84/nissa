#ifndef _GAUGECONF_HPP
#define _GAUGECONF_HPP

#include <sstream>

#include "base/old_field.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  //Starting condition for a gauge conf
  enum start_conf_cond_t{UNSPEC_START_COND,HOT_START_COND,COLD_START_COND};
  
  //Boundary conditions
  enum boundary_cond_t{UNSPEC_BOUNDARY_COND,PERIODIC_BOUNDARY_COND,OPEN_BOUNDARY_COND};
  
  //results of a unitarity check
  struct unitarity_check_result_t
  {
    int64_t nbroken_links;
    double average_diff;
    double max_diff;
    
    unitarity_check_result_t () : nbroken_links(0),average_diff(0.0),max_diff(0.0) {}
  };
  
  /////////////////////////////////////////////////////////////
  
  double average_gauge_energy(const LxField<quad_su3>& conf);
  
  void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis);
  void ac_rotate_vector(void *out,void *in,int axis,size_t bps);
  
  /// Adapt the border condition
  void adapt_theta(LxField<quad_su3>& conf,
		   momentum_t& old_theta,
		   const momentum_t& put_theta,
		   const int& putonbords=false,
		   const int& putonedges=false);
  
  /// Put boundary conditions on the gauge conf
  void put_boundaries_conditions(LxField<quad_su3>& conf,
				 const momentum_t& theta_in_pi,
				 const int& putonbords=false,
				 const int& putonedges=false);
  
  /// Remove boundary conditions on the gauge conf
  void rem_boundaries_conditions(LxField<quad_su3>& conf,
				 const momentum_t& theta_in_pi,
				 const int& putonbords=false,
				 const int& putonedges=false);
  
  void generate_cold_eo_conf(OldEoField<quad_su3>& conf);
  void generate_hot_eo_conf(OldEoField<quad_su3>& conf);
  
  /// Generate an identical conf
  template <typename C>
  void generate_cold_lx_conf(C& conf)
  {
    NISSA_LOC_VOL_LOOP(ivol)
      for(int mu=0;mu<NDIM;mu++)
	su3_put_to_id(conf[ivol][mu]);
    
    set_borders_invalid(conf);
  }
  
  /// Perform a unitarity check on a lx conf
  template <typename C>
  void unitarity_check_lx_conf(unitarity_check_result_t &result,const C& conf)
  {
    //results
    LxField<double> locAvg("locAvg");
    LxField<double> locMax("locMax");
    LxField<int64_t> locNbroken("locNbroken");
    
    PAR(0,locVol,
	CAPTURE(TO_WRITE(locAvg),
		TO_WRITE(locMax),
		TO_WRITE(locNbroken),
		TO_READ(conf)),
	ivol,
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
	});
    
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
