#ifndef _GAUGECONF_HPP
#define _GAUGECONF_HPP

#include <sstream>

#include "operations/su3_paths/gauge_sweeper.hpp"

namespace nissa
{
  //Starting condition for a gauge conf
  enum start_conf_cond_t{UNSPEC_START_COND,HOT_START_COND,COLD_START_COND};
  
  //Boundary conditions
  enum boundary_cond_t{UNSPEC_BOUNDARY_COND,PERIODIC_BOUNDARY_COND,OPEN_BOUNDARY_COND};
  
  //results of a unitarity check
  struct unitarity_check_result_t
  {
    int nbroken_links;
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
    
    int def_itheory(){return 0;}
    int def_each(){return 1;}
    int def_after(){return 0;}
    std::string def_path(){return "gauge_obs";}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false)
    {
      std::ostringstream os;
      
      os<<"MeasPlaqPol\n";
      if(each!=def_each()||full) os<<" Each\t\t=\t"<<each<<"\n";
      if(after!=def_after()||full) os<<" After\t\t=\t"<<after<<"\n";
      if(path!=def_path()||full) os<<" Path\t\t=\t\""<<path.c_str()<<"\"\n";
      
      return os.str();
    }
    
    int is_nonstandard()
    {
      return
	each!=def_each()||
	after!=def_after()||
	path!=def_path();
    }
    
    gauge_obs_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()) {}
  };
  
  /////////////////////////////////////////////////////////////
  
  inline int check_add_square_staple(int *isquare_staples_to_ask,int &nsquare_staple_to_ask,int ivol,int dir,int verse,int iter);
  void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis);
  void ac_rotate_vector(void *out,void *in,int axis,size_t bps);
  void adapt_theta(quad_su3 *conf,double *old_theta,double *put_theta,int putonbords,int putonedges);
  void cool_lx_conf(quad_su3 *conf,gauge_sweeper_t *sweeper);
  void generate_cold_eo_conf(quad_su3 **conf);
  void generate_hot_eo_conf(quad_su3 **conf);
  void generate_cold_lx_conf(quad_su3 *conf);
  void generate_hot_lx_conf(quad_su3 *conf);
  void heatbath_lx_conf(quad_su3 *conf,gauge_sweeper_t *sweeper,double beta,int nhits);
  void overrelax_lx_conf(quad_su3 *conf,gauge_sweeper_t *sweeper,int nhits);
  void put_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
  void rem_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
  void unitarity_check_lx_conf(unitarity_check_result_t &result,quad_su3 *conf);
  void unitarize_lx_conf_orthonormalizing(quad_su3 *conf);
  void unitarize_lx_conf_maximal_trace_projecting(quad_su3 *conf);
  void unitarize_eo_conf_maximal_trace_projecting(quad_su3 **conf);
}

#endif
