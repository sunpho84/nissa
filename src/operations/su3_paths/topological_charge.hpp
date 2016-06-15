#ifndef _TOPOLOGICAL_CHARGE_HPP
#define _TOPOLOGICAL_CHARGE_HPP

#include "operations/smearing/smooth.hpp"

#include "hmc/gauge/topological_action.hpp"

namespace nissa
{
  //parameters to measure topology properties
  struct top_meas_pars_t
  {
    int each;
    int after;
    std::string path;
    smooth_pars_t smooth_pars;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    std::string def_path(){return "Topo";}
    
    int is_nonstandard()
    {
      return
	each!=def_each()||
	after!=def_after()||
	path!=def_path()||
	smooth_pars.is_nonstandard();
    }
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    top_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()) {}
  };
  
  void four_leaves_point(as2t_su3 leaves_summ,quad_su3 *conf,int X);
  void measure_topology_eo_conf(top_meas_pars_t &pars,quad_su3 **unsmoothed_conf_eo,int iconf,bool conf_created);
  void measure_topology_lx_conf(top_meas_pars_t &pars,quad_su3 *unsmoothed_conf,int iconf,bool conf_created,bool preserve_unsmoothed);
  void total_topological_charge_eo_conf(double *tot_charge,quad_su3 **eo_conf);
  void total_topological_charge_lx_conf(double *tot_charge,quad_su3 *lx_conf);
  void topological_staples(quad_su3 *staples,quad_su3 *conf);
}

#endif
