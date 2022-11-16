#ifndef _PLINE_HPP
#define _PLINE_HPP

#include "operations/smearing/smooth.hpp"

namespace nissa
{
  //pars to compute polyakov loop
  struct poly_corr_meas_pars_t
  {
    int each;
    int after;
    std::string path;
    smooth_pars_t smear_pars;
    int dir;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    std::string def_path(){return "luppoli";}
    int def_dir(){return 0;}
    
    int master_fprintf(FILE *fout,int a,bool full=false);
    
    int is_nonstandard()
    {
      return
       each!=def_each()||
       after!=def_after()||
       path!=def_path()||
       dir!=def_dir()||
       smear_pars.is_nonstandard();
    }
    
    poly_corr_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()),
      dir(def_dir())
    {}
  };
  
  void average_and_corr_polyakov_loop_lx_conf(complex tra,FILE *fout,quad_su3 *conf,int mu,int itraj);
  void average_polyakov_loop_lx_conf(complex tra,quad_su3 *conf,int mu);
  void average_polyakov_loop_eo_conf(complex tra,eo_ptr<quad_su3> eo_conf,int mu);
  // void compute_Pline_dag_internal(su3 *pline,quad_su3 *conf,int mu,int xmu_start);
  // void compute_Pline_dag_point(su3 *pline,quad_su3 *conf,int mu,coords glb_x_start);
  // void compute_Pline_dag_wall(su3 *pline,quad_su3 *conf,int mu,int xmu_start);
  // void compute_Wstat_prop_finalize(su3spinspin *prop,quad_su3 *conf,int mu,int xmu_start,su3 *pline);
  // void compute_Wstat_prop_point(su3spinspin *prop,quad_su3 *conf,int mu,coords x_start);
  // void compute_Wstat_prop_wall(su3spinspin *prop,quad_su3 *conf,int mu,int xmu_start);
  // void compute_Wstat_stoch_prop(colorspinspin *prop,quad_su3 *conf,int mu,int xmu_start,color *source);
   // void compute_stoch_Pline_dag(color *pline,quad_su3 *conf,int mu,int xmu_start,color *source);
}

#endif
