#ifndef _ALL_RECTANGLES_HPP
#define _ALL_RECTANGLES_HPP

#include "operations/smearing/smooth.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //parameters to measure all rectangles path
  struct all_rects_meas_pars_t
  {
    int each;
    int after;
    std::string path;
    smooth_pars_t temp_smear_pars;
    smooth_pars_t spat_smear_pars;
    int Tmin,Tmax,Dmin,Dmax;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    std::string def_path(){return "rectangles";}
    int def_Tmin(){return 1;}
    int def_Tmax(){return 9;}
    int def_Dmin(){return 1;}
    int def_Dmax(){return 9;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	each!=def_each() or
	after!=def_after() or
	path!=def_path() or
	Tmin!=def_Tmin() or
	Tmax!=def_Tmax() or
	Dmin!=def_Dmin() or
	Dmax!=def_Dmax() or
	temp_smear_pars.is_nonstandard() or
	spat_smear_pars.is_nonstandard();
    }
    
    all_rects_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()),
      Tmin(def_Tmin()),
      Tmax(def_Tmax()),
      Dmin(def_Dmin()),
      Dmax(def_Dmax()) {}
  };
  
  void measure_all_rectangular_paths(all_rects_meas_pars_t *pars,quad_su3  *conf,int iconf,int create_output_file);
  void measure_all_rectangular_paths_old(all_rects_meas_pars_t *pars,quad_su3  *conf,int iconf,int create_output_file);
  void measure_all_rectangular_paths(all_rects_meas_pars_t *pars,eo_ptr<quad_su3> conf,int iconf,int create_output_file);
}

#endif
