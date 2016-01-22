#ifndef _WATUSSO_HPP
#define _WATUSSO_HPP

#include <vector>

#include "operations/smearing/smooth.hpp"

namespace nissa
{
  //parameters to measure flux tube
  struct watusso_meas_pars_t
  {
    int each;
    int after;
    std::string path;
    smooth_pars_t smear_pars;
    std::vector<int> sizes;
    int dmax;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    std::string def_path(){return "watusso";}
    int def_dmax(){return 10;}
    
    
    void fill_sizes_seq(int min,int step,int max)
    {for(int size=min;size<=max;size+=step) sizes.push_back(size);}
    
    int is_nonstandard()
    {
      return
	each!=def_each()||
	after!=def_after()||
	path!=def_path()||
	sizes.size()||
	dmax!=def_dmax()||
	smear_pars.is_nonstandard();
    }
    
    watusso_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()),
      dmax(def_dmax()) {}
  };
  
  void measure_watusso(watusso_meas_pars_t *pars,quad_su3 **conf,int iconf,int create_output_file);
}

#endif
