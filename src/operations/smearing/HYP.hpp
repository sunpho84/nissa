#ifndef _HYP_HPP
#define _HYP_HPP

#include <sstream>

#include "geometry/geometry_lx.hpp"
#include "routines/ios.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  //parameters to perform hyp
  struct hyp_pars_t
  {
    double alpha0;
    double alpha1;
    double alpha2;
    
    int nlevels;
    
    double def_alpha0() {return 1;}
    double def_alpha1() {return 1;}
    double def_alpha2() {return 0.5;}
    int def_nlevels() {return 1;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false)
    {
      std::ostringstream os;
      
      os<<"Hyp\n";
      if(full or is_nonstandard())
	{
	  if(full or alpha0!=def_alpha0() or alpha1!=def_alpha1() or alpha2!=def_alpha2()) os<<" Alphas\t\t=\t{"<<alpha0<<","<<alpha1<<","<<alpha2<<"}\n";
	  if(full or nlevels!=def_nlevels()) os<<" NLevels\t=\t"<<nlevels<<"\n";
	}
      
      return os.str();
    }
    
    int is_nonstandard()
    {
      return
	alpha0!=def_alpha0() or
	alpha1!=def_alpha1() or
	alpha2!=def_alpha2() or
	nlevels!=def_nlevels();
    }
    hyp_pars_t() :
      alpha0(def_alpha0()),
      alpha1(def_alpha1()),
      alpha2(def_alpha2()),
      nlevels(def_nlevels()) {}
  };
  
  void hyp_smear_conf(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2,const Coords<bool>& dirs=all_dirs);
}

#endif
