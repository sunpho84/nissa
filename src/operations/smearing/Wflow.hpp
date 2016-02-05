#ifndef _WFLOW_HPP
#define _WFLOW_HPP

#include <sstream>

#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //structure to Wilson flow
  struct Wflow_pars_t
  {
    double T;
    double dt;
    double def_T(){return 10;}
    double def_dt(){return 0.2;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false)
    {
      std::ostringstream os;
      
      os<<"WFlow\n";
      if(full||is_nonstandard())
	{
	  if(full||T!=def_T()) os<<" FlowTime\t=\t"<<T<<"\n";
	  if(full||dt!=def_dt()) os<<" InteStep\t=\t"<<dt<<"\n";
	}
      
      return os.str();
    }
    
    int is_nonstandard()
    {
      return
	T!=def_T()||
	dt!=def_dt();
    }
    
    Wflow_pars_t() :
      T(def_T()),
      dt(def_dt()) {}
  };
  
  void Wflow_lx_conf(quad_su3 *conf,double dt,int *dirs=all_dirs);
}

#endif
