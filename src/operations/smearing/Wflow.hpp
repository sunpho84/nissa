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
    int nflows;
    double dt;
    int nrecu;
    double def_nflows(){return 50;}
    double def_dt(){return 0.2;}
    int def_nrecu(){return 5;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false)
    {
      std::ostringstream os;
      
      os<<"WFlow\n";
      if(full or is_nonstandard())
	{
	  if(full or nflows!=def_nflows()) os<<" NFlows\t=\t"<<nflows<<"\n";
	  if(full or dt!=def_dt()) os<<" FlowStep\t=\t"<<dt<<"\n";
	  if(full or nrecu!=def_nrecu()) os<<" NRecu\t=\t"<<nrecu<<"\n";
	}
      
      return os.str();
    }
    
    int is_nonstandard()
    {
      return
	nflows!=def_nflows() or
	dt!=def_dt() or
	nrecu!=def_nrecu();
    }
    
    Wflow_pars_t() :
      nflows(def_nflows()),
      dt(def_dt()),
      nrecu(def_nrecu()) {}
  };
  
  void Wflow_lx_conf(quad_su3 *conf,double dt,int *dirs=all_dirs);
}

#endif
