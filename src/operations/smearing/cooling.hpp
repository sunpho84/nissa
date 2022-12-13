#ifndef _COOLING_HPP
#define _COOLING_HPP

#include <sstream>

#include "routines/ios.hpp"
#include "hmc/gauge/gluonic_action.hpp"

namespace nissa
{
   //structure to cool
  struct cool_pars_t
  {
    gauge_action_name_t gauge_action;
    
    int nsteps;
    
    int overrelax_flag;
    
    double overrelax_exp;
    
    gauge_action_name_t def_gauge_action() const
    {
      return WILSON_GAUGE_ACTION;
    }
    
    int def_nsteps() const
    {
      return 120;
    }
    
    int def_overrelax_flag() const
    {
      return 0;
    }
    
    double def_overrelax_exp() const
    {
      return 0.0;
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      os<<"Cooling\n";
      if(full or is_nonstandard())
	{
	  if(full or gauge_action!=def_gauge_action()) os<<" GaugeAction\t=\t"<<gauge_action_str_from_name(gauge_action).c_str()<<"\n";
	  if(full or nsteps!=def_nsteps()) os<<" NSteps\t\t=\t"<<nsteps<<"\n";
	}
      
      return os.str();
    }
    
    int is_nonstandard() const
    {
      return
	gauge_action!=def_gauge_action() or
	nsteps!=def_nsteps() or
	overrelax_flag!=def_overrelax_flag() or
	overrelax_exp!=def_overrelax_exp();
    }
    
    cool_pars_t() :
      gauge_action(def_gauge_action()),
      nsteps(def_nsteps()),
      overrelax_flag(def_overrelax_flag()),
      overrelax_exp(def_overrelax_exp())
    {
    }
  };
}

#endif
