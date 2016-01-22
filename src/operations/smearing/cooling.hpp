#ifndef _COOLING_HPP
#define _COOLING_HPP

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
    
    gauge_action_name_t def_gauge_action(){return WILSON_GAUGE_ACTION;}
    int def_nsteps(){return 120;}
    int def_overrelax_flag(){return 0;}
    double def_overrelax_exp(){return 0.0;}
    
    int master_fprintf(FILE *fout,bool full=false)
    {
      int nprinted=0;
      
      if(full||is_nonstandard())
	{
	  nprinted+=nissa::master_fprintf(fout,"CoolPars\n");
	  if(full||gauge_action!=def_gauge_action()) nprinted+=nissa::master_fprintf(fout," Action\t\t=\t%s\n",get_gauge_action_name_str(gauge_action).c_str());
	  if(full||nsteps!=def_nsteps()) nprinted+=nissa::master_fprintf(fout," NSteps\t\t=\t%d\n",nsteps);
	}
      
      return nprinted;
    }
    
    int is_nonstandard()
    {
      return
	gauge_action!=def_gauge_action()||
	nsteps!=def_nsteps()||
	overrelax_flag!=def_overrelax_flag()||
	overrelax_exp!=def_overrelax_exp();
    }
    
    cool_pars_t() :
      gauge_action(def_gauge_action()),
      nsteps(def_nsteps()),
      overrelax_flag(def_overrelax_flag()),
      overrelax_exp(def_overrelax_exp()) {}
  };
}

#endif
