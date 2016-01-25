#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "smooth.hpp"

namespace nissa
{
  int smooth_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    if(full||is_nonstandard())
      {
	if(full||method!=def_method()||
	   (method==COOLING&&cool_pars.is_nonstandard())||
	   (method==STOUT&&stout_pars.is_nonstandard())||
	   (method==WFLOW&&Wflow_pars.is_nonstandard())||
	   (method==APE&&ape_pars.is_nonstandard())||
	   (method==HYP&&hyp_pars.is_nonstandard()))
	  {
	    nprinted+=nissa::master_fprintf(fout," /* alternatives: Cooling, Stout, WFlow, Ape, Hyp */\n");
	    nprinted+=nissa::master_fprintf(fout," SmoothMethod\t=\t");
	  }
	switch(method)
	  {
	  case COOLING: if(cool_pars.is_nonstandard()||full) nprinted+=cool_pars.master_fprintf(fout,full);break;
	  case STOUT: if(stout_pars.is_nonstandard()||full) nprinted+=stout_pars.master_fprintf(fout,full);break;
	  case WFLOW: if(Wflow_pars.is_nonstandard()||full) nprinted+=Wflow_pars.master_fprintf(fout,full);break;
	  case APE: if(ape_pars.is_nonstandard()||full) nprinted+=ape_pars.master_fprintf(fout,full);break;
	  case HYP: if(hyp_pars.is_nonstandard()||full) nprinted+=hyp_pars.master_fprintf(fout,full);break;
	  case UNSPEC_SMOOTH_METHOD: crash("unspecified");
	  }
	if(full||meas_each!=def_meas_each()) nprinted+=nissa::master_fprintf(fout," MeasEach\t=\t%lg\n",def_meas_each());
      }
    
    return nprinted;
  }
}
