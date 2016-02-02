#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "theory_pars.hpp"

namespace nissa
{
  std::string theory_pars_t::get_str(int full)
  {
    std::ostringstream os;
    //header
    if(full||is_nonstandard()) os<<"Theory\n";
    //gauge action
    if(full||(gauge_action_name!=def_gauge_action_name()))
      os<<" GaugeAction\t=\t"<<gauge_action_str_from_name(gauge_action_name)<<"\n";
    //beta
    if(full||(beta!=def_beta())) os<<" Beta\t\t=\t"<<beta<<"\n";
    //topotential_pars
    os<<topotential_pars.get_str(full)<<"\n";
    //quarks
    for(size_t i=0;i<quarks.size();i++)	os<<quarks[i].get_str(full)<<"\n";
    //stout pars
    os<<stout_pars.get_str()<<"\n";
    //global em field pars
    os<<em_field_pars.get_str(full)<<"\n";
    
    return os.str();
  }
}
