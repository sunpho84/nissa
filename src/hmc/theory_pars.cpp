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
    if(full||is_nonstandard()) os<<"\nTheory\n";
    //gauge action
    if(full||(gauge_action_name!=def_gauge_action_name()))
      os<<" GaugeAction\t=\t"<<gauge_action_str_from_name(gauge_action_name)<<"\n";
    //beta
    if(full||(beta!=def_beta())) os<<" Beta\t\t=\t"<<beta<<"\n";
    //topotential_pars
    std::string topo_str=topotential_pars.get_str(full);
    os<<topo_str;
    if(topo_str.size()) os<<"\n";
    //quarks
    for(size_t i=0;i<quarks.size();i++)
      {
	std::string quark_str=quarks[i].get_str(full);
	os<<quark_str;
	if(quark_str.size()) os<<"\n";
      }
    //stout pars
    std::string stout_str=stout_pars.get_str();
    os<<stout_str;
    if(stout_str.size()) os<<"\n";
    //global em field pars
    std::string em_field_str=em_field_pars.get_str(full);
    os<<em_field_str;
    if(em_field_str.size()) os<<"\n";
    
    return os.str();
  }
}
