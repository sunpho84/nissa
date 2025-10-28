#ifndef _GLUONIC_ACTION_HPP
#define _GLUONIC_ACTION_HPP

#include "base/debug.hpp"

#include "Symanzik_action.hpp"
#include "Wilson_action.hpp"

namespace nissa
{
  //Gauge action
  enum gauge_action_name_t{UNSPEC_GAUGE_ACTION,WILSON_GAUGE_ACTION,TLSYM_GAUGE_ACTION,IWASAKI_GAUGE_ACTION};
  
  //convert a string into gauge action name
  inline gauge_action_name_t gauge_action_name_from_str(const char *name)
  {
    //database
    const int nact_known=3;
    gauge_action_name_t act_known[nact_known]={WILSON_GAUGE_ACTION,TLSYM_GAUGE_ACTION,IWASAKI_GAUGE_ACTION};
    const char name_known[nact_known][20]={"Wilson","tlSym","Iwasaki"};
    
    //search
    int iact=0;
    while(iact<nact_known && strcasecmp(name,name_known[iact])!=0) iact++;
    
    //check
    if(iact==nact_known) CRASH("unknown gauge action: %s",name);
    
    return act_known[iact];
  }
  
  //convert a gauge action name into a str
  inline std::string gauge_action_str_from_name(gauge_action_name_t name)
  {
    std::string res;
    
    switch(name)
      {
      case WILSON_GAUGE_ACTION:
	res="Wilson";
	break;
      case TLSYM_GAUGE_ACTION:
	res="tlSym";
	break;
      case IWASAKI_GAUGE_ACTION:
	res="Iwasaki";
	break;
      default:
	res="Unknown";
      }
    
    return res;
  }
  
  template <class T>
  double gluonic_action(T&& conf,
			const gauge_action_name_t& gauge_action_name,
			const double& beta)
  {
    double gluon_action;
    
    switch(gauge_action_name)
      {
      case WILSON_GAUGE_ACTION:
	gluon_action=Wilson_action(conf,beta);
	break;
      case TLSYM_GAUGE_ACTION:
	gluon_action=tlSym_action(conf,beta);
	break;
      case IWASAKI_GAUGE_ACTION:
	gluon_action=Iwasaki_action(conf,beta);
	break;
      default:
	CRASH("Unknown action");
      }
    
    return gluon_action;
  }
}

#endif
