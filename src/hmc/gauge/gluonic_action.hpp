#ifndef _GLUONIC_ACTION_HPP
#define _GLUONIC_ACTION_HPP

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "new_types/new_types_definitions.hpp"
#include "routines/ios.hpp"

#include "Symanzik_action.hpp"
#include "Wilson_action.hpp"

namespace nissa
{
  template <class T> void gluonic_action(double *gluon_action,T conf,theory_pars_t *theory_pars,bool stag_phase_present=false)
  {
    verbosity_lv1_master_printf("Computing gauge action\n");
    
    switch(theory_pars->gauge_action_name)
      {
      case WILSON_GAUGE_ACTION:Wilson_action(gluon_action,conf,theory_pars->beta,stag_phase_present);break;
      case TLSYM_GAUGE_ACTION:tlSym_action(gluon_action,conf,theory_pars->beta,stag_phase_present);break;
      case IWASAKI_GAUGE_ACTION:Iwasaki_action(gluon_action,conf,theory_pars->beta,stag_phase_present);break;
      default:crash("Unknown action");
      }
  }
}

#endif
