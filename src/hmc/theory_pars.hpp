#ifndef _THEORY_PARS_HPP
#define _THEORY_PARS_HPP

#include "hmc/gauge/topological_action.hpp"
#include "operations/smearing/stout.hpp"
#include "hmc/backfield.hpp"
#include "hmc/gauge/gluonic_action.hpp"
#include "hmc/quark_pars.hpp"

namespace nissa
{
  //theory content
  struct theory_pars_t
  {
    double beta;
    double def_beta(){return 6;}
    
    std::vector<quad_u1**> backfield;
    std::vector<quark_content_t> quark_content;
    gauge_action_name_t gauge_action_name;
    topotential_pars_t topotential_pars;
    stout_pars_t stout_pars;
    em_field_pars_t em_field_pars;
    
    int is_nonstandard()
    {return beta!=def_beta();}
    
    int nflavs()
    {return quark_content.size();}
    
    theory_pars_t() :
      beta(def_beta()) {}
    
    void allocate_backfield();
    void init_backfield();
    void allocinit_backfield();
  };
}

#endif
