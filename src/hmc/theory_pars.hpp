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
    gauge_action_name_t gauge_action_name;
    double def_beta(){return 6;}
    gauge_action_name_t def_gauge_action_name(){return WILSON_GAUGE_ACTION;}
    
    std::vector<eo_ptr<quad_u1>> backfield;
    std::vector<quark_content_t> quarks;
    topotential_pars_t topotential_pars;
    stout_pars_t stout_pars;
    em_field_pars_t em_field_pars;
    
    std::string get_str(int full=false);
    int master_fprintf(FILE *fout,int full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    
    bool clover_to_be_computed()
    {
      bool to_be_computed=false;
      for(int iflav=0;iflav<nflavs();iflav++) to_be_computed|=ferm_discretiz::include_clover(quarks[iflav].discretiz);
      return to_be_computed;
    }
    
    int is_nonstandard()
    {
      return
	beta!=def_beta() or
	gauge_action_name!=def_gauge_action_name() or
	quarks.size() or
	topotential_pars.is_nonstandard() or
	stout_pars.is_nonstandard() or
	em_field_pars.is_nonstandard();
    }
    
    int nflavs()
    {return quarks.size();}
    
    theory_pars_t() :
      beta(def_beta()),
      gauge_action_name(def_gauge_action_name())
    {}
    
    void allocate_backfield();
    void init_backfield();
    void allocinit_backfield();
    void destroy_backfield();
  };
}

#endif
