#ifndef _ELLESETTETE_HPP
#define _ELLESETTETE_HPP

#include "hmc/theory_pars.hpp"

#include "fermionic_meas.hpp"
#include "stag.hpp"

namespace nissa
{
  using namespace stag;
  struct ellesettete_meas_pars_t : base_fermionic_meas_t
  {
    enum MethodType{ANALYTICAL,NUMERICAL};
    
    int max_order;
    MethodType method;
    double epsilon; // precision for numerical method
    GAMMA_INT taste_channel;

    std::string def_path(){return "ellesettete";}
    MethodType def_method(){return ANALYTICAL;}
    double def_epsilon(){return 1e-6;}
    GAMMA_INT def_taste_channel(){return GAMMA_INT::GAMMA_1;}
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    static const char* getAnalyticalNumericalTag(const MethodType& method)
    {
      switch(method)
	{
	case ANALYTICAL:
	  return "Numerical";
	  break;
	case NUMERICAL:
	  return "Analytical";
	  break;
	default:
	  crash("Unknown method %d",method);
	  break;
	}
      
      return "";
    }

    inline std::string gamma_int_to_str(GAMMA_INT gamma) const
    {
        switch(gamma)
        {
            case GAMMA_INT::GAMMA_IDENTITY: return "GammaID";
            case GAMMA_INT::GAMMA_5: return "Gamma5";
            case GAMMA_INT::GAMMA_1: return "Gamma1";
            case GAMMA_INT::GAMMA_2: return "Gamma2";
            case GAMMA_INT::GAMMA_3: return "Gamma3";
            case GAMMA_INT::GAMMA_5_GAMMA_1: return "Ax1";
            case GAMMA_INT::GAMMA_5_GAMMA_2: return "Ax2";
            case GAMMA_INT::GAMMA_5_GAMMA_3: return "Ax3";
            default: 
              crash("Unknown type or not implemented");
              return "Unknown";
        }
    }

    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	path!=def_path() or
	method!=def_method() or
	epsilon!=def_epsilon() or
  taste_channel!=def_taste_channel();
    }
    
    ellesettete_meas_pars_t() :
      base_fermionic_meas_t(),
      method(def_method()),
      epsilon(def_epsilon()),
      taste_channel(def_taste_channel())
    {
      path=def_path();
    }
    virtual ~ellesettete_meas_pars_t(){}
  };
  
  void measure_ellesettete(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,ellesettete_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
