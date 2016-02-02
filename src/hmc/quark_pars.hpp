#ifndef _QUARK_PARS_HPP
#define _QUARK_PARS_HPP

#include <stdio.h>
#include <string>

#include "routines/ios.hpp"

namespace nissa
{
  //quark content
  struct quark_content_t
  {
    std::string name;
    int deg;
    bool is_stag;
    double mass;
    double kappa;
    double re_pot;
    double im_pot;
    double charge;
    
    std::string def_name(){return "quark";}
    int def_deg(){return 1;}
    bool def_is_stag(){return true;}
    double def_mass(){return 0.1;}
    double def_kappa(){return 0.125;}
    double def_re_pot(){return 0;}
    double def_im_pot(){return 0;}
    double def_charge(){return 0;}
    
    int master_fprintf(FILE *fout,int full) {return nissa::master_fprintf(fout,get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	deg!=def_deg()||
	is_stag!=def_is_stag()||
	mass!=def_mass()||
	kappa!=def_kappa()||
	re_pot!=def_re_pot()||
	im_pot!=def_im_pot()||
	charge!=def_charge();
    }
    
    quark_content_t() :
      deg(def_deg()),
      is_stag(def_is_stag()),
      mass(def_mass()),
      kappa(def_kappa()),
      re_pot(def_re_pot()),
      im_pot(def_im_pot()),
      charge(def_charge()) {}
  };
}

#endif
