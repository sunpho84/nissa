#ifndef _MAGNETIZATION_HPP
#define _MAGNETIZATION_HPP

#include "stag.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct magnetization_meas_pars_t : base_fermionic_meas_t
  {
    std::string def_path(){return "magnetization";}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard()||
	path!=def_path();
    }
    
    magnetization_meas_pars_t() :
      base_fermionic_meas_t()
    {path=def_path();}
    virtual ~magnetization_meas_pars_t(){}
  };
  
  void magnetization(complex *magn,quad_su3 **conf,quark_content_t *quark,color **rnd,color **chi,complex *point_magn,coords *arg,int mu,int nu);
  void magnetization(complex *magn,quad_su3 **conf,int quantization,quad_u1 **u1b,quark_content_t *quark,double residue,color **rnd);
  void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,magnetization_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
