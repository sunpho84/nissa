#ifndef _MESONS_HPP
#define _MESONS_HPP

#include "stag.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  //form the mask for x (-1)^[x*(s^<+n^>)]
  inline int form_stag_op_pattern(int ispin,int itaste)
  {
    int res=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	int p=0;
	for(int nu=0;nu<mu;nu++) p+=(itaste>>nu)&1;
	for(int nu=mu+1;nu<NDIM;nu++) p+=(ispin>>nu)&1;
	p&=1;
	
	res+=(p<<mu);
      }
    
    return res;
  }
  inline int form_stag_meson_pattern_with_g5g5(int ispin,int itaste)
  {
    //add g5*g5
    ispin^=15;
    itaste^=15;
    return form_stag_op_pattern(ispin,itaste);
  }
  
  void put_stag_phases(color **source,int mask);
  
  struct meson_corr_meas_pars_t : base_fermionic_meas_t
  {
    std::vector<std::pair<int,int> > mesons;
    
    std::string def_path(){return "meson_corrs";}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	base_fermionic_meas_t::is_nonstandard() or
	mesons.size() or
	path!=def_path();
    }
    
    meson_corr_meas_pars_t() :
      base_fermionic_meas_t()
    {path=def_path();}
    virtual ~meson_corr_meas_pars_t(){}
  };
  
  void measure_meson_corr(quad_su3 **conf,theory_pars_t &tp,meson_corr_meas_pars_t &pars,int iconf,int conf_created);
}

#endif
