#ifndef _STAG_HPP
#define _STAG_HPP

#include "hmc/theory_pars.hpp"

namespace nissa
{
  struct base_fermionic_meas_t
  {
    int each;
    int after;
    std::string path;
    double residue;
    int itheory;
    int ncopies;
    int nhits;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    virtual std::string def_path()=0;
    double def_residue(){return 1e-12;}
    int def_itheory(){return 0;}
    int def_ncopies(){return 1;}
    int def_nhits(){return 1;}
    
    int measure_is_due(int ext_itheory,int iconf)
    {return (itheory==ext_itheory)&&(each>0)&&(iconf%each==0)&&(iconf>=after);}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	each!=def_each()||
	after!=def_after()||
	residue!=def_residue()||
	itheory!=def_itheory()||
	ncopies!=def_ncopies()||
	nhits!=def_nhits();
    }
    
    base_fermionic_meas_t() :
      each(def_each()),
      after(def_after()),
      residue(def_residue()),
      itheory(def_itheory()),
      ncopies(def_ncopies()),
      nhits(def_nhits()) {}
    
    ~base_fermionic_meas_t(){};
  };
  
  void compute_fw_bw_der_mel(complex *res_fw_bw,color **left,quad_su3 **conf,int mu,color **right,complex *point_result);
  void mult_Minv(color **prop,quad_su3 **conf,quad_u1 **u1b,double m,double residue,color **source);
  void mult_Minv(color **prop,quad_su3 **conf,theory_pars_t *pars,int iflav,double residue,color **source);
}

#endif
