#ifndef _FERMION_MEAS_HPP
#define _FERMION_MEAS_HPP

#include "routines/ios.hpp"
#include "base/random.hpp"

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
    rnd_t rnd_type;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    virtual std::string def_path()=0;
    double def_residue(){return 1e-12;}
    int def_itheory(){return 0;}
    int def_ncopies(){return 1;}
    int def_nhits(){return 1;}
    rnd_t def_rnd_type(){return RND_Z2;}
    
    int measure_is_due(int ext_itheory,int iconf)
    {
      return
	(itheory==ext_itheory) and
	(each>0) and
	(iconf%each==0) and
	(iconf>=after);
    }
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	each!=def_each() or
	after!=def_after() or
	residue!=def_residue() or
	itheory!=def_itheory() or
	ncopies!=def_ncopies() or
	rnd_type!=def_rnd_type() or
	nhits!=def_nhits();
    }
    
    base_fermionic_meas_t() :
      each(def_each()),
      after(def_after()),
      residue(def_residue()),
      itheory(def_itheory()),
      ncopies(def_ncopies()),
      nhits(def_nhits()),
      rnd_type(def_rnd_type())
    {}
    
    ~base_fermionic_meas_t(){};
  };
}

#endif
