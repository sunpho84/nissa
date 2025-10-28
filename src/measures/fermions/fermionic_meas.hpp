#ifndef _FERMION_MEAS_HPP
#define _FERMION_MEAS_HPP

#include <sstream>

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
    
    int def_each() const
    {
      return 1;
    }
    
    int def_after() const
    {
      return 0;
    }
    
    virtual std::string def_path() const =0;
    
    double def_residue() const
    {
      return 1e-12;
    }
    
    int def_itheory() const
    {
      return 0;
    }
    
    int def_ncopies() const
    {
      return 1;
    }
    
    int def_nhits() const
    {
      return 1;
    }
    
    rnd_t def_rnd_type() const
    {
      return RND_Z2;
    }
    
    int measure_is_due(const int& ext_itheory,
		       const int& iconf) const
    {
      return
	(itheory==ext_itheory) and
	(each>0) and
	(iconf%each==0) and
	(iconf>=after);
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool full=false) const
    {
      std::ostringstream os;
      
      if(each!=def_each() or full) os<<" Each\t\t=\t"<<each<<"\n";
      if(after!=def_after() or full) os<<" After\t\t=\t"<<after<<"\n";
      if(path!=def_path() or full) os<<" Path\t\t=\t\""<<path.c_str()<<"\"\n";
      if(residue!=def_residue() or full) os<<" Residue\t=\t"<<residue<<"\n";
      if(ncopies!=def_ncopies() or full) os<<" NCopies\t=\t"<<ncopies<<"\n";
      if(itheory!=def_itheory() or full) os<<" ITheory\t=\t"<<itheory<<"\n";
      if(rnd_type!=def_rnd_type() or full) os<<" NoiseType\t=\t"<<rnd_t_str[rnd_type]<<"\n";
      if(nhits!=def_nhits() or full) os<<" NHits\t\t=\t"<<nhits<<"\n";
      
      return os.str();
    }
    
    bool is_nonstandard() const
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
    {
    }
  };
}

#endif
