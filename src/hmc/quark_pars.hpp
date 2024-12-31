#ifndef _QUARK_PARS_HPP
#define _QUARK_PARS_HPP

#include <stdio.h>
#include <string>
#include <string.h>
#include <sstream>

#include "base/debug.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  namespace ferm_discretiz
  {
    const int nknown=3;
    
    enum name_t{ROOT_STAG,ROOT_TM_CLOV,OVERLAP};
    
    const name_t list[nknown]={ROOT_STAG,ROOT_TM_CLOV,OVERLAP};
    
    const char text[nknown][20]={"RootStag","RootTmClov","Overlap"};
    
    //change from name to string
    inline std::string str_from_name(const name_t& name)
    {
      std::string res;
      
      switch(name)
	{
	case ROOT_STAG:
	  res="RootStag";
	  break;
	case ROOT_TM_CLOV:
	  res="RootTMClov";
	  break;
	case OVERLAP:
	  res="Overlap";
	  break;
	}
      
      return res;
    }
    
    //string into name
    inline name_t name_from_str(const char *in)
    {
      //search
      int iname=0;
      while(iname<nknown and strcasecmp(in,text[iname])!=0)
	iname++;
      
      //check
      if(iname==nknown)
	CRASH("unknown fermion discretiz action: %s",in);
      
      return list[iname];
    }
    
    //determine if staggered or not
    inline bool is_stag(const name_t& name)
    {
      bool res=false;
      
      switch(name)
	{
	case ROOT_STAG:
	  res=true;
	  break;
	case ROOT_TM_CLOV:
	  res=false;
	  break;
	case OVERLAP:
	  res=false;
	  break;
	}
      
      return res;
    }
    
    //determine if clover or not
    inline bool include_clover(const name_t& name)
    {
      bool res=false;
      
      switch(name)
	{
	case ROOT_STAG:
	  res=false;
	  break;
	case ROOT_TM_CLOV:
	  res=true;
	  break;
	case OVERLAP:
	  res=false;
	  break;
	}
      
      return res;
    }
    
    //root needed to have 1 quarks
    inline int root_needed(const name_t& name)
    {
      int res=1;
      
      switch(name)
	{
	case ROOT_STAG:
	  res=4;
	  break;
	case ROOT_TM_CLOV:
	  res=2;
	  break;
	case OVERLAP:
	  res=1;
	  break;
	}
      
      return res;
    }
  }
  
  //quark content
  struct quark_content_t
  {
    std::string name;
    
    int deg;
    
    ferm_discretiz::name_t discretiz;
    
    double mass;
    
    double mass_overlap;
    
    double kappa;
    
    double cSW;
    
    double re_pot;
    
    double im_pot;
    
    double charge;
    
    std::string def_name() const
    {
      return "quark";
    }
    
    int def_deg() const
    {
      return 1;
    }
    
    ferm_discretiz::name_t def_discretiz() const
    {
      return ferm_discretiz::ROOT_STAG;
    }
    
    double def_mass() const
    {
      return 0.1;
    }
    
    double def_mass_overlap() const
    {
      return 1.0;
    }
    
    double def_kappa() const
    {
      return 0.125;
    }
    
    double def_cSW() const
    {
      return 0;
    }
    
    double def_re_pot() const
    {
      return 0;
    }
    
    double def_im_pot() const
    {
      return 0;
    }
    
    double def_charge() const
    {
      return 0;
    }
    
    
    int master_fprintf(FILE *fout,
		       const int& full) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"Quark\t\t=\t\""<<name.c_str()<<"\"\n";
      if(full or deg!=def_deg()) os<<" Degeneracy\t=\t"<<deg<<"\n";
      if(full or discretiz!=def_discretiz()) os<<" Discretiz\t=\t"<<ferm_discretiz::str_from_name(discretiz)<<"\n";
      if(full or mass!=def_mass()) os<<" Mass\t\t=\t"<<mass<<"\n";
      if(full or mass_overlap!=def_mass_overlap()) os<<" MassOverlap\t\t=\t"<<mass_overlap<<"\n";
      if(full or kappa!=def_kappa()) os<<" Kappa\t\t=\t"<<kappa<<"\n";
      if(full or cSW!=def_cSW()) os<<" cSW\t\t=\t"<<cSW<<"\n";
      if(full or re_pot!=def_re_pot()) os<<" RePotCh\t=\t"<<re_pot<<"\n";
      if(full or im_pot!=def_im_pot()) os<<" ImPotCh\t=\t"<<im_pot<<"\n";
      if(full or charge!=def_charge()) os<<" ElecCharge\t=\t"<<charge<<"\n";
      
      return os.str();
    }
    
    bool is_nonstandard() const
    {
      return
	deg!=def_deg() or
	discretiz!=def_discretiz() or
	mass!=def_mass() or
	mass_overlap!=def_mass_overlap() or
	kappa!=def_kappa() or
	cSW!=def_cSW() or
	re_pot!=def_re_pot() or
	im_pot!=def_im_pot() or
	charge!=def_charge();
    }
    
    quark_content_t() :
      deg(def_deg()),
      discretiz(def_discretiz()),
      mass(def_mass()),
      mass_overlap(def_mass_overlap()),
      kappa(def_kappa()),
      cSW(def_cSW()),
      re_pot(def_re_pot()),
      im_pot(def_im_pot()),
      charge(def_charge())
    {
    }
  };
}

#endif
