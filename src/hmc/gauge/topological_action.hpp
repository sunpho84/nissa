#ifndef _TOPOLOGICAL_ACTION_HPP
#define _TOPOLOGICAL_ACTION_HPP

#include "new_types/metadynamics.hpp"
#include "operations/smearing/stout.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //parameters to add topological potential
  struct topotential_pars_t : meta_pars_t
  {
    int flag;
    double theta;
    stout_pars_t stout_pars;
    int def_flag(){return 0;}
    double def_theta(){return 0;}
    
    //methods inside opearations/su3_paths/topological_charge.cpp
    void store_if_needed(quad_su3 **conf,int iconf);
    
    int master_fprintf(FILE *fout,bool full=false)
    {
          int nprinted=0;
	  const char name_known[3][10]={"None","","Meta"};
	  if(full||flag!=def_flag()) nprinted+=nissa::master_fprintf(fout,"TopoPotential\t=\t%s\n",name_known[flag]);
	  switch(flag)
	    {
	    case 0:break;
	    case 1:nprinted+=nissa::master_fprintf(fout,"Theta\t\t%lg\n",theta);break;
	    case 2:
	      nprinted+=meta_pars_t::master_fprintf(fout);
	      stout_pars.master_fprintf(fout);
	      break;
	    }
	  
	  return nprinted;
    }
    
    int is_nonstandard()
    {
      return
	flag!=def_flag()||
	theta!=def_theta();
    }
    
    topotential_pars_t() :
      meta_pars_t(),
      flag(def_flag()),
      theta(def_theta()){}
  };

  double topodynamical_potential(double Q,topotential_pars_t &pars);
  void save_topodynamical_potential(topotential_pars_t &pars);
  void load_topodynamical_potential(topotential_pars_t &pars,bool mandatory);
  double topotential_action(quad_su3 **ext_conf,topotential_pars_t &pars);
}

#endif
