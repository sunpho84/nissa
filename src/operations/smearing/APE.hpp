#ifndef _APE_HPP
#define _APE_HPP

#include <sstream>

#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //parameters to ape smear
  struct ape_pars_t
  {
    int nlevels;
    double alpha;
    
    int def_nlevels(){return 20;}
    double def_alpha(){return 0.5;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false)
    {
      std::ostringstream os;
      
      os<<"Ape\n";
      if(full or is_nonstandard())
	{
	  if(full or alpha!=def_alpha()) os<<" Alpha\t\t=\t"<<alpha<<"\n";
	  if(full or nlevels!=def_nlevels()) os<<" NLevels\t=\t"<<nlevels<<"\n";
	}
      
      return os.str();
    }
    
    int is_nonstandard()
    {
      return
	nlevels!=def_nlevels() or
	alpha!=def_alpha();
    }
    
    ape_pars_t() :
      nlevels(def_nlevels()),
      alpha(def_alpha()) {}
  };
  
  void ape_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep,bool *dir=all_dirs,int min_staple_dir=0);
  inline void ape_single_dir_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep,int mu,int min_staple_dir=0)
  {ape_smear_conf(smear_conf,origi_conf,alpha,nstep,only_dir[mu],min_staple_dir);}
  inline void ape_perp_dir_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep,int mu,int min_staple_dir=0)
  {ape_smear_conf(smear_conf,origi_conf,alpha,nstep,all_other_dirs[mu],min_staple_dir);}
  inline void ape_temporal_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep)
  {ape_single_dir_smear_conf(smear_conf,origi_conf,alpha,nstep,0);}
  inline void ape_spatial_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep)
  {ape_perp_dir_smear_conf(smear_conf,origi_conf,alpha,nstep,0,1);}
}

#endif
