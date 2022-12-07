#ifndef _APE_HPP
#define _APE_HPP

#include <sstream>

#include <base/field.hpp>
#include <geometry/geometry_lx.hpp>
#include <new_types/su3.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  //parameters to ape smear
  struct ape_pars_t
  {
    int nlevels;
    
    double alpha;
    
    int def_nlevels() const
    {
      return 20;
    }
    
    double def_alpha() const
    {
      return 0.5;
    }
    
    int master_fprintf(FILE *fout,
		       const bool full) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
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
    
    int is_nonstandard() const
    {
      return
	nlevels!=def_nlevels() or
	alpha!=def_alpha();
    }
    
    ape_pars_t() :
      nlevels(def_nlevels()),
      alpha(def_alpha())
    {
    }
  };
  
  void ape_smear_conf(LxField<quad_su3>& smear_conf,
		      LxField<quad_su3> origi_conf,
		      const double& alpha,
		      const int& nstep,
		      const which_dir_t& dirs=all_dirs,
		      const int& min_staple_dir=0);
  
  inline void ape_single_dir_smear_conf(LxField<quad_su3>& smear_conf,
					const LxField<quad_su3>& origi_conf,
					const double& alpha,
					const int& nstep,
					const int& mu,
					const int& min_staple_dir=0)
  {
    ape_smear_conf(smear_conf,origi_conf,alpha,nstep,only_dir[mu],min_staple_dir);
  }
  
  inline void ape_perp_dir_smear_conf(LxField<quad_su3>& smear_conf,
				      const LxField<quad_su3>& origi_conf,
				      const double& alpha,
				      const int& nstep,
				      const int& mu,
				      const int& min_staple_dir=0)
  {
    ape_smear_conf(smear_conf,origi_conf,alpha,nstep,all_other_dirs[mu]);
  }
  
  inline void ape_temporal_smear_conf(LxField<quad_su3>& smear_conf,
				      const LxField<quad_su3>& origi_conf,
				      const double& alpha,
				      const int& nstep)
  {
    ape_single_dir_smear_conf(smear_conf,origi_conf,alpha,nstep,0);
  }
  
  inline void ape_spatial_smear_conf(LxField<quad_su3>& smear_conf,
				     const LxField<quad_su3>& origi_conf,
				     const double& alpha,
				     const int& nstep)
  {
    ape_perp_dir_smear_conf(smear_conf,origi_conf,alpha,nstep,0,1);
  }
}

#endif
