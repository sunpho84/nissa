#ifndef _APE_HPP
#define _APE_HPP

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
    
    int master_fprintf(FILE *fout,bool full=false)
    {
      int nprinted=0;
      
      if(full||is_nonstandard())
	{
	  nprinted+=nissa::master_fprintf(fout,"Ape\n");
	  if(full||alpha!=def_alpha()) nprinted+=nissa::master_fprintf(fout," Alpha\t\t=\t%lg\n",alpha);
	  if(full||nlevels!=def_nlevels()) nprinted+=nissa::master_fprintf(fout," NLevels\t=\t%d\n",nlevels);
	}
      
      return nprinted;
    }
    
    int is_nonstandard()
    {
      return
	nlevels!=def_nlevels()||
	alpha!=def_alpha();
    }
    
    ape_pars_t() :
      nlevels(def_nlevels()),
      alpha(def_alpha()) {}
  };
  
  void ape_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep,int ndirs=NDIM,int *dir=all_dirs);
  void ape_spatial_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep);
  void ape_temporal_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep);
  inline void ape_single_dir_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep,int mu)
  {ape_smear_conf(smear_conf,origi_conf,alpha,nstep,1,&mu);}
  inline void ape_perp_dir_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep,int mu)
  {ape_smear_conf(smear_conf,origi_conf,alpha,nstep,3,perp_dir[mu]);}
}

#endif
