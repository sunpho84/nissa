#ifndef _HYP_HPP
#define _HYP_HPP

#include "routines/ios.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  //parameters to perform hyp
  struct hyp_pars_t
  {
    double alpha0;
    double alpha1;
    double alpha2;
    int nlevels;
    
    double def_alpha0() {return 1;}
    double def_alpha1() {return 1;}
    double def_alpha2() {return 0.5;}
    int def_nlevels() {return 1;}
    
    int master_fprintf(FILE *fout,bool full=false)
    {
      int nprinted=0;
      
      if(full||is_nonstandard())
	{
	  nprinted+=nissa::master_fprintf(fout,"HypPars\n");
	  if(full||alpha0!=def_alpha0()||alpha1!=def_alpha1()||alpha2!=def_alpha2()) nprinted+=nissa::master_fprintf(fout," Alphas\t\t=\t{%lg,%lg,%lg}\n",alpha0,alpha1,alpha2);
	  if(full||nlevels!=def_nlevels()) nprinted+=nissa::master_fprintf(fout," NLevels\t=\t%d\n",nlevels);
	}
      
      return nprinted;
    }
    
    int is_nonstandard()
    {
      return
	alpha0!=def_alpha0()||
	alpha1!=def_alpha1()||
	alpha2!=def_alpha2()||
	nlevels!=def_nlevels();
    }
    hyp_pars_t() :
      alpha0(def_alpha0()),
      alpha1(def_alpha1()),
      alpha2(def_alpha2()),
      nlevels(def_nlevels()) {}
  };
  
  void hyp_smear_conf(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2);
  void hyp_smear_conf_dir(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2,int req_mu);
}

#endif
