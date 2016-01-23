#ifndef _STOUT_HPP
#define _STOUT_HPP

#include "routines/ios.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  //structure to store data for stouting
  struct stout_link_staples
  {
    su3 C;
    su3 Omega;
    su3 Q;
  };
  
  //parameters to stout
  struct stout_pars_t
  {
    int nlevels;
    double rho;
    
    int def_nlevels(){return 0;}
    double def_rho(){return 0;}
    
    int master_fprintf(FILE *fout,bool full=false)
    {
      int nprinted=0;
      if(full||is_nonstandard())
	{
	  nprinted+=nissa::master_fprintf(fout,"Stout\n");
	  if(full||nlevels!=def_nlevels()) nprinted+=nissa::master_fprintf(fout," NLevels\t=\t%d\n",nlevels);
	  if(full||rho!=def_rho()) nprinted+=nissa::master_fprintf(fout," Rho\t\t=\t%lg\n",rho);
	}
      
      return nprinted;
    }
    
    int is_nonstandard()
    {
      return
	nlevels!=def_nlevels()||
	rho!=def_rho();
    }
    
    stout_pars_t() :
      nlevels(def_nlevels()),
      rho(def_rho()) {}
  };
  
  void stouted_force_compute_Lambda(su3 Lambda,su3 U,su3 F,hermitian_exp_ingredients *ing);
  //eo
  void stout_smear_whole_stack(quad_su3 ***out,quad_su3 **in,stout_pars_t *stout_pars);
  void stout_smear(quad_su3 **ext_out,quad_su3 **ext_in,stout_pars_t *stout_pars);
  void stout_smear_single_level(quad_su3 **out,quad_su3 **ext_in, double rho);
  void stout_smear_compute_staples(stout_link_staples *out,quad_su3 **conf,int p,int A,int mu,double rho);
  void stout_smear_compute_weighted_staples(su3 staples,quad_su3 **conf,int p,int A,int mu,double rho);
  void stout_smear_conf_stack_allocate(quad_su3 ****out,quad_su3 **in,int nlev);
  void stout_smear_conf_stack_free(quad_su3 ****out,int nlev);
  void stouted_force_compute_Lambda(su3 Lambda,su3 U,su3 F,hermitian_exp_ingredients *ing);
  void stouted_force_remap(quad_su3 **F,quad_su3 ***sme_conf,stout_pars_t *stout_pars);
  void stouted_force_remap_step(quad_su3 **F,quad_su3 **conf,double rho);
  //lx
  void stout_smear_whole_stack(quad_su3 **out,quad_su3 *in,stout_pars_t *stout_pars);
  void stout_smear(quad_su3 *ext_out,quad_su3 *ext_in,stout_pars_t *stout_pars);
  void stout_smear_single_level(quad_su3 *out,quad_su3 *ext_in,double rho);
  void stout_smear_compute_staples(stout_link_staples *out,quad_su3 *conf,int p,int A,int mu,double rho);
  void stout_smear_compute_weighted_staples(su3 staples,quad_su3 *conf,int p,int A,int mu,double rho);
  void stout_smear_conf_stack_allocate(quad_su3 ***out,quad_su3 *in,int nlev);
  void stout_smear_conf_stack_free(quad_su3 ***out,int nlev);
  void stouted_force_compute_Lambda(su3 Lambda,su3 U,su3 F,hermitian_exp_ingredients *ing);
  void stouted_force_remap(quad_su3 *F,quad_su3 **sme_conf,stout_pars_t *stout_pars);
  void stouted_force_remap_step(quad_su3 *F,quad_su3 *conf,double rho);
}

#endif
