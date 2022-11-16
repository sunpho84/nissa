#ifndef _STOUT_HPP
#define _STOUT_HPP

#include <sstream>

#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "routines/ios.hpp"
#include "new_types/su3_op.hpp"

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
    
    int def_nlevels()
    {
       return 0;
    }
    double def_rho(){return 0;}
    
    int master_fprintf(FILE *fout,int full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(int full=false)
    {
      std::ostringstream os;
      if(full or is_nonstandard())
	{
	  os<<"Stout\n";
	  if(full or nlevels!=def_nlevels()) os<<" NLevels\t=\t"<<nlevels<<"\n";
	  if(full or rho!=def_rho()) os<<" Rho\t\t=\t"<<rho<<"\n";
	}
      return os.str();
    }
    
    int is_nonstandard()
    {
      return
	nlevels!=def_nlevels() or
	rho!=def_rho();
    }
    
    stout_pars_t() :
      nlevels(def_nlevels()),
      rho(def_rho()) {}
  };
  
  CUDA_HOST_AND_DEVICE void stouted_force_compute_Lambda(su3 Lambda,su3 U,su3 F,hermitian_exp_ingredients *ing);
  //eo
  void stout_smear_whole_stack(eo_ptr<quad_su3> *out,eo_ptr<quad_su3> in,const stout_pars_t& stout_pars,const Coords<bool>& dirs=all_dirs);
  void stout_smear(eo_ptr<quad_su3> ext_out,eo_ptr<quad_su3> ext_in,const stout_pars_t& stout_pars,const Coords<bool>& dirs=all_dirs);
  void stout_smear_single_level(eo_ptr<quad_su3> out,eo_ptr<quad_su3> ext_in,const double& rho,const Coords<bool>& dirs=all_dirs);
  CUDA_HOST_AND_DEVICE void stout_smear_compute_staples(stout_link_staples *out,eo_ptr<quad_su3> conf,int p,const LocLxSite& A,const Dir& mu,const double& rho);
  CUDA_HOST_AND_DEVICE void stout_smear_compute_weighted_staples(su3 staples,eo_ptr<quad_su3> conf,int p,const LocLxSite& A,const Dir& mu,const double& rho);
  void stout_smear_conf_stack_allocate(eo_ptr<quad_su3> **out,eo_ptr<quad_su3> in,int nlev);
  void stout_smear_conf_stack_free(eo_ptr<quad_su3> **out,int nlev);
  void stouted_force_compute_Lambda(su3 Lambda,su3 U,su3 F,hermitian_exp_ingredients *ing);
  void stouted_force_remap(eo_ptr<quad_su3> F,eo_ptr<quad_su3> *sme_conf,const stout_pars_t& stout_pars);
  void stouted_force_remap_step(eo_ptr<quad_su3> *F,eo_ptr<quad_su3> *conf,const double& rho);
  //lx
  void stout_smear_whole_stack(quad_su3 **out,quad_su3 *in,const stout_pars_t& stout_pars,const Coords<bool>& dirs=all_dirs);
  void stout_smear(quad_su3 *ext_out,quad_su3 *ext_in,const stout_pars_t& stout_pars,const Coords<bool>& dirs=all_dirs);
  void stout_smear_single_level(quad_su3 *out,quad_su3 *ext_in,const double& rho,const Coords<bool>& dirs=all_dirs);
  CUDA_HOST_AND_DEVICE void stout_smear_compute_staples(stout_link_staples *out,quad_su3 *conf,int p,const LocLxSite& A,const Dir& mu,const double& rho);
  CUDA_HOST_AND_DEVICE void stout_smear_compute_weighted_staples(su3 staples,quad_su3 *conf,int p,const LocLxSite& A,const Dir& mu,const double& rho);
  void stout_smear_conf_stack_allocate(quad_su3 ***out,quad_su3 *in,int nlev);
  void stout_smear_conf_stack_free(quad_su3 ***out,int nlev);
  void stouted_force_compute_Lambda(su3 Lambda,su3 U,su3 F,hermitian_exp_ingredients *ing);
  void stouted_force_remap(quad_su3 *F,quad_su3 **sme_conf,const stout_pars_t& stout_pars);
  void stouted_force_remap_step(quad_su3 *F,quad_su3 *conf,const double& rho);
}

#endif
