#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "dirac_operators/momenta/MFACC.hpp"
#include "geometry/geometry_eo.hpp"
#include "inverters/momenta/cg_invert_MFACC.hpp"
#include "hmc/gauge/MFACC_fields.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"

namespace nissa
{
  //compute the action of the momenta
  double momenta_action(eo_ptr<quad_su3> H)
  {
    //summ the square of H
    double glb_action_eo[2];
    for(int eo=0;eo<2;eo++)
      double_vector_glb_scalar_prod(&(glb_action_eo[eo]),(double*)(H[eo]),(double*)(H[eo]),sizeof(quad_su3)/sizeof(double)*loc_volh);
    
    return (glb_action_eo[EVN]+glb_action_eo[ODD])/2;
  }
  
  //lx version
  double momenta_action(quad_su3 *H)
  {
    //summ the square of H
    double glb_action_lx;
    double_vector_glb_scalar_prod(&(glb_action_lx),(double*)(H),(double*)(H),sizeof(quad_su3)/sizeof(double)*loc_vol);
    
    return glb_action_lx/2;
  }
  
  //MFACC accelerated action
  double momenta_action_with_FACC(quad_su3 *conf,double kappa,int niter,double residue,quad_su3 *H)
  {
    quad_su3 *H_temp=nissa_malloc("H_temp",loc_vol+bord_vol,quad_su3);
    inv_MFACC_cg(H_temp,NULL,conf,kappa,niter,residue,H);
    
    //compute according to HGH
    double act;
    double_vector_glb_scalar_prod(&act,(double*)(H),(double*)(H_temp),sizeof(quad_su3)/sizeof(double)*loc_vol);
    nissa_free(H_temp);
    
    return act/2;
  }
  
  //compute the action for the Fourier acceleration-related momenta
  void MFACC_momenta_action(double* tot_action,su3** pi,int naux_fields,quad_su3* conf,double kappa)
  {
    //allocate temporary field where to store output
    su3 *V=nissa_malloc("V",loc_vol,su3);
    
    double glb_action_id[naux_fields];
    for(int id=0;id<naux_fields;id++)
      {
        //apply the kernel
        apply_MFACC(V,conf,kappa,pi[id]);
        double_vector_glb_scalar_prod(&(glb_action_id[id]),(double*)(pi[id]),(double*)V,sizeof(su3)/sizeof(double)*loc_vol);
      }
    
    (*tot_action)=0;
    for(int id=0;id<naux_fields;id++) (*tot_action)+=glb_action_id[id]/2;
    
    nissa_free(V);
  }
}
