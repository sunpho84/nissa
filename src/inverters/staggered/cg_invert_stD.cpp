#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "cg_invert_evn_stD.hpp"

#include "base/quda_bridge.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"

//this is the famous trick to invert the full D matrix using e/o precond: sol[ODD]=1/m*(source[ODD]-Doe*sol[EVN])

namespace nissa
{
  void inv_stD_cg(eo_ptr<color> sol,color *guess,eo_ptr<quad_su3> conf,double m,int niter,double residue,eo_ptr<color> source)
  {
    bool solved=false;
    
    if(checkIfQudaAvailableAndRequired() and not solved)
      solved=quda_iface::solve_stD(sol,conf,m,niter,residue,source);
    
    if(not solved)
      {
	inv_evn_stD_cg(sol[EVN],guess,conf,m,niter,residue,source);
	apply_st2Doe(sol[ODD],conf,sol[EVN]);
	double_vector_linear_comb((double*)(sol[ODD]),(double*)(source[ODD]),1/m,(double*)(sol[ODD]),-0.5/m,locVolh*6);
      }
    
    //check solution
    eo_ptr<color> residueVec={nissa_malloc("temp_evn",locVolh,color),nissa_malloc("temp_odd",locVolh,color)};
    apply_stD(residueVec,conf,m,sol);
    
    /// Source L2 norm
    double sourceNorm2=0.0;
    
    /// Residue L2 norm
    double residueNorm2=0.0;
    
    for(int eo=0;eo<2;eo++)
      {
	double_vector_subtassign((double*)residueVec[eo],(double*)source[eo],locVolh*sizeof(color)/sizeof(double));
	
	sourceNorm2+=double_vector_glb_norm2(source[eo],locVolh);
        residueNorm2+=double_vector_glb_norm2(residueVec[eo],locVol);
      }
    
    master_printf("check solution, residue: %lg/%lg=%lg, target one: %lg\n",residueNorm2,sourceNorm2,residueNorm2/sourceNorm2,residue);
    
    for(int eo=0;eo<2;eo++)
      nissa_free(residueVec[eo]);
  }
}
