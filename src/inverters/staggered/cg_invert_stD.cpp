#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "cg_invert_evn_stD.hpp"

#include "base/old_field.hpp"
// #include "base/quda_bridge.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"

//this is the famous trick to invert the full D matrix using e/o precond: sol[ODD]=1/m*(source[ODD]-Doe*sol[EVN])

namespace nissa
{
  void inv_stD_cg(EoField<color>& sol,
		  const std::optional<EvnField<color>>& guess,
		  const EoField<quad_su3>& conf,
		  const double& m,
		  const int& niter,
		  const double& residue,
		  const EoField<color>& source)
  {
    // bool solved=false;
    
    // if(checkIfQudaAvailableAndRequired() and not solved)
    //   {
    // 	source.rese
    // 	for(int eo=0;eo<2;eo++)
    // 	  vector_reset(source[eo]);
    // 	if(is_master_rank())
    // 	  source[EVN][0][0][0]=1.0;
	
    // 	//test
    // 	eo_ptr<color> test={nissa_malloc("temp_evn",locVolh,color),nissa_malloc("temp_odd",locVolh,color)};
    // 	apply_stD(test,conf,m,source);
	
    // 	solved=quda_iface::solve_stD(sol,conf,m,niter,residue,source);
    // 	for(int eo=0;eo<2;eo++)
    // 	  for(int ivol=0;ivol<locVolh;ivol++)
    // 	    for(int ic=0;ic<NCOL;ic++)
    // 	      for(int ri=0;ri<2;ri++)
    // 		{
    // 		  const double& e=sol[eo][ivol][ic][ri];
    // 		  const double& f=test[eo][ivol][ic][ri];
    // 		  if(fabs(e) or fabs(f))
    // 		    printf("%d %d %d %d %lg %lg\n",eo,ivol,ic,ri,e,f);
    // 		}
    // 	for(int eo=0;eo<2;eo++)
    // 	  nissa_free(test[eo]);
    //   }
    
    // if(not solved)
    //   {
    inv_evn_stD_cg(sol.evenPart,guess,conf,m,niter,residue,source);
    apply_st2Doe(sol.oddPart,conf,sol.evenPart);
    crash("reimplement");
    // double_vector_linear_comb((double*)(sol[ODD].data),(double*)(source[ODD].data),1/m,(double*)(sol[ODD].data),-0.5/m,locVolh*6);
      // }
    
    // //check solution
    // eo_ptr<color> residueVec={nissa_malloc("temp_evn",locVolh,color),nissa_malloc("temp_odd",locVolh,color)};
    // apply_stD(residueVec,conf,m,sol);
    
    // /// Source L2 norm
    // double sourceNorm2=0.0;
    
    // /// Residue L2 norm
    // double residueNorm2=0.0;
    
    // for(int eo=0;eo<2;eo++)
    //   {
    // 	double_vector_subtassign((double*)residueVec[eo],(double*)source[eo],locVolh*sizeof(color)/sizeof(double));
	
    // 	sourceNorm2+=double_vector_glb_norm2(source[eo],locVolh);
    //     residueNorm2+=double_vector_glb_norm2(residueVec[eo],locVolh);
    //   }
    
    // master_printf("check solution, residue: %lg/%lg=%lg, target one: %lg\n",residueNorm2,sourceNorm2,residueNorm2/sourceNorm2,residue);
    
    // for(int eo=0;eo<2;eo++)
    //   nissa_free(residueVec[eo]);
  }
}
