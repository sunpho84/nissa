#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "cg_invert_evn_stD.hpp"

#ifdef USE_QUDA
# include "base/quda_bridge.hpp"
#endif

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
    enum SolverType{NATIVE_SOLVER
#ifdef USE_QUDA
      ,QUDA_SOLVER
#endif
    };
    
    /// Take note of whther we want to use an external solver
    SolverType solverType=NATIVE_SOLVER;
    
#ifdef USE_QUDA
    if(use_quda)
	solverType=QUDA_SOLVER;
#endif
    
    switch(solverType)
      {
	// Quda
#ifdef USE_QUDA
      case QUDA_SOLVER:
	master_printf("Using QUDA here, %s\n",__FUNCTION__);
	quda_iface::solve_stD(sol,conf,m,niter,residue,source);
	break;
#endif
	
      case NATIVE_SOLVER:
	inv_evn_stD_cg(sol[EVN],guess,conf,m,niter,residue,source);
	apply_st2Doe(sol[ODD],conf,sol[EVN]);
	double_vector_linear_comb((double*)(sol[ODD]),(double*)(source[ODD]),1/m,(double*)(sol[ODD]),-0.5/m,locVolh*6);
	break;
      }
    
  }
}
