#ifndef _QUDA_BRIDGE_HPP
#define _QUDA_BRIDGE_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_QUDA
# include <quda.h>
#endif

#include <complex>
#include <map>

#include "base/multiGridParams.hpp"
#include "io/checksum.hpp"
#include "routines/ios.hpp"
#include "geometry/geometry_eo.hpp"

#ifndef EXTERN_QUDA_BRIDGE
 #define EXTERN_QUDA_BRIDGE extern
 #define INIT_QUDA_BRIDGE_TO(cond)
#else
 #define INIT_QUDA_BRIDGE_TO(cond) cond
#endif

#ifdef USE_QUDA

namespace nissa
{
  // namespace Robbery
  // {
  //   enum ROB_MG{param_coarse,coarse,coarse_solver,solver,evecs,evals,diracCoarseSmoother,Y_d,Yhat_d,gauge};
    
  //   /// Allow to rob the param_coarse
  //   template struct Rob<param_coarse,quda::MG,&quda::MG::param_coarse>;
    
  //   /// Allow to rob the coarse
  //   template struct Rob<coarse,quda::MG,&quda::MG::coarse>;
    
  //   /// Allow to rob the coarse_solver
  //   template struct Rob<coarse_solver,quda::MG,&quda::MG::coarse_solver>;
    
  //   /// Allow to rob the solver
  //   template struct Rob<solver,quda::PreconditionedSolver,&quda::PreconditionedSolver::solver>;
    
  //   /// Allow to rob the evecs
  //   template struct Rob<evecs,quda::Solver,&quda::Solver::evecs>;
    
  //   /// Allow to rob the evals
  //   template struct Rob<evals,quda::Solver,&quda::Solver::evals>;
    
  //   /// Allow to rob the diracCoarseSmoother of a MG
  //   template struct Rob<diracCoarseSmoother,quda::MG,&quda::MG::diracCoarseSmoother>;
    
  //   /// Allow to rob the Y_d of a DiracCoarse
  //   template struct Rob<Y_d,quda::DiracCoarse,&quda::DiracCoarse::Y_d>;
    
  //   /// Allow to rob the Yhat_h of a DiracCoarse
  //   template struct Rob<Yhat_d,quda::DiracCoarse,&quda::DiracCoarse::Yhat_d>;
  // }
}

#endif

namespace quda_iface
{
  /// Tags needed to define a setup
  using SetupID=
    std::tuple<std::string,nissa::Checksum>;

#ifdef USE_QUDA
  // /// Restore or take copy of raw data, taking care of the direction of the request
  // inline void restoreOrTakeCopyOfData(void* host,
  // 				      void* device,
  // 				      const size_t& size,
  // 				      const bool takeCopy)
  // {
  //   qudaMemcpy(takeCopy?host:device,
  // 	       takeCopy?device:host,
  // 	       size,takeCopy?qudaMemcpyDeviceToHost:qudaMemcpyHostToDevice);
  // }
#endif
  
//   /// Store a full setup of the multigrid
//   struct QudaSetup
//   {
//     /// Allocated memory on B and eVecs
//     size_t allocatedMemory;
    
//     /// B vectors needed to restrain/prolongate
//     std::vector<std::vector<char*>> B;
    
//     /// Y vectors needed to implement operator
//     std::vector<char*> Y;
    
//     /// Yhat vectors needed to implement operator
//     std::vector<char*> Yhat;
    
//     /// Eigenvectors
//     std::vector<char*> eVecs;
    
//     /// Eigenvalues
//     std::vector<std::complex<double>> eVals;
    
// #ifdef USE_QUDA
    
//     /// Implants this setup into Quda
//     void restore()
//     {
//       restoreOrTakeCopy(false);
//     }
    
//     /// Explant the current Quda setup into this setup
//     void takeCopy()
//     {
//       restoreOrTakeCopy(true);
//     }
    
//     /// Restore or take copy of the B vectors for a given level
//     void restoreOrTakeCopyOfB(const bool takeCopy,
// 			      std::vector<quda::ColorSpinorField*>& Bdev,
// 			      const size_t lev);
    
//     /// Restore or take copy of the eigenvectors
//     void restoreOrTakeCopyOfEig(const bool takeCopy,
// 				quda::Solver* csv);
    
//     /// Restore or take copy of the Y and Yhat
//     void restoreOrTakeCopyOfAllY(const bool takeCopy);
    
//     /// Unified method to take copy or restore
//     void restoreOrTakeCopy(const bool takeCopy=false);
    
//     /// Reset the setup
//     void reset()
//     {
//       master_printf("Resetting stored setup\n");
//       for(auto& Bl : B)
// 	for(auto& Bli : Bl)
// 	nissa::nissa_free(Bli);
//       B.clear();
      
//       for(auto& ei : eVecs)
// 	nissa::nissa_free(ei);
//       eVecs.clear();
      
//       eVals.clear();
      
//       allocatedMemory=0;
//     }
    
//     /// Destructor
//     ~QudaSetup()
//     {
//       reset();
//     }
    
// #endif
    
//     QudaSetup()=default;
    
//     QudaSetup(const QudaSetup&)=default;
    
//     QudaSetup(QudaSetup&&)=default;
//   };
  
//   EXTERN_QUDA_BRIDGE std::map<SetupID,QudaSetup> qudaSetups;
}

namespace quda_iface
{
  using su3_ptr=nissa::su3*;
  using quda_conf_t=nissa::my_array<su3_ptr,NDIM>;
  
#ifdef USE_QUDA
  EXTERN_QUDA_BRIDGE QudaGaugeParam  gauge_param;
  EXTERN_QUDA_BRIDGE QudaInvertParam inv_param;
  
  EXTERN_QUDA_BRIDGE QudaMultigridParam quda_mg_param;
  EXTERN_QUDA_BRIDGE QudaInvertParam inv_mg_param;
  EXTERN_QUDA_BRIDGE QudaEigParam mg_eig_param[QUDA_MAX_MG_LEVEL];
  
  EXTERN_QUDA_BRIDGE MPI_Comm cart_comm;
#endif
  
  EXTERN_QUDA_BRIDGE void* quda_mg_preconditioner INIT_QUDA_BRIDGE_TO(=nullptr);
  
  EXTERN_QUDA_BRIDGE bool inited INIT_QUDA_BRIDGE_TO(=false);
  
  /// Conf used to remap
  EXTERN_QUDA_BRIDGE quda_conf_t quda_conf INIT_QUDA_BRIDGE_TO({});
}

#include "new_types/su3.hpp"

namespace nissa
{
  EXTERN_QUDA_BRIDGE int use_quda INIT_QUDA_BRIDGE_TO(=true);
  
  /// If Quda is available, check if requested
  inline bool checkIfQudaAvailableAndRequired()
  {
#ifdef USE_QUDA
    if(use_quda)
      return true;
    else
#endif
      return false;
  }
}

#ifdef USE_QUDA
# define QUDA_API
# define QUDA_ESCAPE_IF_NOT_AVAILABLE(ARGS...)
#else
# define QUDA_API inline
# define QUDA_ESCAPE_IF_NOT_AVAILABLE(ARGS...) {crash("Quda not available!"); ARGS}
#endif

#ifndef USE_QUDA
  QUDA_API void freeGaugeQuda() QUDA_ESCAPE_IF_NOT_AVAILABLE();
  QUDA_API void plaqQuda(double*) QUDA_ESCAPE_IF_NOT_AVAILABLE();
  QUDA_API void loadGaugeQuda(void*,void*) QUDA_ESCAPE_IF_NOT_AVAILABLE();
#endif

namespace quda_iface
{
  using namespace nissa;

  QUDA_API void initialize() QUDA_ESCAPE_IF_NOT_AVAILABLE();
  QUDA_API void finalize() QUDA_ESCAPE_IF_NOT_AVAILABLE();
  QUDA_API void apply_tmD(spincolor *out,quad_su3 *conf,double kappa,double csw,double mu,spincolor *in) QUDA_ESCAPE_IF_NOT_AVAILABLE();
  
  QUDA_API bool solve_tmD(LxField<spincolor>& sol,
			  const LxField<quad_su3>& conf,
			  const double& kappa,
			  const double& csw,
			  const double& mu,
			  const int& niter,
			  const double& residue,
			  const LxField<spincolor>& source)
    QUDA_ESCAPE_IF_NOT_AVAILABLE(return 0;);
  
  QUDA_API bool solve_stD(eo_ptr<color> sol,eo_ptr<quad_su3> conf,const double& mass,const int& niter,const double& residue,eo_ptr<color> source) QUDA_ESCAPE_IF_NOT_AVAILABLE(return 0;);
  
  void remap_nissa_to_quda(quda_conf_t& out,
			   const LxField<quad_su3>& in);
  
  /// Load a gauge conf
  template <typename T>
  double load_conf(const T& nissa_conf)
  {
    master_printf("freeing the QUDA gauge conf\n");
    freeGaugeQuda();
    
    remap_nissa_to_quda(quda_conf,nissa_conf);
    master_printf("loading to QUDA the gauge conf\n");
#ifdef USE_QUDA
    loadGaugeQuda((void*)&quda_conf[0],&gauge_param);
#endif
    
    double plaq[3]={}; //to avoid warning when quda not available
    plaqQuda(plaq);
    
    return plaq[0];
  }
}

#undef QUDA_API
#undef QUDA_ESCAPE_IF_NOT_AVAILABLE

#undef INIT_QUDA_BRIDGE_TO
#undef EXTERN_QUDA_BRIDGE

#endif
