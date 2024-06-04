#ifndef _QUDA_BRIDGE_HPP
#define _QUDA_BRIDGE_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_QUDA
# include <quda.h>
# include <multigrid.h>
#endif

#include <complex>

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
  namespace Robbery
  {
    enum ROB_MG{param_coarse,coarse,coarse_solver,solver,evecs,evals};
    
    /// Allow to rob the param_coarse
    template struct Rob<param_coarse,quda::MG,&quda::MG::param_coarse>;
    
    /// Allow to rob the param_coarse
    template struct Rob<coarse,quda::MG,&quda::MG::coarse>;
    
    /// Allow to rob the coarse_solver
    template struct Rob<coarse_solver,quda::MG,&quda::MG::coarse_solver>;
    
    /// Allow to rob the solver
    template struct Rob<solver,quda::PreconditionedSolver,&quda::PreconditionedSolver::solver>;
    
    /// Allow to rob the evecs
    template struct Rob<evecs,quda::Solver,&quda::Solver::evecs>;
    
    /// Allow to rob the evals
    template struct Rob<evals,quda::Solver,&quda::Solver::evals>;
  }
}
#endif

namespace quda_iface
{
  using SetupID=std::tuple<std::string,nissa::checksum>;
  
  struct QudaSetup
  {
    std::vector<std::vector<char*>> B;
    
    std::vector<char*> eVecs;
    
    std::vector<std::complex<double>> eVals;
    
    /// Implants this setup into Quda
    void restore()
    {
      restoreOrTakeCopy(false);
    }
    
    /// Explant the current Quda setup into this setup
    void takeCopy()
    {
      restoreOrTakeCopy(true);
    }
    
    /// Unified method to take copy or restore
    void restoreOrTakeCopy(const bool takeCopy=false);
    
    /// Reset the setup
    void reset()
    {
      master_printf("Resetting stored setup\n");
      for(auto& Bi : B)
	nissa::nissa_free(Bi);
      B.clear();
      
      for(auto& ei : eVecs)
	nissa::nissa_free(ei);
      eVecs.clear();
      
      eVals.clear();
    }
    
    /// Destructor
    ~QudaSetup()
    {
      reset();
    }
    
    QudaSetup(const QudaSetup&)=default;
    
    QudaSetup(QudaSetup&&)=default;
  };
  
  EXTERN_QUDA_BRIDGE std::map<SetupID,QudaSetup> qudaSetups;
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
  QUDA_API void remap_nissa_to_quda(spincolor *out,spincolor *in) QUDA_ESCAPE_IF_NOT_AVAILABLE();
  QUDA_API void remap_quda_to_nissa(spincolor *out,spincolor *in) QUDA_ESCAPE_IF_NOT_AVAILABLE();
  QUDA_API void remap_nissa_to_quda(quda_conf_t out,quad_su3 *in) QUDA_ESCAPE_IF_NOT_AVAILABLE();
  QUDA_API void remap_nissa_to_quda(quda_conf_t out,eo_ptr<quad_su3> in) QUDA_ESCAPE_IF_NOT_AVAILABLE();
  
  QUDA_API bool solve_tmD(spincolor *sol,quad_su3 *conf,const double& kappa,const double& csw,const double& mu,const int& niter,const double& residue,spincolor *source) QUDA_ESCAPE_IF_NOT_AVAILABLE(return 0;);
  QUDA_API bool solve_stD(eo_ptr<color> sol,eo_ptr<quad_su3> conf,const double& mass,const int& niter,const double& residue,eo_ptr<color> source) QUDA_ESCAPE_IF_NOT_AVAILABLE(return 0;);
  
  /// Load a gauge conf
  template<typename T>
  double load_conf(T nissa_conf)
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
