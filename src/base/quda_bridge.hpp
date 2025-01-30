#ifndef _QUDA_BRIDGE_HPP
#define _QUDA_BRIDGE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#ifdef USE_QUDA
# include <quda.h>
# ifdef READY_TO_DEL
#  include <dirac_quda.h>
#  include <invert_quda.h>
#  include <color_spinor_field.h>
# endif
#endif

#include <complex>
#include <map>

#include "base/multiGridParams.hpp"
#include "io/checksum.hpp"
#include "routines/ios.hpp"
#include "geometry/geometry_eo.hpp"

#ifndef EXTERN_QUDA_BRIDGE
# define EXTERN_QUDA_BRIDGE extern
# define INIT_QUDA_BRIDGE_TO(cond)
#else
# define INIT_QUDA_BRIDGE_TO(cond) cond
#endif

namespace quda_iface
{
  using su3_ptr=nissa::su3*;
  using quda_conf_t=nissa::MyArray<su3_ptr,NDIM>;
  
#ifdef USE_QUDA
  EXTERN_QUDA_BRIDGE QudaGaugeParam  gauge_param;
  EXTERN_QUDA_BRIDGE QudaInvertParam inv_param;
  
  EXTERN_QUDA_BRIDGE QudaMultigridParam quda_mg_param;
  EXTERN_QUDA_BRIDGE QudaInvertParam inv_mg_param;
  EXTERN_QUDA_BRIDGE QudaEigParam mg_eig_param[QUDA_MAX_MG_LEVEL];
  
  EXTERN_QUDA_BRIDGE MPI_Comm cart_comm;
  
# ifdef READY_TO_DEL
  EXTERN_QUDA_BRIDGE quda::Dirac* D;
  EXTERN_QUDA_BRIDGE quda::Dirac* DSloppy;
  EXTERN_QUDA_BRIDGE quda::Dirac* DPre;
  EXTERN_QUDA_BRIDGE quda::DiracMatrix* M;
  EXTERN_QUDA_BRIDGE quda::DiracMatrix* MSloppy;
  EXTERN_QUDA_BRIDGE quda::DiracMatrix* MPre;
  EXTERN_QUDA_BRIDGE quda::SolverParam* solverParam;
  EXTERN_QUDA_BRIDGE quda::Solver* solver;
  EXTERN_QUDA_BRIDGE quda::ColorSpinorField* b;
  EXTERN_QUDA_BRIDGE quda::ColorSpinorField* x;
  EXTERN_QUDA_BRIDGE std::map<std::string,quda::TimeProfile> profilers;
# endif
  
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
    MASTER_PRINTF("freeing the QUDA gauge conf\n");
    freeGaugeQuda();
    
    remap_nissa_to_quda(quda_conf,nissa_conf);
    MASTER_PRINTF("loading to QUDA the gauge conf\n");
#ifdef USE_QUDA
    loadGaugeQuda((void*)&quda_conf[0],&gauge_param);
#endif
    
    double plaq[3]={}; //to avoid warning when quda not available
    plaqQuda(plaq);
    
    return plaq[0];
  }
  
  /// Since the solver might have not deleted the eigenvectors, try to flag them so maybe they will be deleted
  void maybeFlagTheMultigridEigenVectorsForDeletion();
}

#undef QUDA_API
#undef QUDA_ESCAPE_IF_NOT_AVAILABLE

#undef INIT_QUDA_BRIDGE_TO
#undef EXTERN_QUDA_BRIDGE

#endif
