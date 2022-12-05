#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <base/field.hpp>
#include <communicate/borders.hpp>
#include <free_theory/free_theory_types_routines.hpp>

namespace nissa
{
  gauge_info create_tlSym_gauge_info(const double& alpha,
				     const momentum_t& bc,
				     const double c1)
  {
    gauge_info out;
    
    out.bc=bc;
    out.alpha=alpha;
    out.c1=c1;
    
    return out;
  }
  
  gauge_info create_Wilson_gauge_info(const double& alpha,
				      const momentum_t& bc)
  {
    return create_tlSym_gauge_info(alpha,bc,0);
  }
  
  tm_quark_info create_twisted_quark_info(const double& kappa,
					  const double& mass,
					  const momentum_t& bc,
					  const int& r,
					  const double& zmp)
  {
    tm_quark_info out;
    
    out.bc=bc;
    out.kappa=kappa;
    out.zmp=zmp;
    out.mass=mass;
    out.r=r;
    
    return out;
  }
  
  tm_quark_info create_Wilson_quark_info(const double& kappa,
					 const momentum_t& bc)
  {
    return create_twisted_quark_info(kappa,0,bc,0);
  }
  
  void get_spin_from_spinspin(LxField<spin>& out,
			      const LxField<spinspin>& in,
			      const int& id_so)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int id_si=0;id_si<NDIRAC;id_si++)
	complex_copy(out[ivol][id_si],in[ivol][id_si][id_so]);
    NISSA_PARALLEL_LOOP_END;
    
    out.invalidateHalo();
  }
  
  void put_spin_into_spinspin(LxField<spinspin>& out,
			      const LxField<spin>& in,
			      const int& id_so)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int id_si=0;id_si<NDIRAC;id_si++)
	complex_copy(out[ivol][id_si][id_so],in[ivol][id_si]);
    NISSA_PARALLEL_LOOP_END;
    
    out.invalidateHalo();
  }
}
