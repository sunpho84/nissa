#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/field.hpp>
#include <communicate/borders.hpp>
#include <free_theory/free_theory_types_routines.hpp>

namespace nissa
{
  gauge_info create_tlSym_gauge_info(const gauge_info::which_gauge_t& which_gauge,
				     const Momentum& bc,
				     const double c1)
  {
    gauge_info out;
    
    out.bc=bc;
    out.which_gauge=which_gauge;
    out.c1=c1;
    
    return out;
  }
  
  gauge_info create_Wilson_gauge_info(const gauge_info::which_gauge_t& which_gauge,
				      const Momentum& bc)
  {
    return create_tlSym_gauge_info(which_gauge,bc,0);
  }
  
  TmQuarkInfo create_twisted_quark_info(const double& kappa,
					const double& mass,
					const Momentum& bc,
					const int& r,
					const double& zmp)
  {
    TmQuarkInfo out;
    
    out.bc=bc;
    out.kappa=kappa;
    out.zmp=zmp;
    out.mass=mass;
    out.r=r;
    
    return out;
  }
  
  TmQuarkInfo create_Wilson_quark_info(const double& kappa,
				       const Momentum& bc)
  {
    return create_twisted_quark_info(kappa,0,bc,0);
  }
  
  void get_spin_from_spinspin(LxField<spin>& out,
			      const LxField<spinspin>& in,
			      const int& id_so)
  {
    PAR(0,locVol,
	CAPTURE(id_so,
		TO_WRITE(out),
		TO_READ(in)),ivol,
	{
	  for(int id_si=0;id_si<NDIRAC;id_si++)
	    complex_copy(out[ivol][id_si],in[ivol][id_si][id_so]);
	});
  }
  
  void put_spin_into_spinspin(LxField<spinspin>& out,
			      const LxField<spin>& in,
			      const int& id_so)
  {
    PAR(0,locVol,
	CAPTURE(id_so,
		TO_WRITE(out),
		TO_READ(in)),ivol,
	{
	  for(int id_si=0;id_si<NDIRAC;id_si++)
	    complex_copy(out[ivol][id_si][id_so],in[ivol][id_si]);
	});
  }
}
