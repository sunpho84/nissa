#pragma once

#include "base/field.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

#include <dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_portable.hpp>

namespace nissa
{
  //Refers to the doc: "doc/eo_inverter.lyx" for explenations
  
  // template <typename O,
  // 	    typename I>
  //put g5
  void tmDkern_eoprec_eos_put_together_and_include_gamma5(OddField<spincolor>& out,
							  const OddField<spincolor>& temp)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVolh)
      for(int id=0;id<NDIRAC/2;id++)
	for(int ic=0;ic<NCOL;ic++)
	  for(int ri=0;ri<2;ri++)
	    { //gamma5 is explicitely implemented
	      out[ivol][id  ][ic][ri]=+temp[ivol][id  ][ic][ri]-out[ivol][id  ][ic][ri]*0.25;
	      out[ivol][id+NDIRAC/2][ic][ri]=-temp[ivol][id+NDIRAC/2][ic][ri]+out[ivol][id+2][ic][ri]*0.25;
	    }
    NISSA_PARALLEL_LOOP_END;
    
    out.invalidateHalo();
  }
  
  //implement Koo defined in equation (7)
  void tmDkern_eoprec_eos(OddField<spincolor>& out,
			  EvnField<spincolor>& temp,
			  const EoField<quad_su3>& conf,
			  const double& kappa,
			  const double& mu,
			  const OddField<spincolor>& in)
  {
    inv_tmDee_or_oo_eos(temp,kappa,mu,out);
    tmn2Deo_or_tmn2Doe_eos(out,conf,temp);
    tmn2Deo_or_tmn2Doe_eos(out.castSitesCoverage<EVEN_SITES>(),conf,in);
    
    tmDee_or_oo_eos(temp,kappa,mu,in);
    
    tmDkern_eoprec_eos_put_together_and_include_gamma5(out,temp.castSitesCoverage<ODD_SITES>());
  }
  
  //square of Koo
  void tmDkern_eoprec_square_eos(OddField<spincolor>& out,
				 OddField<spincolor>& temp1,
				 EvnField<spincolor>& temp2,
				 const EoField<quad_su3>& conf,
				 const double& kappa,
				 const double& mu,
				 const OddField<spincolor>& in)
  {
    tmDkern_eoprec_eos(temp1,temp2,conf,kappa,-mu, in   );
    tmDkern_eoprec_eos(out,  temp2,conf,kappa,+mu, temp1);
  }
}
