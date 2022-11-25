#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <math.h>

#include <base/debug.hpp>
#include <base/vectors.hpp>
#include <communicate/borders.hpp>
#include <free_theory/twisted_free_Dirac_eoprec_operator.hpp>
#include <geometry/geometry_eo.hpp>
#include <new_types/complex.hpp>
#include <threads/threads.hpp>

namespace nissa
{
  // Implement ee or oo part of Dirac operator, equation(3)
  template <typename O,
	    typename I>
  void tmDee_or_oo_eos(O&& out,
		       const tm_quark_info& qu,
		       const I& in)
  {
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      {
	const complex z={1/(2*qu.kappa),qu.mass};
	
	for(int id=0;id<NDIRAC/2;id++)
	  unsafe_complex_prod(out[X][id],in[X][id],z);
	for(int id=NDIRAC/2;id<NDIRAC;id++)
	  unsafe_complex_conj2_prod(out[X][id],in[X][id],z);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //implement Koo defined in equation (7)
  void tmDkern_eoprec_eos(OddField<spin>& out,
			  EvnField<spin>& temp,
			  const tm_quark_info& qu,
			  const OddField<spin>& in)
  {
    tmn2Deo_or_tmn2Doe_eos(out,in,qu.bc);
    inv_tmDee_or_oo_eos(temp,qu,out);
    tmn2Deo_or_tmn2Doe_eos(out,temp,qu.bc);
    inv_tmDee_or_oo_eos(temp,qu,out);
    tmDee_or_oo_eos(temp,qu,in);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVolh)
      for(int id=0;id<2;id++)
	for(int ri=0;ri<2;ri++)
	  { //gamma5 is explicitely implemented
	    out[ivol][id  ][ri]=+temp[ivol][id  ][ri]-out[ivol][id  ][ri]*0.25;
	    out[ivol][id+2][ri]=-temp[ivol][id+2][ri]+out[ivol][id+2][ri]*0.25;
	  }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //square of Koo
  void tmDkern_eoprec_square_eos(OddField<spin>& out,
				 OddField<spin>& temp1,
				 EvnField<spin> &temp2,
				 const tm_quark_info& qu,
				 const OddField<spin>& in)
  {
    tm_quark_info mqu=qu;
    mqu.mass*=-1;
    
    tmDkern_eoprec_eos(temp1,temp2,mqu, in   );
    tmDkern_eoprec_eos(out,  temp2,qu,  temp1);
  }
}
