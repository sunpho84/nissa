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
		       const TmQuarkInfo& qu,
		       const I& in)
  {
    PAR(0,locVolh,
	CAPTURE(TO_WRITE(out),
		TO_READ(in),
		qu),X,
	{
	  const complex z={1/(2*qu.kappa),qu.mass};
	  
	  for(int id=0;id<NDIRAC/2;id++)
	    unsafe_complex_prod(out[X][id],in[X][id],z);
	  for(int id=NDIRAC/2;id<NDIRAC;id++)
	    unsafe_complex_conj2_prod(out[X][id],in[X][id],z);
	});
  }
  
  //implement Koo defined in equation (7)
  void tmDkern_eoprec_eos(OddField<spin>& out,
			  EvnField<spin>& temp,
			  const TmQuarkInfo& qu,
			  const OddField<spin>& in)
  {
    tmn2Deo_or_tmn2Doe_eos(out,in,qu.bc);
    inv_tmDee_or_oo_eos(temp,qu,out);
    tmn2Deo_or_tmn2Doe_eos(out,temp,qu.bc);
    inv_tmDee_or_oo_eos(temp,qu,out);
    tmDee_or_oo_eos(temp,qu,in);
    
    PAR(0,locVolh,
	CAPTURE(TO_WRITE(out),
		TO_WRITE(temp)),ivol,
	{
	  for(int id=0;id<2;id++)
	    for(int ri=0;ri<2;ri++)
	      { //gamma5 is explicitely implemented
		out[ivol][id  ][ri]=+temp[ivol][id  ][ri]-out[ivol][id  ][ri]*0.25;
		out[ivol][id+2][ri]=-temp[ivol][id+2][ri]+out[ivol][id+2][ri]*0.25;
	      }
	});
  }
  
  //square of Koo
  void tmDkern_eoprec_square_eos(OddField<spin>& out,
				 OddField<spin>& temp1,
				 EvnField<spin> &temp2,
				 const TmQuarkInfo& qu,
				 const OddField<spin>& in)
  {
    TmQuarkInfo mqu=qu;
    mqu.mass*=-1;
    
    tmDkern_eoprec_eos(temp1,temp2,mqu, in   );
    tmDkern_eoprec_eos(out,  temp2,qu,  temp1);
  }
}
