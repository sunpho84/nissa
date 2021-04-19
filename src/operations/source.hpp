#ifndef _SOURCE_HPP
#define _SOURCE_HPP

#include "geometry/geometry_lx.hpp"
#include "new_types/coords.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //select a timeslice
  template<class prop_type> void select_propagator_timeslice(prop_type *prop_out,prop_type *prop_in,const GlbCoord timeslice)
  {
    if(prop_out!=prop_in) vector_copy(prop_out,prop_in);
    
    //put to zero everywhere but on the slice
    if(timeslice>=0 and timeslice<glbTimeSize)
      {
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  if(glbCoordOfLoclx(ivol,timeDirection)!=timeslice)
	    memset(prop_out[ivol.nastyConvert()],0,sizeof(prop_type));
	NISSA_PARALLEL_LOOP_END;
	set_borders_invalid(prop_out);
      }
  }
}

#endif
