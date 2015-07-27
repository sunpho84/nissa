#ifndef _SOURCE_H
#define _SOURCE_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  //select a timeslice
  template<class prop_type> void select_propagator_timeslice(prop_type *prop_out,prop_type *prop_in,int timeslice)
  {
    GET_THREAD_ID();
    
    if(prop_out!=prop_in) vector_copy(prop_out,prop_in);
    
    //put to zero everywhere but on the slice
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      if(timeslice>=0 && timeslice<glb_size[0] && glb_coord_of_loclx[ivol][0]!=timeslice)
	memset(prop_out[ivol],0,sizeof(prop_type));
    set_borders_invalid(prop_out);
  }
}

#endif
