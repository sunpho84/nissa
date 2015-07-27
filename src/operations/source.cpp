#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>

#include "new_types/new_types_definitions.hpp"

#include "base/global_variables.hpp"
#include "base/vectors.hpp"

namespace nissa
{
  //select a timeslice
  template<class prop_type> void select_propagator_timeslice_internal(prop_type *prop_out,prop_type *prop_in,int timeslice)
  {
    if(prop_out!=prop_in) vector_copy(prop_out,prop_in);
    
    //put to zero everywhere but on the slice
    NISSA_LOC_VOL_LOOP(ivol)
      if(timeslice>=0 && timeslice<glb_size[0] && glb_coord_of_loclx[ivol][0]!=timeslice)
	memset(prop_out[ivol],0,sizeof(prop_type));
  }
  
  void select_propagator_timeslice(colorspinspin *prop_out,colorspinspin *prop_in,int timeslice)
  {select_propagator_timeslice_internal(prop_out,prop_in,timeslice);}
  void select_propagator_timeslice(su3spinspin *prop_out,su3spinspin *prop_in,int timeslice)
  {select_propagator_timeslice_internal(prop_out,prop_in,timeslice);}
}
