#ifndef _SOURCE_HPP
#define _SOURCE_HPP

#include <base/field.hpp>
#include "threads/threads.hpp"

namespace nissa
{
  //select a timeslice
  template <typename T>
  void select_propagator_timeslice(LxField<T>& prop_out,
				   const LxField<T>& prop_in,
				   const int& timeslice)
  {
    if(prop_out!=prop_in)
      prop_out=prop_in;
    
    //put to zero everywhere but on the slice
    if(timeslice>=0 and timeslice<glbSize[0])
      {
	prop_out.forEachSiteDeg([&timeslice](double& p,
				   const int& ivol,
				   const int& deg)
	{
	  if(glbCoordOfLoclx[ivol][0]!=timeslice)
	    p=0;
	});
	
	set_borders_invalid(prop_out);
      }
  }
}

#endif
