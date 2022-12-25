#ifndef _SOURCE_HPP
#define _SOURCE_HPP

#include <base/field.hpp>
#include <threads/threads.hpp>

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
	FOR_EACH_SITE_DEG_OF_FIELD(prop_out,
				   CAPTURE(timeslice,
					   TO_WRITE(prop_out)),
				   ivol,deg,
	{
	  if(glbCoordOfLoclx[ivol][0]!=timeslice)
	    prop_out(ivol,deg)=0;
	});
      }
  }
}

#endif
