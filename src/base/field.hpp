#ifndef _FIELD_HPP
#define _FIELD_HPP

#include <type_traits>
#include <cstdint>
#include <new_types/su3.hpp>

namespace nissa
{
  template <typename T>
  struct field
  {
    typedef std::remove_all_extents_t<T> Base;
    
    static constexpr int nstatic_comps=sizeof(T)/sizeof(Base);
    
    Base *data;
    
    int64_t l;
    
    bool border_is_valid;
    
    /// Calculate the index - no more components to parse
    template <typename I>
    int64_t indexInternal(const int64_t outer,
			  const I* dummy) const
    {
      return outer;
    }
    
    template <typename I,
	      typename...Tail>
    int64_t indexInternal(const int64_t outer,
			  const I* dummy,
			  const int64_t thisComp,
			  Tail&&...innerComps) const
    {
      using Inner=std::remove_reference_t<decltype(**dummy)>;
      
      const int64_t thisSize=sizeof(I)/sizeof(Inner);
      
      const int64_t thisVal=outer*thisSize+thisComp;
      
      return indexInternal(thisVal,(Inner*)nullptr,std::forward<Tail>(innerComps)...);
    }
    
    template <typename...I>
    int64_t index(I&&...i) const
    {
      return indexInternal(0,(T*)nullptr,std::forward<I>(i)...);
    }
    
    template <typename...I>
    Base& operator()(I&&...i)
    {
      return data[index(std::forward<I>(i)...)];
    }
  };
  
  void new_index_test();
}

#endif
