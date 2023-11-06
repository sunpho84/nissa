#ifndef _TYPEPOSITION_HPP
#define _TYPEPOSITION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tuples/typePosition.hpp

#include <cstddef>
#include <tuple>

namespace nissa
{
  namespace impl
  {
    /// Determine the position of a type in a tuple
    ///
    /// Forward declaration
    template <typename T,
	      typename TP>
    struct _TypePositionInTuple;
    
    /// Determine the position of a type in a tuple
    ///
    /// Internal implementation
    template <typename T,
	      typename...TP>
    struct _TypePositionInTuple<T,std::tuple<TP...>>
    {
      /// Computes the position
      static constexpr size_t value()
      {
	///Store types equality
	constexpr bool is[]{std::is_same_v<T,TP>...};
	
	size_t i=0;
	while(i<sizeof...(TP) and is[i]==0)
	  i++;
	
	return i;
      }
    };
  }
  
  /// Determine the position of a type in a tuple
  ///
  ///Tuple size is returned in case no element is present
  template <typename T,
	    typename TP>
  static constexpr size_t typePositionInTuple=
    impl::_TypePositionInTuple<T,TP>::value();
  
  /////////////////////////////////////////////////////////////////
  
  namespace impl
  {
    /// Determine whether elements of tuple tp are consecutive in tuple TP
    ///
    /// Forward declaration
    template <typename tp,
	      typename TP>
    struct _TypesAreConsecutiveInTuple;
    
    /// Determine whether elements of tuple tp are consecutive in tuple TP
    ///
    /// Internal implementation
    template <typename...tp,
	      typename TP>
    struct _TypesAreConsecutiveInTuple<std::tuple<tp...>,
				       TP>
    {
      /// Number of elements of tuple<tp...>
      static constexpr size_t np=
	sizeof...(tp);
      
      /// Number of elements of tuple TP
      static constexpr size_t NP=
	std::tuple_size_v<TP>;
      
      /// Determine consecutivity
      static constexpr bool value()
      {
	if constexpr(np>0)
	  {
	    constexpr size_t pos[]{typePositionInTuple<tp,TP>...};
	    
	    bool ok=(pos[0]!=NP);
	    
	    for(size_t i=1;i<np;i++)
	      ok&=(pos[i]!=NP and pos[i]==pos[i-1]+1);
	    
	    return ok;
	  }
	
	return true;
      }
    };
  }
  
  /// Determine whether elements of tuple tp are consecutive in tuple TP
  template <typename tp,
	    typename TP>
  constexpr bool typesAreConsecutiveInTuple=
    impl::_TypesAreConsecutiveInTuple<tp,TP>::value();
}

#endif
