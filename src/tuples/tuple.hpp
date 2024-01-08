#ifndef _TUPLE_HPP
#define _TUPLE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <utility>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  namespace impl
  {
    /// Provides the I-th element of the tuple, of type T
    template <std::size_t I,
	      typename T>
    struct _TupleElementProvider
    {
      /// Stored value
      T t;
      
#define PROVIDE_GET(ATTRIB)						\
      /*! Get the component via index*/					\
      INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB			\
      ATTRIB T& _get(std::integral_constant<std::size_t,I>) ATTRIB	\
      {									\
	return t;							\
      }									\
      									\
      /*! Get the component via type */					\
      INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB			\
      ATTRIB T& _get(std::remove_reference_t<T>*) ATTRIB		\
      {									\
	return t;							\
      }
      
      PROVIDE_GET(const);
      
      PROVIDE_GET(/* non const */);
      
#undef PROVIDE_GET
    };
    
    /// Tuple of objects
    ///
    /// Internal implementation, forward declaration
    template <typename I,
	      typename...T>
    struct _Tuple;
    
    /// Tuple of objects
    ///
    /// Internal implementation
   template <std::size_t...I,
	      typename...T>
    struct _Tuple<std::index_sequence<I...>,T...> :
      _TupleElementProvider<I,T>...
    {
      /// Import implementation of get
      using _TupleElementProvider<I,T>::_get...;
      
#define PROVIDE_GET(ATTRIB)						\
      /*! Get the component via index*/					\
      template <std::size_t J>						\
	constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB			\
	decltype(auto) get() ATTRIB					\
      {									\
	return this->_get(std::integral_constant<std::size_t,J>());	\
      }									\
									\
      /*! Get the component via type */					\
      template <typename M>						\
	constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB			\
	decltype(auto) get() ATTRIB					\
      {									\
	return this->_get((std::decay_t<M>*)nullptr);			\
      }
      
      PROVIDE_GET(const);
      
      PROVIDE_GET(/* non const */);
      
#undef PROVIDE_GET
      
      /// Get J-th type of the tuple
      template <std::size_t J>
	using ElementType=
	decltype(tupleGet<J>(std::declval<_Tuple>()));
      
#define PROVIDE_APPLY_TO(ATTRIB)		\
      /*! Pass the tuple arguments to the passed function */	\
      template <typename F>			\
	decltype(auto) applyTo(F&& f) ATTRIB	\
      {						\
	f(_TupleElementProvider<I,T>::t...);	\
      }
      
      PROVIDE_APPLY_TO(const);
      
      PROVIDE_APPLY_TO(/* non const */);
      
#undef PROVIDE_APPLY_TO
    };
  }
  
  /// Tuple of objects
  template <typename...T>
  using Tuple=
    impl::_Tuple<std::make_index_sequence<sizeof...(T)>,T...>;
  
  /// Gets the I-th element of the tuple T
  template <std::size_t I,
	    typename T>
  constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
  decltype(auto) get(T&& t)
  {
    return t[std::integral_constant<std::size_t,I>()];
  }
  
  /// Gets the element of type M of the tuple T
  template <typename M,
	    typename T>
  constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
  decltype(auto) get(T&& t)
  {
    return t[(std::decay_t<M>*)nullptr];
  }
  
  /// Tie the custom tuple
  template <typename...Args>
  constexpr auto tie(Args&... args)
  {
    return Tuple<Args&...>{args...};
  }
}

#endif
