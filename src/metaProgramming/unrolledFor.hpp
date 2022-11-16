#ifndef _UNROLLED_FOR_HPP
#define _UNROLLED_FOR_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

/// \file unrolledFor.hpp
///
/// \brief Provides functions to unroll loops

#include <tuple>
#include <utility>

#include <base/cuda.hpp>
#include <metaProgramming/inliner.hpp>

namespace nissa
{
  namespace resources
  {
    /// Wraps the function to be called
    ///
    /// Return an integer, to allow variadic expansion of the
    /// unrolling without any recursion
    template <typename F,
	      typename...Args>
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE constexpr
    int call(const F& f,          ///< Function to be called
	     Args&&...args)       ///< Calling arguments
    {
      f(std::forward<Args>(args)...);
      
      return 0;
    }
    
    /// Unroll a loop
    ///
    /// Actual implementation
    template <int offset,
	      int...Is,
	      typename F>
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE constexpr
    void unrollForInternal(std::integer_sequence<int,Is...>,
			     const F& f)
    {
      /// Dummy initialized list, discarded at compile time
      ///
      /// The attribute avoids compiler warning.
      [[ maybe_unused ]]
      auto list=
	{call(f,Is+offset)...,0};
    }
  }
  
  /// Unroll a loop, wrapping the actual implementation
  template <int min,
	    int max,
	    typename F>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE constexpr
  void unrollFor(const F& f)
  {
    resources::unrollForInternal<min>(std::make_integer_sequence<int,max-min>{},f);
  }
  
  /// Unroll a loop in the range [min,max)
  ///
  /// The called function must take as a first argument the iteration
  template <int min,                               // Minimal value
	    int max,                               // Maximal (excluded) value
	    typename F,                            // Type of the function
	    typename...Args>                       // Types of the arguments to pass
  INLINE_FUNCTION void unrollFor2(F f,             ///< Function to cal
				    Args&&...args)   ///< Arguments to pass
  {
    /// Current value
    constexpr int cur=
      min;
    
    /// Next value
    constexpr int next=
      cur+1;
    
    // Exec the function
    f(cur,std::forward<Args>(args)...);
    
    // Iterates on the loop
    if constexpr(next<max)
      unrollFor2<next,max,F,Args...>(f,std::forward<Args>(args)...);
  }
  
  
#ifdef COMPILING_FOR_DEVICE
  
  /// Uses nvcc builtin
  ///
  /// \todo move it to a dedicated macro to call the proper bultin
# define UNROLL_FOR(TYPE,NAME,MIN,MAX)			\
  PRAGMA(unroll (MAX-MIN))				\
  for(TYPE NAME=MIN;I<MAX;I++)				\
    {
  
# define UNROLL_FOR_END			\
  }
  
#else
  
  /// Create an unrolled for
  ///
  /// Hides the complexity
# define UNROLL_FOR(TYPE,NAME,MIN,MAX)			\
  unrollFor2<MIN,MAX>([&](const TYPE& NAME) INLINE_ATTRIBUTE {
  
  /// Finish an unrolled for
# define UNROLL_FOR_END })
  
#endif
}

#endif
