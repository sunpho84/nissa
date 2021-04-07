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
    INLINE_FUNCTION CUDA_HOST_DEVICE
    int call(F&& f,          ///< Function to be called
	     Args&&...args)  ///< Calling arguments
    {
      f(std::forward<Args>(args)...);
      
      return 0;
    }
    
    /// Unroll a loop
    ///
    /// Actual implementation
    template <int...Is,
	      typename F>
    INLINE_FUNCTION CUDA_HOST_DEVICE
    void unrolledForInternal(std::integer_sequence<int,Is...>,F&& f)
    {
      /// Dummy initialized list, discarded at compile time
      ///
      /// The attribute avoids compiler warning.
      [[ maybe_unused ]]
	auto list=
	{call(std::forward<F>(f),Is)...,0};
    }
  }
  
  /// Unroll a loop, wrapping the actual implementation
  template <int N,
	    typename F>
  INLINE_FUNCTION
  void unrolledFor(F&& f)
  {
    resources::unrolledForInternal(std::make_integer_sequence<int,N>{},std::forward<F>(f));
  }
  
#ifdef COMPILING_FOR_DEVICE
  
  /// Uses nvcc builtin
  ///
  /// \todo move it to a dedicated macro to call the proper bultin
#define UNROLLED_FOR(I,N)	\
  PRAGMA(unroll N)				\
  for(int I=0;I<N;I++)				\
    {
  
# define UNROLLED_FOR_END			\
  }
  
#else
  
  /// Create an unrolled for
  ///
  /// Hides the complexity
# define UNROLLED_FOR(I,N)			\
  unrolledFor<N>([&](const auto& I) INLINE_ATTRIBUTE {
  
  /// Finish an unrolled for
# define UNROLLED_FOR_END })
  
#endif
}

#endif
