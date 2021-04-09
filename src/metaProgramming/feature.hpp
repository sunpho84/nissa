#ifndef _FEATURE_HPP
#define _FEATURE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file feature.hpp
///
/// \brief Implements static polymorphism

#include <base/cuda.hpp>
#include <metaProgramming/nonConstMethod.hpp>

namespace nisssa
{
  /// Define a feature
#define DEFINE_FEATURE(FEATURE_NAME)			\
  template <typename Defeated>				\
  struct FEATURE_NAME					\
  {							\
    PROVIDE_DEFEAT_METHOD(Defeated);			\
  }
  
  /// Provides the method to cast to the featuring class
#define PROVIDE_DEFEAT_METHOD(T)			\
  /*! Cast to the base type, with const attribute */	\
  CUDA_HOST_DEVICE					\
  operator const T&() const				\
  {							\
    return						\
      *static_cast<const T*>(this);			\
  }							\
  							\
  /*! Cast to the base type */						\
  /*!                       */						\
  /*! Cannot be achieved with the preprocessor macro, since */		\
  /*! the name of the method is weird */				\
  CUDA_HOST_DEVICE							\
  constexpr operator T&()						\
  {									\
    return *static_cast<T*>(this);					\
  }									\
  									\
  /*! Cast to the featuring class */					\
  CUDA_HOST_DEVICE							\
  constexpr const T& deFeat() const					\
  {									\
    return *this;							\
  }									\
  									\
  PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(deFeat,CUDA_HOST_DEVICE)
  
  /// Import method from the feature class
#define IMPORT_FEATURE_METHOD(A...)				\
  /*! Calls A in the base class */				\
  template <typename...Args>					\
  CUDA_HOST_DEVICE						\
  decltype(auto) A(Args&&...args) const				\
  {								\
    return (*this)().A(std::forward<Args>(args)...);		\
  }
}

#endif
