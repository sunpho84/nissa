#ifndef _OPERATOREXISTS_HPP
#define _OPERATOREXISTS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

/// \file operatorExists.hpp

namespace nissa
{
  /// Check that the binary operator exists
#define DEFINE_BINARY_OPERATOR_IMPLEMENTATION_CHECK(CHECK_NAME,         \
                                                    STRUCT_NAME,        \
                                                    OPERATOR)           \
  /*! Structure used to check if the operator is implemented */         \
  template <typename S,                                                 \
            typename T>                                                 \
  struct STRUCT_NAME                                                    \
  {                                                                     \
    /*! Path followed when the operator is implemented */               \
    template <typename U,                                               \
              typename V>                                               \
    static auto test(U*)->decltype((*std::declval<U*>())		\
                                     OPERATOR                           \
                                     std::declval<V>());                \
                                                                        \
    /*! Default case in which the binary operation cannot be performed */ \
    template <typename,                                                 \
              typename>                                                 \
      static auto test(...)->std::false_type;                           \
                                                                        \
    /*! Result of the check */                                          \
    static constexpr bool res=                                          \
      not std::is_same_v<std::false_type,decltype(test<S,T>(0))>;	\
  };                                                                    \
                                                                        \
  /*! Check that operator OPERATOR is implemented */                    \
  template <typename S,                                                 \
            typename T>                                                 \
  [[ maybe_unused ]]                                                    \
  constexpr bool CHECK_NAME=                                            \
    STRUCT_NAME<S,T>::res
}

#endif
