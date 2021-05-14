#ifndef _UNIVERSALREFERENCE_HPP
#define _UNIVERSALREFERENCE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// Universal constructor unprioritize
#define UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR ...

/// Forward the member of the Type, see https://en.cppreference.com/w/cpp/utility/forward
#define FORWARD_MEMBER_VAR(TYPE,OBJ,MEMBER_VAR)		\
  std::forward<decltype(std::forward<TYPE>(OBJ).MEMBER_VAR)>(OBJ.MEMBER_VAR)

 /// Provides template par list to unprioritize default SFINAE
  ///
  /// Use as last argument of a function overloaded by a other
  /// implementations using SFINAE to detect the proper version to be
  /// used. This has to be used in conjunction with the other macros
  /// UNPRIORITIZE_DEFAULT_VERSION_ARGS and
  /// UNPRIORITIZE_DEFAULT_VERSION_ARGS_CHECK as in this example
  ///
  /// \code
  /// template <typename T,
  ///           UNPRIORITIZE_DEFAULT_VERSION_TEMPLATE_PARS>
  /// int tell(T a,UNPRIORITIZE_DEFAULT_VERSION_ARGS)
  /// {
  ///   UNPRIORITIZE_DEFAULT_VERSION_ARGS_CHECK;
  ///   return 1;
  /// }
  ///
  /// template <typename T,
  ///           std::enable_if_t<sizeof(T)==4,void*> =nullptr>
  /// decltype(auto) tell(T&& a)
  /// {
  ///   return a+1;
  /// }
  ///
  /// int main()
  /// {
  ///    tell(1); //returns 2
  ///
  ///    return 0;
  /// }
  ///\endcode
#define UNPRIORITIZE_DEFAULT_VERSION_TEMPLATE_PARS                     \
  typename...DummyTypes       /* Fake list of types                  */
  
  /// Provide empty list of args, used to unprioritize default version
#define UNPRIORITIZE_DEFAULT_VERSION_ARGS      \
  DummyTypes...     /*< Fake list of args */
  
  /// Check that no extra arg is passed
#define UNPRIORITIZE_DEFAULT_VERSION_ARGS_CHECK        \
  static_assert(sizeof...(DummyTypes)==0,"Dummy parameters actually catched!")
  
namespace nissa
{
}

#endif
