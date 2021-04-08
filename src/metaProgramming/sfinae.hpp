#ifndef _SFINAE_HPP
#define _SFINAE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// Provides a SFINAE to be used in template par list
///
/// This follows
/// https://stackoverflow.com/questions/32636275/sfinae-with-variadic-templates
/// as in this example
/// \code
/// template <typename D,
///           ENABLE_THIS_TEMPLATE_IF(std::is_same<D,int>::value)>
/// void foo(D i) {} // fails if D is not int
/// \endcode
#define ENABLE_THIS_TEMPLATE_IF(...)			\
  std::enable_if_t<(__VA_ARGS__),void*> =nullptr

#endif
