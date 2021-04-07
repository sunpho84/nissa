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
///           SFINAE_ON_TEMPLATE_ARG(IsSame<D,int>)>
/// void foo(D i) {} // fails if D is not int
/// \endcode
#define SFINAE_ON_TEMPLATE_ARG(...)	\
  std::enable_if_t<(__VA_ARGS__),void*> =nullptr

#endif
