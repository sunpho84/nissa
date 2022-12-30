#ifndef _CRTP_HPP
#define _CRTP_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file crtp.hpp

namespace grill
{
#define DE_CRTPFY(TYPE,PTR)			\
  (*static_cast<TYPE*>(PTR))
}

#endif
