#ifndef _DEMANGLE_HPP
#define _DEMANGLE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#if CAN_DEMANGLE
# include <cxxabi.h>
#endif

#include <string>

namespace nissa
{
  /// Demangle a string
  ///
  /// If the compiler has no abi functionality, the original string is
  /// returned.
  inline std::string demangle(const std::string& what)
  {
#if CAN_DEMANGLE
    /// Returned status of demangle
    int status=1;
    
    /// Demangled
    char* name=
      abi::__cxa_demangle(what.c_str(),0,0,&status);
    
    /// Copy the result
    std::string out=
      (status==0)?
      name:
      (what+" (failed demangle)");
    
    // Free if succeded
    if(status==0)
      free(name);
    
    return
      out;
#else
    return what;
#endif
  }
}

#endif
