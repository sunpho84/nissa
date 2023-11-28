#ifndef _TUPLEDESCRIBE_HPP
#define _TUPLEDESCRIBE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tuples/tupleDescribe.hpp

#include <tuple>
#include <string>
#include <sstream>

#include <base/debug.hpp>
#include <tuples/forEachIndexOfTuple.hpp>

namespace nissa
{
  /// Describe tuple TP
  template <typename TP>
  std::string tupleDescribe(const TP& t)
  {
    std::ostringstream os;
    os<<"(";
    forEachIndexOfTuple<TP>([&os,&t]<size_t I>(){os<<std::string((I==0)?"":",")<<demangle(typeid(std::tuple_element<I,TP>).name())<<"="<<(std::string)std::get<I>(t);});
    os<<"(";
    
    return os.str();
  }
}

#endif
