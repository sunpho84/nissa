#ifndef _QCD_HPP
#define _QCD_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <expr/comp.hpp>

namespace nissa
{
#define NCOL 3
  
  DECLARE_TRANSPOSABLE_COMP(Color,int,NCOL,color);
  
  constexpr Color nCol=NCOL;
}

#endif
