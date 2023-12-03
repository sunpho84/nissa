#ifndef _QCD_HPP
#define _QCD_HPP

#include <expr/comp.hpp>

namespace nissa
{
  constexpr int nDirac=4;
  
  constexpr int nCol=3;
  
  DECLARE_TRANSPOSABLE_COMP(Spin,int,nDirac,spin);
  DECLARE_TRANSPOSABLE_COMP(Color,int,nCol,color);
  DECLARE_TRANSPOSABLE_COMP(Fuf,int,1,fuf);
}

#endif
