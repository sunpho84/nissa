#ifndef _QCD_HPP
#define _QCD_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <expr/comp.hpp>

namespace nissa
{
#define NCOL 3

#define NSPIN 4
  
  DECLARE_TRANSPOSABLE_COMP(Spin,int,NSPIN,spin);
  DECLARE_TRANSPOSABLE_COMP(Color,int,NCOL,color);
  DECLARE_TRANSPOSABLE_COMP(Fuf,int,1,fuf);
  
  constexpr Color nCol=NCOL;
  constexpr Spin nDirac=NSPIN;
}

#endif
