#include <nissa.hpp>

using namespace nissa;

void in_main(int narg,char** arg)
{
  Tensor<OfComps<Dir,ColorRow>> rt;
  FOR_ALL_DIRS(mu)
    FOR_ALL_ROW_COLORS(cr)
      rt(mu,cr)=cr()+NCOL*mu();
  
  Tensor<OfComps<Dir,ColorCln>> ct;
  FOR_ALL_DIRS(mu)
    FOR_ALL_CLN_COLORS(cc)
      ct(mu,cc)=cc()+NCOL*mu();
  
  /// Scalar product on color index
  Tensor<OfComps<Dir>> res;
  res=ct*rt;
}

int main(int narg,char** arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
