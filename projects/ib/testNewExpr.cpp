#include <nissa.hpp>

using namespace nissa;

void in_main(int narg,char** arg)
{
  Tensor<OfComps<Dir,ColorRow>> v;
  FOR_ALL_DIRS(mu)
    FOR_ALL_ROW_COLORS(cr)
      v(mu,cr)=cr()+NCOL*mu();
  
  /// Scalar product on color index
  Tensor<OfComps<Dir>> res1;
  res1=dag(v)*v;
  
  master_printf("%lg\n",res1(Dir(0))); //5
  master_printf("%lg\n",res1(Dir(1))); //50
  
  /// Outer product on color index
  Tensor<OfComps<Dir,ColorRow,ColorCln>> res2;
  res2=v*dag(v);
  
  /// Comp by comp product on color index
  Tensor<OfComps<Dir,ColorRow>> res3;
  res3=v*v;
}

int main(int narg,char** arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
