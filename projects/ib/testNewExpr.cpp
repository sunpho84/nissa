#include <nissa.hpp>

using namespace nissa;

void in_main(int narg,char** arg)
{
  Tensor<OfComps<Dir,ColorRow>> v;
  FOR_ALL_DIRS(mu)
    FOR_ALL_ROW_COLORS(cr)
      v(mu,cr)=cr()+NCOL*mu();
  
  /// Scalar product on color index
  Tensor<OfComps<Dir>> res;
  res=transp(v)*v;
  
  master_printf("%lg\n",res(Dir(0))); //5
  master_printf("%lg\n",res(Dir(1))); //50
}

int main(int narg,char** arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
