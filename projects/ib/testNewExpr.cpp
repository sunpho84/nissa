#include <nissa.hpp>

using namespace nissa;

void test(Tensor<OfComps<Dir,ColorRow,ComplId>>& v,
	  Tensor<OfComps<Dir,ColorRow,ComplId>>& z)
{
  
  /// Scalar product on color index
  Tensor<OfComps<Dir>> res1;
  ASM_BOOKMARK_BEGIN("scalProd");
  res1=real(dag(v)*z);
  ASM_BOOKMARK_END("scalProd");
  
  master_printf("%lg\n",res1(Dir(0))); //55
  //master_printf("%lg\n",res1(Dir(1))); //50
  
  /// Outer product on color index
  Tensor<OfComps<Dir,ColorRow,ColorCln,ComplId>> res2;
  ASM_BOOKMARK_BEGIN("outerProd");
  res2=v*dag(z);
  ASM_BOOKMARK_END("outerProd");
  
  master_printf("%lg\n",res2(Dir(0),ColorRow(0),ColorCln(0),Re)); //1
  
  /// Comp by comp product on color index
  Tensor<OfComps<Dir,ColorRow,ComplId>> res3;
  ASM_BOOKMARK_BEGIN("compByCompProd");
  res3=v*z;
  ASM_BOOKMARK_END("compByCompProd");
  
  master_printf("%lg\n",res3(Dir(0),ColorRow(0),Re)); //-1
}

void in_main(int narg,char** arg)
{
  Tensor<OfComps<Dir,ColorRow,ComplId>> v;
  FOR_ALL_DIRS(mu)
    FOR_ALL_ROW_COLORS(cr)
    FOR_REIM_PARTS(ri)
    v(mu,cr,ri)=ri()+2*(cr()+NCOL*mu());
  
  test(v,v);
}

int main(int narg,char** arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
