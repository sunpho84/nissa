#include <nissa.hpp>

using namespace nissa;

void test(Tensor<OfComps<LocLxSite,ColorRow,ComplId>>& v,
	  Tensor<OfComps<LocLxSite,ColorRow,ComplId>>& z)
{
  /// Scalar product on color index
  Tensor<OfComps<LocLxSite>> res1(locVol);
  ASM_BOOKMARK_BEGIN("scalProd");
  res1=real(dag(v)*z);
  ASM_BOOKMARK_END("scalProd");
  
  master_printf("%lg\n",res1(LocLxSite(0))); //55
  
  /// Outer product on color index
  Tensor<OfComps<LocLxSite,ColorRow,ColorCln,ComplId>> res2(locVol);
  ASM_BOOKMARK_BEGIN("outerProd");
  res2=v*dag(z);
  ASM_BOOKMARK_END("outerProd");
  
  master_printf("%lg\n",res2(LocLxSite(0),ColorRow(0),ColorCln(0),Re)); //1
  
  /// Comp by comp product on color index
  Tensor<OfComps<LocLxSite,ColorRow,ComplId>> res3(locVol);
  ASM_BOOKMARK_BEGIN("compByCompProd");
  res3=v*z;
  ASM_BOOKMARK_END("compByCompProd");
  
  master_printf("%lg\n",res3(LocLxSite(0),ColorRow(0),Re)); //-1
}

void test2(Tensor<OfComps<LocLxSite,ColorRow,ColorCln>>& v,
	   Tensor<OfComps<LocLxSite,ColorRow,ColorCln>>& z)
{
  /// Matrix product
  Tensor<OfComps<LocLxSite>> res1;
  res1.allocate(locVol);
  ASM_BOOKMARK_BEGIN("scalProd");
  Tensor<OfComps<LocLxSite,ColorRow,ColorCln>> c;
  c=v*z;
  ASM_BOOKMARK_END("scalProd");
}

DECLARE_ROW_OR_CLN_COMPONENT(Spin,int,NDIRAC);


void in_main(int narg,char** arg)
{
  init_grid(4,4);
  
  Tensor<OfComps<LocLxSite,ColorRow,ComplId>> v;
  v.allocate(locVol);
  
  NISSA_LOC_VOL_LOOP(site)
    FOR_ALL_ROW_COLORS(cr)
    FOR_REIM_PARTS(ri)
    v(site,cr,ri)=ri()+2*(cr()+NCOL*site());
  
  test(v,v);
}

int main(int narg,char** arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
