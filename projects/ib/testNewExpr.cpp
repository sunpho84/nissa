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

template <bool IsComplProd,
	  typename PComp,
	  typename...FirstComps,
	  typename...ContractedComps,
	  size_t...I>
constexpr int process(TensorComps<FirstComps...>,
			 TensorComps<ContractedComps...>,
			 std::index_sequence<I...>)
{
  constexpr bool isContracted=
    (std::is_same_v<typename ContractedComps::Transp,PComp>||...);
  
  constexpr bool isComplInComplProd=
	      std::is_same_v<ComplId,PComp> and IsComplProd;
  
  constexpr int pos=
	      ((std::is_same_v<PComp,FirstComps>?(I+1):0)+...)-1;
  
  constexpr int res=
	      (isContracted or isComplInComplProd)?
	      -1:
              pos;
  
  return
    res;
}

template <bool IsComplProd,
	  typename...PComps,
	  typename...FirstComps,
	  typename...ContractedComps>
constexpr auto process(TensorComps<PComps...>,
		       TensorComps<FirstComps...>,
		       TensorComps<ContractedComps...>)
{
  return std::integer_sequence<int,
			       process<IsComplProd,PComps>(TensorComps<FirstComps...>{},
                                                           TensorComps<ContractedComps...>{},
                                                           std::make_index_sequence<sizeof...(FirstComps)>{})...>{};
}

template <typename P>
constexpr auto process()
{
  using PC=
    typename P::Comps;
  
  using C1=
    typename P::template NestedExpr<0>::Comps;
  
  using CC=
    typename P::ContractedComps;
  
  return process<P::isComplProd>(PC{},C1{},CC{});
}

void test3(Tensor<OfComps<SpinRow,ColorRow,ColorCln,SpinCln,ComplId,LocLxSite>>& v,
	   Tensor<OfComps<ColorRow,ColorCln,SpinRow,ComplId,LocLxSite>>& z)
{
  const auto p=v*z;
  using P=decltype(p);

  const std::integer_sequence<int, 0, 1, -1, -1, 5> pr = process<P>();
}

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
