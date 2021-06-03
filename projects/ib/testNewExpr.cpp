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

namespace nissa
{
  namespace internal
  {
    /// Computes the Component melding barriers starting from the knowledge of in-out comps and the in barriers
    ///
    /// Internal implementation, forward declaration
    template <typename TcOut,
	      typename TcIn1,
	      typename TcIn2,
	      typename MBsIn1,
	      typename MBsIn2,
	      typename OutPositionsInIn1,
	      typename OutPositionsInIn2,
	      typename II=std::make_index_sequence<std::tuple_size_v<TcOut>-1>>
    struct _GetTensorCompsMeldBarriersForProd;
    
    /// Computes the Component melding barriers starting from the knowledge of in-out comps and the in barriers
    ///
    /// Internal implementation
    template <typename...TcOut,
	      typename...TcIns1,
	      typename...TcIns2,
	      size_t...MBsIns1,
	      size_t...MBsIns2,
	      size_t...OutPositionsInIn1,
	      size_t...OutPositionsInIn2,
	      size_t...IIs>
    struct _GetTensorCompsMeldBarriersForProd<TensorComps<TcOut...>,
					      TensorComps<TcIns1...>,
					      TensorComps<TcIns2...>,
					      TensorCompsMeldBarriers<MBsIns1...>,
					      TensorCompsMeldBarriers<MBsIns2...>,
					      std::index_sequence<OutPositionsInIn1...>,
					      std::index_sequence<OutPositionsInIn2...>,
					      std::index_sequence<IIs...>>
    {
      /// Positions of output components as found in the input one into an array (for whatever reason, a c-array is not welcome)
      static constexpr std::array<size_t,sizeof...(OutPositionsInIn1)> outPositionsInIn1=
	{OutPositionsInIn1...};
      
      /// Positions of output components as found in the input one into an array (for whatever reason, a c-array is not welcome)
      static constexpr std::array<size_t,sizeof...(OutPositionsInIn2)> outPositionsInIn2=
	{OutPositionsInIn2...};
      
      /// Put the positions of input barriers into an array
      static constexpr size_t mBsIns1[]=
	{MBsIns1...};
      
      /// Put the positions of input barriers into an array
      static constexpr size_t mBsIns2[]=
	{MBsIns2...};
      
      /// Check if a barrier must be included between IPrev and IPrev+1
      template <size_t IPrev>
      static constexpr bool insertBarrier=
	(outPositionsInIn1[IPrev]+1!=outPositionsInIn1[IPrev+1] and not
	 (outPositionsInIn1[IPrev]==outPositionsInIn1[IPrev] and outPositionsInIn1[IPrev]==sizeof...(TcIns1))) or
	((outPositionsInIn1[IPrev+1]==MBsIns1)||...) or
	(outPositionsInIn2[IPrev]+1!=outPositionsInIn2[IPrev+1] and not
	 (outPositionsInIn2[IPrev]==outPositionsInIn2[IPrev] and outPositionsInIn2[IPrev]==sizeof...(TcIns2))) or
	((outPositionsInIn2[IPrev+1]==MBsIns2)||...);
      
      /// If a barrier is needed, returns a tuple with the integral constant, otherwise an empty one
      template <size_t IPrev>
      using OptionalBarrier=
	std::conditional_t<insertBarrier<IPrev>,
			   std::tuple<std::integral_constant<size_t,IPrev+1>>,
			   std::tuple<>>;
      
      /// Put together the possible barriers in a single tuple
      using BarriersInATuple=
	TupleCat<OptionalBarrier<IIs>...>;
      
      /// Resulting type obtained flattening the tuple of optional barriers
      using type=
	TupleOfIntegralConstantsToIntegerSequence<BarriersInATuple,size_t>;
    };
  }
}

template <typename PComp,
	  typename...Comps,
	  typename...ExcludedComps,
	  size_t...I>
constexpr size_t process(TensorComps<Comps...>,
			 TensorComps<ExcludedComps...>,
			 std::index_sequence<I...>)
{
  constexpr bool isExcluded=
    (std::is_same_v<ExcludedComps,PComp>||...);
  
  constexpr size_t pos=
	      firstOccurrenceOfType<PComp>(Comps{}...);
  
  constexpr size_t res=
	      isExcluded?
	      sizeof...(I):
              pos;
  
  return
    res;
}

template <typename...PComps,
	  typename...FirstComps,
	  typename...ContractedComps>
constexpr auto process(TensorComps<PComps...>,
		       TensorComps<FirstComps...>,
		       TensorComps<ContractedComps...>)
{
  return std::index_sequence<process<PComps>(TensorComps<FirstComps...>{},
    TensorComps<ContractedComps...>{},
    std::make_index_sequence<sizeof...(FirstComps)>{})...>{};
}

template <typename P>
constexpr auto process()
{
  using PC=
    typename P::Comps;
  
  using F1=
    typename P::template NestedExpr<0>;
  
  using C1=
    typename F1::Comps;
  
  using M1=
    typename F1::CompsMeldBarriers;
  
  using F2=
    typename P::template NestedExpr<1>;
  
  using C2=
    typename F2::Comps;
  
  using M2=
    typename F2::CompsMeldBarriers;
  
  using CC2=
    typename P::ContractedComps;
  
  using CC1=
    TransposeTensorComps<CC2>;
  
  using P1=
    decltype(process(PC{},C1{},CC1{}));
  
  using P2=
    decltype(process(PC{},C2{},CC2{}));

  using MM=typename internal::_GetTensorCompsMeldBarriersForProd<PC,C1,C2,M1,M2,P1,P2>::type;
  
  return MM{};
}


void test3(Tensor<OfComps<SpinRow,ColorRow,ColorCln,SpinCln,ComplId,LocLxSite>>& v,
	   Tensor<OfComps<ColorRow,ColorCln,SpinRow,ComplId,LocLxSite>>& z)
{
  const auto p=v*z;
  using P=decltype(p);

  const auto pf=process<P>();
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
