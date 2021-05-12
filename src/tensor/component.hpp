#ifndef _COMPONENT_SIZE_HPP
#define _COMPONENT_SIZE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file componentSize.hpp
///
/// \brief Implements the size of a tensor component

#include <stdint.h>

#include <base/cuda.hpp>
#include <metaProgramming/feature.hpp>
#include <metaProgramming/inliner.hpp>
#include <metaProgramming/nonConstMethod.hpp>
#include <metaProgramming/sfinae.hpp>
#include <metaProgramming/tuple.hpp>
#include <metaProgramming/typeConversion.hpp>
#include <metaProgramming/unrolledFor.hpp>

namespace nissa
{
  /////////////////////////////////////////////////////////////////
  
  // Size
  
  /// Dynamic size
  constexpr int DYNAMIC=-1;
  
  /// Specify the size at compile time
  template <typename I,
	    I SIZE=DYNAMIC>
  struct TensorCompSize
  {
    /// Value beyond end
    CUDA_HOST_DEVICE
    static constexpr I sizeAtCompileTime()
    {
      return SIZE;
    }
  };
  
  namespace predicate
  {
    /// Predicate returning whether the size is known ow not at compile time
    template <bool Comp=true> // Value to be compared
    struct SizeIsKnownAtCompileTime
    {
      /// Internal implementation
      template <typename T>
      struct t
      {
	/// Predicate result
	static constexpr bool value=
	  (T::sizeIsKnownAtCompileTime==Comp);
      };
    };
  }
  
  /////////////////////////////////////////////////////////////////
  
  // Signature
  
  /// Define the signature for a componenent convertible to TYPE of given NAME and SIZE
  ///
  /// Common part
#define _COMPONENT_SIGNATURE_BODY(NAME,TYPE,LENGTH)		\
    public TensorCompSize<TYPE,LENGTH>				\
  {								\
    /*! Type used for the index */				\
    using Index=						\
      TYPE;							\
  }
  
  /// Define the signature for a componenent convertible to TYPE of given NAME and SIZE
  ///
  /// The specified name is suffixed with "Signature", to allow the
  /// Definition of the actual component with the expected name
#define DECLARE_COMPONENT_SIGNATURE(NAME,TYPE,LENGTH)		\
  /*! Signature for the NAME component */			\
  struct NAME ## Signature :					\
  _COMPONENT_SIGNATURE_BODY(NAME,TYPE,LENGTH)
  
  // Signature
  
  /// Define the signature for a componenent convertible to TYPE of given NAME and SIZE
  ///
  /// Templated case
#define DECLARE_TEMPLATED_COMPONENT_SIGNATURE(NAME,LENGTH)	\
  /*! Signature for the NAME component */			\
  template <typename T>						\
  struct NAME ## Signature :					\
  _COMPONENT_SIGNATURE_BODY(NAME,T,LENGTH)
  
  /////////////////////////////////////////////////////////////////
  
  // Row, column and transposition
  
  /// Row or column
  enum RwCl{ROW,CLN,ANY};
  
  /// Transposed of a row or column
  ///
  /// Forward declaration
  template <RwCl>
  RwCl transpRwCl;
  
  /// Transposed of a row
  template <>
  inline constexpr RwCl transpRwCl<ROW> =
    CLN;
  
  /// Transposed of a column
  template <>
  inline constexpr RwCl transpRwCl<CLN> =
    ROW;
  
  /// Transposed of any
  template <>
  inline constexpr RwCl transpRwCl<ANY> =
    ANY;
  
  /////////////////////////////////////////////////////////////////
  
  // Component
  
  DEFINE_FEATURE(TensorCompFeat);
  
  /// Tensor component defined by signature type S
  template <typename S,
	    RwCl RC=ROW,
	    int Which=0>
  struct TensorComp :
    public TensorCompFeat<TensorComp<S,RC,Which>>
  {
    /// Transposed component
    using Transp=
      TensorComp<S,transpRwCl<RC>,Which>;
    
    /// Non column version
    using NonCol=
      TensorComp<S,((RC==ANY)?ANY:ROW),Which>;
    
    /// Signature type
    using Signature=
      S;
    
    /// Value type
    using Index=
      typename S::Index;
    
    /// Value
    Index i;
    
    /// Row or column
    static constexpr
    RwCl rC=
      RC;
    
    /// Index of the component
    static constexpr
    int which=
      Which;
    
    /// Size at compile time
    CUDA_HOST_DEVICE
    static constexpr Index sizeAtCompileTime()
    {
      return Signature::sizeAtCompileTime();
    }
    
    /// Check if the size is known at compile time
    static constexpr
    bool sizeIsKnownAtCompileTime=
      sizeAtCompileTime()!=DYNAMIC;
    
    /// Returns the size at compile time, with assert
    static constexpr Index sizeAtCompileTimeAssertingNotDynamic()
    {
      static_assert(sizeIsKnownAtCompileTime,"Size not known at compile time!");
      
      return sizeAtCompileTime();
    }
    
    /// Default constructor
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    TensorComp()
    {
    }
    
    /// Init from value
    template <typename T=Index,
	      ENABLE_THIS_TEMPLATE_IF(isSafeNumericConversion<Index,T>)>
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    TensorComp(T&& i) : i(i)
    {
    }
    
    /// Convert to actual reference with const attribute
    INLINE_FUNCTION constexpr CUDA_HOST_DEVICE
    const Index& toPod() const
    {
      return i;
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(toPod,CUDA_HOST_DEVICE);
    
    /// Convert to actual reference with const attribute
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    explicit operator const Index&() const
    {
      return toPod();
    }
    
    /// Convert to actual reference
    INLINE_FUNCTION constexpr CUDA_HOST_DEVICE
    explicit operator Index&()
    {
      return toPod();
    }
    
    /// Convert to actual reference with const attribute
    INLINE_FUNCTION constexpr CUDA_HOST_DEVICE
    const Index& operator()() const
    {
      return toPod();
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(operator(),CUDA_HOST_DEVICE);
    
    /// Convert to actual reference with const attribute, to be remoed
    INLINE_FUNCTION constexpr CUDA_HOST_DEVICE
    const Index& nastyConvert() const
    {
      return toPod();
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(nastyConvert,CUDA_HOST_DEVICE);
    
    /// Transposed index
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    auto transp()
      const
    {
      return
	Transp{i};
    }
    
    /// Dagger index, alias for transposed
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    auto dag()
      const
    {
      return
	transp();
    }
    
    /// Forbid assignement to a temporary
    TensorComp& operator=(const TensorComp& oth) && = delete;
    
    /// Assignment operator
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    TensorComp& operator=(const Index& oth) &
    {
      i=oth;
      
      return
	*this;
    }
    
    /// Assignment operator of a TensComp
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    TensorComp& operator=(const TensorComp& oth) &
    {
      return
	(*this)=oth.i;
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_POSTFIX_OP(OP)				\
    /*! Self OP operator */				\
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr		\
    TensorComp operator OP(int)				\
    {							\
      /*! Keep result */				\
      TensorComp res=i;					\
							\
      i OP;						\
      							\
      return						\
	res;						\
    }
    
    PROVIDE_POSTFIX_OP(++)
    PROVIDE_POSTFIX_OP(--)
    
#undef PROVIDE_POSTFIX_OP
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_SELFOP_WITH_OTHER(OP)			\
    /*! Assignment OP operator */			\
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr		\
    TensorComp& operator OP ## =(const TensorComp& oth)	\
      {							\
	i OP ## =oth.i;					\
							\
	return						\
	*this;						\
    }
    
    PROVIDE_SELFOP_WITH_OTHER(+)
    PROVIDE_SELFOP_WITH_OTHER(-)
    PROVIDE_SELFOP_WITH_OTHER(*)
    PROVIDE_SELFOP_WITH_OTHER(/)
    PROVIDE_SELFOP_WITH_OTHER(%)
    
#undef PROVIDE_SELFOP_WITH_OTHER
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_SELF_OP(OP)				\
    /*! Self OP operator */				\
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr		\
    TensorComp operator OP()				\
    {							\
      return OP i;						\
    }
    
    PROVIDE_SELF_OP(+)
    PROVIDE_SELF_OP(-)
    
#undef PROVIDE_SELF_OP
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_BINARY_OP(OP,RESULT_TYPE...)				\
    /*! Combine with other */						\
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr				\
    friend RESULT_TYPE operator OP(const TensorComp& first,const TensorComp& second) \
    {									\
      return								\
	first.i OP second.i;						\
    }
    
    PROVIDE_BINARY_OP(+,TensorComp<S,RC,Which>);
    PROVIDE_BINARY_OP(-,TensorComp<S,RC,Which>);
    PROVIDE_BINARY_OP(*,TensorComp<S,RC,Which>);
    PROVIDE_BINARY_OP(/,TensorComp<S,RC,Which>);
    PROVIDE_BINARY_OP(%,TensorComp<S,RC,Which>);
    
    PROVIDE_BINARY_OP(==,bool);
    PROVIDE_BINARY_OP(!=,bool);
    PROVIDE_BINARY_OP(<,bool);
    PROVIDE_BINARY_OP(<=,bool);
    PROVIDE_BINARY_OP(>,bool);
    PROVIDE_BINARY_OP(>=,bool);
    
#undef PROVIDE_BINARY_OP
  };
  
  /// Declare a component with no special feature
  ///
  /// The component has no row/column tag or index, so it can be
  /// included only once in a tensor
#define DECLARE_COMPONENT(NAME,TYPE,SIZE)			\
  DECLARE_COMPONENT_SIGNATURE(NAME,TYPE,SIZE);			\
  								\
  /*! NAME component */						\
  using NAME=							\
    TensorComp<NAME ## Signature,ANY,0>
  
  /// Declare a template component
  ///
  /// The component has no row/column tag or index, so it can be
  /// included only once in a tensor
#define DECLARE_TEMPLATED_COMPONENT(NAME,SIZE)			\
  DECLARE_TEMPLATED_COMPONENT_SIGNATURE(NAME,SIZE);		\
  								\
  /*! NAME component */						\
  template <typename T>						\
  using NAME=							\
    TensorComp<NAME ## Signature<T>,ANY,0>
  
  /// Declare a component which can be included more than once
  ///
  /// The component has a row/column tag, and an additional index, so
  /// it can be included twice in a tensor
#define DECLARE_ROW_OR_CLN_COMPONENT(NAME,TYPE,SIZE)		\
  DECLARE_COMPONENT_SIGNATURE(NAME,TYPE,SIZE);			\
  								\
  /*! NAME component */						\
  template <RwCl RC=ROW,					\
	    int Which=0>					\
  using NAME ## RC=TensorComp<NAME ## Signature,RC,Which>;	\
								\
  /*! Row kind of NAME component */				\
  using NAME ## Row=NAME ## RC<ROW,0>;				\
								\
  /*! Column kind of NAME component */				\
  using NAME ## Cln=NAME ## RC<CLN,0>;				\
  								\
  /*! Default NAME component is Row */				\
  using NAME=NAME ## Row/*;*/					\
  								\
    /*DECLARE_COMPONENT_FACTORY(FACTORY ## Row,NAME ## Row);*/	\
								\
    /*DECLARE_COMPONENT_FACTORY(FACTORY ## Cln,NAME ## Cln);*/	\
								\
    /*DECLARE_COMPONENT_FACTORY(FACTORY,NAME)*/
  
  /////////////////////////////////////////////////////////////////
  
#define FOR_ALL_COMPONENT_VALUES(TYPE,NAME)				\
  for(TYPE NAME=0;NAME<TYPE::sizeAtCompileTimeAssertingNotDynamic();NAME++)
  
#define UNROLLED_FOR_COMPONENT_VALUES_IN_RANGE(TYPE,NAME,MIN,MAX,CORE...) \
  unrolledFor<MIN,MAX>([&](const TYPE& NAME) INLINE_ATTRIBUTE \
  {									\
    CORE;								\
  })
  
#define FOR_ALL_COMPONENT_VALUES_STARTING_AT(TYPE,NAME,MIN,CORE...)	\
  FOR_COMPONENT_VALUES_IN_RANGE(TYPE,NAME,MIN,TYPE::sizeAtCompileTimeAssertingNotDynamic(),CORE)
  
#define UNROLLED_FOR_ALL_COMPONENT_VALUES_STARTING_AT(TYPE,NAME,MIN,CORE...)	\
  UNROLLED_FOR_COMPONENT_VALUES_IN_RANGE(TYPE,NAME,MIN,TYPE::sizeAtCompileTimeAssertingNotDynamic(),CORE)
  
#define FOR_ALL_COMPONENT_VALUES(TYPE,NAME,CORE...)			\
  FOR_ALL_COMPONENT_VALUES_STARTING_AT(TYPE,NAME,0,CORE)
  
#define UNROLLED_FOR_ALL_COMPONENT_VALUES(TYPE,NAME,CORE...)			\
  UNROLLED_FOR_ALL_COMPONENT_VALUES_STARTING_AT(TYPE,NAME,0,CORE)
  
  /////////////////////////////////////////////////////////////////
  
  /// Collection of components
  template <typename...Tc>
  using TensorComps=
    std::tuple<Tc...>;
  
  /// Alias to make it easier to understand tensor instantiation
  template <typename...Tc>
  using OfComps=
    TensorComps<Tc...>;
  
  /// Returns the number of components of a tensComp
  template <typename T>
  constexpr int nOfComps=
    std::tuple_size<typename T::Comps>::value;
  
  namespace impl
  {
    /// Provides the result of filtering from a list of components the Row or Column
    ///
    /// Forward definition
    template <RwCl F,
	      typename TC>
    struct _TensorCompsFilterRwCl;
    
    /// Cannot use directly the TupleFilter, because of some template template limitation
    template <RwCl F,
	      typename...Tc>
    struct _TensorCompsFilterRwCl<F,TensorComps<Tc...>>
    {
      /// Predicate to filter
      ///
      /// Forward definition
      template <typename T>
      struct Filter;
      
      /// Predicate to filter
      template <typename S,
		RwCl RC,
		int Which>
      struct Filter<TensorComp<S,RC,Which>>
      {
	/// Predicate result, counting whether the type match
	static constexpr
	bool value=
	  (RC==F);
      };
      
      /// Returned type
      typedef TupleFilter<Filter,TensorComps<Tc...>> type;
    };
  }
  
  /// Filter all Row components
  template <typename TC>
  using TensorCompsFilterRow=
    typename impl::_TensorCompsFilterRwCl<RwCl::ROW,TC>::type;
  
  /// Filter all Col components
  template <typename TC>
  using TensorCompsFilterCln=
    typename impl::_TensorCompsFilterRwCl<RwCl::CLN,TC>::type;
  
  /// Filter all Any components
  template <typename TC>
  using TensorCompsFilterAny=
    typename impl::_TensorCompsFilterRwCl<RwCl::ANY,TC>::type;
  
  /// Gets the dynamic components of a tensComps
  template <typename TC>
  constexpr decltype(auto) getDynamicCompsOfTensorComps(TC&& tc)
  {
    return tupleFilter<predicate::SizeIsKnownAtCompileTime<false>::t>(std::forward<TC>(tc));
  }
  
  /// Gets the dynamic component types of a TensorComps
  template <typename TC>
  using GetDynamicCompsOfTensorComps=
    decltype(getDynamicCompsOfTensorComps(TC{}));
  
  /// Gets the fixed size components of a tensComps
  template <typename TC>
  constexpr decltype(auto) getFixedSizeCompsOfTensorComps(TC&& tc)
  {
    return tupleFilter<predicate::SizeIsKnownAtCompileTime<true>::t>(std::forward<TC>(tc));
  }
  
  /// Gets the fixed size component types of a TensorComps
  template <typename TC>
  using GetFixedSizeCompsOfTensorComps=
    decltype(getFixedSizeCompsOfTensorComps(TC{}));
  
  /////////////////////////////////////////////////////////////////
  
  namespace impl
  {
    /// Transposes a list of components
    ///
    /// Actual implementation, forward declaration
    template <typename TC>
    struct _TransposeTensorComps;
    
    /// Transposes a list of components
    ///
    /// Actual implementation
    template <typename...TC>
    struct _TransposeTensorComps<TensorComps<TC...>>
    {
      /// Returns a given components, or its transposed if it is missing
      template <typename C,
		typename TranspC=typename C::Transp>
      using ConditionallyTransposeComp=
	std::conditional_t<tupleHasType<TensorComps<TC...>,TranspC,1>,C,TranspC>;
      
      /// Resulting type
      using type=
	TensorComps<ConditionallyTransposeComp<TC>...>;
    };
  }
  
  /// Transposes a list of components
  ///
  /// - If a component is not of ROW/CLN case, it is left unchanged
  /// - If a ROW/CLN component is matched with a CLN/ROW one, it is left unchanged
  /// - If a ROW/CLN component is not matched, it is transposed
  ///
  /// \example
  ///
  /// using T=TensorComps<Complex,ColorRow,ColorCln,SpinRow>
  /// using U=TransposeTensorcomps<T>; //TensorComps<Complex,ColorRow,ColorCln,SpinCln>
  template <typename TC>
  using TransposeTensorComps=
    typename impl::_TransposeTensorComps<TC>::type;
}

#endif
