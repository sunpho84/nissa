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
    static constexpr I sizeAtCompileTime=
      SIZE;
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
  /// The specified name is suffixed with "Signature", to allow the
  /// Definition of the actual component with the expected name
#define DECLARE_COMPONENT_SIGNATURE(NAME,TYPE,LENGTH)		\
  /*! Signature for the NAME component */			\
  struct NAME ## Signature :					\
    public TensorCompSize<TYPE,LENGTH>				\
  {								\
    /*! Type used for the index */				\
    using Index=						\
      TYPE;							\
  }
  
  /////////////////////////////////////////////////////////////////
  
  // Row, column and transposition
  
  /// Row or column
  enum RwCl{ROW,CLN,ANY};
  
  /// Transposed of a row or column
  ///
  /// Forward declaration
  template <RwCl>
  RwCl transp;
  
  /// Transposed of a row
  template <>
  inline constexpr RwCl transp<ROW> =
    CLN;
  
  /// Transposed of a column
  template <>
  inline RwCl transp<CLN> =
    ROW;
  
  /// Transposed of any
  template <>
  inline constexpr RwCl transp<ANY> =
    ANY;
  
  /////////////////////////////////////////////////////////////////
  
  // Component
  
  DEFINE_FEATURE(TensorCompFeat);
  
  /// Tensor component defined by base type S
  template <typename S,
	    RwCl RC=ROW,
	    int Which=0>
  struct TensorComp :
    public TensorCompFeat<TensorComp<S,RC,Which>>
  {
    /// Transposed component
    using Transp=
      TensorComp<S,transp<RC>,Which>;
    
    /// Base type
    using Base=
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
    static constexpr Base sizeAtCompileTime=
      Base::sizeAtCompileTime;
    
    /// Check if the size is known at compile time
    static constexpr
    bool sizeIsKnownAtCompileTime=
      sizeAtCompileTime!=DYNAMIC;
    
    /// Returns the size at compile time, with assert
    static constexpr Base sizeAtCompileTimeAssertingNotDynamic()
    {
      static_assert(sizeIsKnownAtCompileTime,"Size not known at compile time!");
      
      return sizeAtCompileTime;
    }
    
    /// Init from value
    template <typename T=Index,
	      ENABLE_THIS_TEMPLATE_IF(isSafeNumericConversion<Index,T>)>
    INLINE_FUNCTION CUDA_HOST_DEVICE
    constexpr TensorComp(T&& i=0) : i(i)
    {
    }
    
    /// Convert to actual reference with const attribute
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
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
    
    /// Assignment operator
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    TensorComp& operator=(const Index& oth)
    {
      i=oth;
      
      return
	*this;
    }
    
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
    
#undef PROVIDE_SELFOP_WITH_OTHER
    
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

  /// Promotes the argument i to a COMPONENT, through a function with given NAME
#define DECLARE_COMPONENT_FACTORY(NAME,COMPONENT...)		\
  template <typename T>						\
  static INLINE_FUNCTION constexpr CUDA_HOST_DEVICE		\
  COMPONENT NAME(T&& i)						\
  {								\
    return							\
      COMPONENT(i);						\
  }
  
  /// Declare a component with no special feature
  ///
  /// The component has no row/column tag or index, so it can be
  /// included only once in a tensor
#define DECLARE_COMPONENT(NAME,TYPE,SIZE/*,FACTORY*/)		\
  DECLARE_COMPONENT_SIGNATURE(NAME,TYPE,SIZE);			\
  								\
  /*! NAME component */						\
  using NAME=							\
    TensorComp<NAME ## Signature,ANY,0>/*;*/			\
								\
    /*DECLARE_COMPONENT_FACTORY(FACTORY,NAME)*/
  
  /// Declare a component which can be included more than once
  ///
  /// The component has a row/column tag, and an additional index, so
  /// it can be included twice in a tensor
#define DECLARE_ROW_OR_CLN_COMPONENT(NAME,TYPE,SIZE/*,FACTORY*/)	\
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
  for(TYPE NAME=0;NAME<TYPE::sizeAtCompileTimeAssertingIfDynamic();NAME++)
  
  /////////////////////////////////////////////////////////////////
  
  /// Collection of components
  template <typename...Tc>
  using TensorComps=
    std::tuple<Tc...>;
  
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
    struct _TensorCompsTransp;
    
    /// Transposes a list of components
    ///
    /// Actual implementation
    template <typename...S,
	      RwCl...RC,
	      int...Which>
    struct _TensorCompsTransp<TensorComps<TensorComp<S,RC,Which>...>>
    {
      /// Resulting type
      using type=
	TensorComps<TensorComp<S,transp<RC>,Which>...>;
    };
  }
  
  /// Transposes a list of components
  template <typename TC>
  using TensorCompsTransp=
    typename impl::_TensorCompsTransp<TC>::type;
}

#endif
