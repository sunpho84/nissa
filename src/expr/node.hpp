#ifndef _NODE_HPP
#define _NODE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/nodes/node.hpp
///
/// \brief Declare base node of the syntactic tree

#include <type_traits>

#include <base/debug.hpp>
#include <expr/bindComps.hpp>
#include <expr/compsMerger.hpp>
// #include <expr/comps/compLoops.hpp>
// #include <expr/assign/deviceAssign.hpp>
// #include <expr/assign/directAssign.hpp>
// #include <expr/assign/executionSpace.hpp>
// #include <expr/nodes/scalar.hpp>
#include <expr/assignDispatcher.hpp>
#include <expr/dynamicTensDeclaration.hpp>
#include <expr/nodeDeclaration.hpp>
// #include <expr/assign/simdAssign.hpp>
// #include <expr/assign/threadAssign.hpp>
// #include <ios/logger.hpp>
#include <metaprogramming/crtp.hpp>
#include <tuples/tupleHasType.hpp>
#include <tuples/tupleFilter.hpp>

namespace nissa
{
  namespace impl
  {
    /// Casts an expression to fund if possible
    template <typename E,
	      typename...C>
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    decltype(auto) castExprToFund(E&& t,
				  CompsList<C...>*)
    {
      if constexpr(((C::sizeAtCompileTime==1) and ...))
	return t.eval((C)0 ...);
    }
  }
  
#define THIS Node<T,CompsList<Ci...>>
  
  /// Node, the base type to be evaluated as an expression
  template <typename T,
	    typename...Ci>
  struct THIS :
    NodeFeat,
    Crtp<THIS,T>,
    MemberCompSubscribeProvider<T,Ci>...,
    DynamicCompsProvider<CompsList<Ci...>>
  {
    using This=THIS;
#undef THIS
    
#define PROVIDE_AUTOMATIC_CAST_TO_FUND(ATTRIB)			\
    /*! Provide automatic cast to fund if needed. How cool! */	\
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE		\
    operator decltype(auto)() ATTRIB				\
    {								\
      return impl::castExprToFund(DE_CRTPFY(ATTRIB T,this),	\
				  (typename T::Comps*)nullptr);	\
    }
    
    PROVIDE_AUTOMATIC_CAST_TO_FUND(const);
    
    PROVIDE_AUTOMATIC_CAST_TO_FUND(/* non const */);
    
#undef PROVIDE_AUTOMATIC_CAST_TO_FUND
    
    using Crtp<This,T>::operator~;
    
    // /// Define the move-assignment operator
    // INLINE_FUNCTION
    // Node& operator=(Node&& oth)
    // {
    //   return this->operator=<T>(std::forward<Node>(oth));
    // }
    
    /// Returns whether can assign: this is actually used when no other assignability is defined
    constexpr bool canAssign() const
    {
      return false;
    }
    
    /// Unless explicitly overloaded
    static constexpr bool canAssignAtCompileTime=false;
    
#define PROVIDE_MERGE_COMPS(ATTRIB,WHAT_TO_PASS)		\
    /*! Provides possibility to merge a list of components  */	\
      template <typename MCL>					\
	constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE		\
	decltype(auto) mergeComps() ATTRIB			\
      {								\
	return nissa::mergeComps<MCL>(WHAT_TO_PASS);		\
      }
    
    PROVIDE_MERGE_COMPS(const&,~*this);
    
    PROVIDE_MERGE_COMPS(/* non const */&,~*this);
    
    PROVIDE_MERGE_COMPS(/* non const */&&,std::move(~*this));
    
#undef PROVIDE_MERGE_COMPS
    
    /// Used to check that the derived type satisfy the Node criterion
    constexpr Node()
    {
      static_assert(assertIsNode<T>,
		    "Incomplete node type");
    }
    
    /// Assert assignability
    template <DerivedFromNode U>
      INLINE_FUNCTION
      constexpr void assertCanAssign(const U& rhs)
    {
      //static_assert(tuplesContainsSameTypes<typename T::Comps,typename U::Comps>,"Cannot assign two expressions which differ for the components");
      
      static_assert(T::canAssignAtCompileTime,"Trying to assign to a non-assignable expression");
      
      auto& lhs=~*this;
      
      if constexpr(not T::canAssignAtCompileTime)
	crash("Trying to assign to a non-assignable expression");
      
      if constexpr(U::hasDynamicComps)
	if(lhs.getDynamicSizes()!=rhs.getDynamicSizes())
	  crash("Dynamic comps not agreeing");
      
      //static_assert(T::execSpace==U::execSpace or
      //U::execSpace==ExecSpace::HOST_DEVICE,"Cannot assign among different execution space, first change one of them");
    }
    
    /// Assign from another expression
    template <typename OP=DirectAssign,
	      DerivedFromNode Rhs>
    constexpr INLINE_FUNCTION
      T& assign(const Rhs& u)
    {
      this->assertCanAssign(u);
      
// #if ENABLE_SIMD
//       if constexpr(T::canSimdify and
// 		   ((Rhs::canSimdify and std::is_same_v<typename T::SimdifyingComp,typename Rhs::SimdifyingComp>)
// 		    or not tupleHasType<typename Rhs::Comps,typename T::SimdifyingComp>))
// 	simdAssign(lhs,rhs);
//       else
// #endif
// #if ENABLE_DEVICE_CODE
// 	if constexpr(Rhs::execSpace==ExecSpace::DEVICE)
// 	  deviceAssign(lhs,rhs);
// 	else
// #endif
// #if ENABLE_THREADS
// 	  if constexpr(Rhs::nDynamicComps==1)
// 	    threadAssign(lhs,rhs);
// 	  else
// #endif
// 	    directAssign(lhs,rhs);
      
      compsLoop<typename T::Comps>([this,&u](const auto&...lhsComps)
      {
	const auto lhsCompsTup=std::make_tuple(lhsComps...);
	
	const auto rhsCompsTup=tupleGetSubset<typename Rhs::Comps>(lhsCompsTup);
	
	OP::dispatch((~*this)(lhsComps...),std::apply(~u,rhsCompsTup));
      },(~*this).getDynamicSizes());
      
      return ~*this;
    }
    
#define PROVIDE_ASSIGN_VARIATION(SYMBOL,OP)				\
    									\
    /*! Assign from another expression */				\
      template <DerivedFromNode Rhs>					\
	requires(not std::is_same_v<T,Rhs>)				\
	constexpr INLINE_FUNCTION					\
	T& operator SYMBOL(const Rhs& u)				\
      {									\
	return (~*this).template assign<OP>(u);				\
      }									\
      									\
    									\
      /*! Assign from a scalar */					\
      template <typename Oth>						\
	requires(std::is_arithmetic_v<Oth>)				\
	constexpr INLINE_FUNCTION					\
	T& operator SYMBOL(const Oth& value)				\
      {									\
	return (~*this) SYMBOL scalar(value);				\
      }
    
    PROVIDE_ASSIGN_VARIATION(=,DirectAssign);
    PROVIDE_ASSIGN_VARIATION(+=,SumAssign);
    PROVIDE_ASSIGN_VARIATION(-=,SubtAssign);
    
#undef PROVIDE_ASSIGN_VARIATION
    
    /// Define the assignment operator with the same expression type
    constexpr INLINE_FUNCTION
    Node& operator=(const Node& oth)
    {
      return (~*this).assign(~oth);
    }
    
    /// Gets a writeable reference
    constexpr INLINE_FUNCTION
    auto getWritable() &
    {
      return (~*this).getRef();
    }
    
    /// Gets a read-only reference
    constexpr INLINE_FUNCTION
    auto getReadable() const&
    {
      return (~*this).getRef();
    }
    
    auto getReadable() && = delete;
    
    /// Returns the size of the component
    template <typename Comp>
      constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
      Comp getCompSize() const
    {
      if constexpr(Comp::sizeIsKnownAtCompileTime)
	return Comp::sizeAtCompileTime;
      else
	return std::get<Comp>((~*this).getDynamicSizes());
    }
    
    template <typename MCL>
      constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
      auto getMergedCompsSize() const
    {
      return
	invokeWithTypesOfTuple<MCL>([this]<typename...MC>()
				    {
				      return (this->template getCompSize<MC>()()*...*1);
				    });
    }
    
    /// Returns the expression as a dynamic tensor
    auto fillDynamicTens() const
    {
      DynamicTens<typename T::Comps,typename T::Fund,T::execSpace> res((~*this).getDynamicSizes());
      
      res=~*this;
      
      return res;
    }
    
    /// Close the expression into an appropriate tensor or field
    template <typename C>
    INLINE_FUNCTION
    auto closeAs(C&& c) const
    {
      auto getStorage=
	[](auto&&...args)
	{
	  return typename std::decay_t<C>::ClosingType(args...);
	};
      
      auto res=
	std::apply(getStorage,c.getEquivalentStoragePars());
      
      res=~*this;
      
      return res;
    }
    
#define PROVIDE_CALL(ATTRIB,WHAT_TO_PASS)				\
    template <DerivedFromComp...C>					\
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE			\
      decltype(auto) operator()(const C&...cs) ATTRIB			\
    {									\
      using Comps=typename T::Comps;					\
									\
      using SubsComps=std::tuple<C...>;					\
									\
      static_assert(tupleHaveTypes<Comps,SubsComps>);			\
									\
      /*! Leftover components */					\
      using ResidualComps=						\
	TupleFilterAllTypes<Comps,SubsComps>;				\
      									\
      if constexpr(std::tuple_size_v<ResidualComps> ==0)		\
	return (constIf<not T::canAssignAtCompileTime>(~*this).eval(cs...)); \
      else								\
	return bindComps(WHAT_TO_PASS,std::make_tuple(cs...));		\
    }
    
    PROVIDE_CALL(const&,~*this);
    
    PROVIDE_CALL(/* non const */&,~*this);
    
    PROVIDE_CALL(/* non const */&&,std::move(~*this));
    
#undef PROVIDE_CALL
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_SUBSCRIBE(ATTRIB,WHAT_TO_PASS)				\
    template <DerivedFromComp C>					\
      constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE			\
      decltype(auto) operator[](const C& c) ATTRIB			\
    {									\
      return (WHAT_TO_PASS)(c);						\
    }
    
    PROVIDE_SUBSCRIBE(const &,~*this);
    
    PROVIDE_SUBSCRIBE(/* non const */&,~*this);
    
    PROVIDE_SUBSCRIBE(/* non const */&&,std::move(~*this));
    
#undef PROVIDE_SUBSCRIBE
    
    /// Computes the squared norm
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    auto norm2() const
    {
      std::decay_t<typename T::Fund> s2=0.0;
      
      compsLoop<typename T::Comps>([t=this->getReadable(),
				    &s2] // CUDA_DEVICE // seems to be making conflict...
				   (const auto&...c) INLINE_ATTRIBUTE
      {
	s2+=sqr(t(c...));
      },(~*this).getDynamicSizes());
      
      return s2;
    }
  };
  
#define PROVIDE_CATCH_OP_WITH_ARITHMETIC(OP)				\
									\
  /* Catch operator OP between a node and an arithmetic type */		\
  template <DerivedFromNode T,						\
	    typename Oth>						\
  requires(std::is_arithmetic_v<Oth>)					\
  constexpr INLINE_FUNCTION						\
  auto operator OP(const T& node,					\
		   const Oth& value)					\
  {									\
    return node OP scalar(value);					\
  }									\
  									\
  /* Catch operator OP between an arithmetic type and a node */		\
  template <DerivedFromNode T,						\
	    typename Oth>						\
  requires(std::is_arithmetic_v<Oth>)					\
  constexpr INLINE_FUNCTION						\
  auto operator OP(const Oth& value,					\
		   const T& node)					\
  {									\
    return scalar(value) OP node;					\
  }
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(+);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(-);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(*);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(/);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(%);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(==);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(!=);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(and);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(or);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(xor);
  
#undef PROVIDE_CATCH_OP_WITH_ARITHMETIC
}

#endif
