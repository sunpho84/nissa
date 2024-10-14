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
#include <expr/fieldDeclaration.hpp>
#include <expr/assignDispatcher.hpp>
#include <expr/dynamicTensDeclaration.hpp>
#include <expr/nodeDeclaration.hpp>
#include <expr/stackTensDeclaration.hpp>
#include <metaprogramming/crtp.hpp>
#include <tuples/tupleHasType.hpp>
#include <tuples/tupleFilter.hpp>

namespace nissa
{
  namespace impl
  {
    /// Implements the cast of an expression to its Fund
    ///
    /// Forward declaration
    template <DerivedFromNode E,
	      typename C=typename E::Comps>
    struct _CastToFund;
    
    /// Implements the cast of an expression to its Fund
    template <DerivedFromNode E,
	      DerivedFromComp...C>
    struct _CastToFund<E,CompsList<C...>>
    {
      /// Detect possibility to cast
      static constexpr bool value=
	((C::sizeAtCompileTime==1) and ...);
      
      /// Exec the cast
      template <DerivedFromNode T>
      static constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
      decltype(auto) exec(T&& t)
	requires(value)
      {
	return t.eval((C)0 ...);
      }
    };
  }
  
  /// Check if an expression can be cast to its Fund
  template <DerivedFromNode E>
  constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
  bool exprCanBeCastToFund()
  {
    return impl::_CastToFund<E>::value;
  }
  
  /// Casts an expression to fund if possible
  template <DerivedFromNode E>
  constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
  decltype(auto) castExprToFund(E&& t)
  {
    return impl::_CastToFund<std::decay_t<E>>::exec(std::forward<E>(t));
  }
  
  /////////////////////////////////////////////////////////////////
  
  template <typename T,
	    typename C>
  INLINE_FUNCTION constexpr
  auto closeExprToStackTens(const Node<T,C>& node);
  
  /////////////////////////////////////////////////////////////////
  
#define THIS Node<T,CompsList<Ci...>>
  
  /// Node, the base type to be evaluated as an expression
  template <typename T,
	    DerivedFromComp...Ci>
  struct THIS :
    NodeFeat,
    Crtp<THIS,T>,
    MemberCompSubscribeProvider<T,Ci>...,
    DynamicCompsProvider<CompsList<Ci...>>
  {
    using This=THIS;
#undef THIS
    
    using Crtp<This,T>::operator~;
    
#define PROVIDE_AUTOMATIC_CAST_TO_FUND(ATTRIB)			\
    /*! Provide automatic cast to fund if needed. How cool! The trick
        is not to constrain the method, but to return something only
        if the cast is possible. Returning void, invalidates the
        cast */							\
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB		\
    operator decltype(auto)() ATTRIB				\
    {								\
      if constexpr(exprCanBeCastToFund<ATTRIB T>())		\
	return castExprToFund(DE_CRTPFY(ATTRIB T,this));	\
    }
    
    PROVIDE_AUTOMATIC_CAST_TO_FUND(const);
    
    PROVIDE_AUTOMATIC_CAST_TO_FUND(/* non const */);
    
#undef PROVIDE_AUTOMATIC_CAST_TO_FUND
    
#define PROVIDE_EXPLICIT_CAST_TO_FUND(REF,VAL,ATTRIB)			\
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB			\
      operator ATTRIB _Fund REF () ATTRIB					\
       requires(std::is_lvalue_reference_v<decltype(DE_CRTPFY(ATTRIB T,this).eval(((Ci)0)...))> ==VAL and \
 	       exprCanBeCastToFund<ATTRIB T>())				\
     {									\
       return castExprToFund(DE_CRTPFY(ATTRIB T,this));			\
     }
    
    // PROVIDE_EXPLICIT_CAST_TO_FUND(&,true,const);
    // PROVIDE_EXPLICIT_CAST_TO_FUND(/*&*/,false,const);
    // PROVIDE_EXPLICIT_CAST_TO_FUND(&,true,/* not const */);
    // PROVIDE_EXPLICIT_CAST_TO_FUND(/*&*/,false,/* not const */);
    
#undef PROVIDE_EXPLICIT_CAST_TO_FUND
    
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
    
#define PROVIDE_REINTERPRET_FUND(ATTRIB)			\
    /* Reinterprets the Fund, when lvalue reference */		\
    template <typename NFund>					\
      constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION		\
    ATTRIB auto& reinterpretFund() ATTRIB&			\
    {								\
      static_assert(sizeof(NFund)==sizeof(typename T::Fund),	\
		    "different sizes");				\
      								\
      return *(ATTRIB typename T::template			\
	       ReinterpretFund<NFund>*)this;			\
    }
    
    PROVIDE_REINTERPRET_FUND(const);
    
    PROVIDE_REINTERPRET_FUND(/* non const */);
    
#undef PROVIDE_REINTERPRET_FUND
    
    /// Reinterprets the Fund, when rvalue reference
    template <typename NFund>
      constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
      auto reinterpretFund() &&
    {
      static_assert(sizeof(NFund)==sizeof(typename T::Fund),
                      "different sizes");
      
      using Res=
	typename T::template ReinterpretFund<NFund>;
      
      return
	(Res)(std::move(*(Res*)this));
    }
    
    /// Gets a copy
    constexpr INLINE_FUNCTION
    auto getCopy() const
    {
      return
	(~*this).
	getSubExprs().
	applyTo([&self=~*this]<typename...A>(A&&...a)
		{
		  return self.recreateFromExprs(std::forward<A>(a)...);
		});
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Unless explicitly overloaded
    static constexpr bool canAssignAtCompileTime=false;
    
#define PROVIDE_MERGE_COMPS(ATTRIB,WHAT_TO_PASS)		\
    /*! Provides possibility to merge a list of components  */	\
      template <typename MCL>					\
	constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB		\
	decltype(auto) mergeComps() ATTRIB			\
      {								\
	return nissa::mergeComps<MCL>(WHAT_TO_PASS);		\
      }
    
    PROVIDE_MERGE_COMPS(const&,~*this);
    
    PROVIDE_MERGE_COMPS(/* non const */&,~*this);
    
    PROVIDE_MERGE_COMPS(/* non const */&&,std::move(~*this));
    
#undef PROVIDE_MERGE_COMPS
    
    /// Assert assignability
    template <DerivedFromNode U>
      INLINE_FUNCTION
      constexpr void assertCanAssign(const U& rhs)
    {
      // static_assert(T::execSpace.hasUniqueExecSpace(),"lhs must have a unique execution space");
      
      static_assert(U::execSpace.isCompatibleWith(T::execSpace),"incompatible execution space of rhs with lhs");
      
      //static_assert(tuplesContainsSameTypes<typename T::Comps,typename U::Comps>,"Cannot assign two expressions which differ for the components");
      
      static_assert(T::canAssignAtCompileTime,"Trying to assign to a non-assignable expression");
      
      auto& lhs=~*this;
      
      if constexpr(not T::canAssignAtCompileTime)
	CRASH("Trying to assign to a non-assignable expression");
      
      if constexpr(U::hasDynamicComps)
	if(lhs.getDynamicSizes()!=rhs.getDynamicSizes())
	  CRASH("Dynamic comps not agreeing");
      
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
	
	OP::dispatch((~*this)(lhsComps...),(typename T::Fund)std::apply(~u,rhsCompsTup));
      },(~*this).getDynamicSizes());
      
      return ~*this;
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Set of possible closing types
    enum class ClosingType{Fund,StackTens,DynamicTens,Field};
    
    /// Check if can be cast to Fund
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
      static bool _canCloseToFund()
    {
      return exprCanBeCastToFund<T>();
    }
    
    /// Check if can be closed to a StackTens
    INLINE_FUNCTION constexpr
      static bool _canCloseToStackTens()
    {
      return std::tuple_size_v<typename T::DynamicComps> ==0;
    }
    
    /// Check if can be closed to a Field
    INLINE_FUNCTION constexpr
      static bool _canCloseToField()
    {
      return tupleHasType<typename T::Comps,LocLxSite>;
    }
    
    /// Check if can be closed to a DynamicTens
    INLINE_FUNCTION constexpr
      static  bool _canCloseToDynamicTens()
    {
      return true;
    }
    
    /// Gets the optimal closing type
    INLINE_FUNCTION constexpr
      static ClosingType getClosingType()
    {
      if constexpr(_canCloseToFund())
	return ClosingType::Fund;
      else if constexpr(_canCloseToStackTens())
	return ClosingType::StackTens;
      else if constexpr(_canCloseToField())
	return ClosingType::Field;
      else
	return ClosingType::DynamicTens;
    }
    
    /// Closes to Fund
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    auto closeToFund() const
      requires(T::_canCloseToFund())
    {
      return (typename T::Fund)*this;
    }
    
    /// Closes to a StackTens
    INLINE_FUNCTION constexpr
      auto closeToStackTens() const
      requires(T::_canCloseToStackTens())
    {
      return closeExprToStackTens(*this);
    }
    
    /// Closes to a Field
    INLINE_FUNCTION constexpr
      auto closeToField() const
      requires(T::_canCloseToField());
    
    /// Closes to a DynamicTens
    INLINE_FUNCTION constexpr
    auto closeToDynamicTens() const
      requires(T::_canCloseToDynamicTens());
    
#define PROVIDE_CLOSE_TO(TYPE)						\
    /* Dispatch the correct close, needed as passing as templated would
       not be specializable */						\
      INLINE_FUNCTION constexpr						\
	auto _closeTo(std::integral_constant<ClosingType,		\
		      ClosingType::TYPE>) const				\
      {									\
	return closeTo ## TYPE();					\
      }
    
    PROVIDE_CLOSE_TO(Fund);
    PROVIDE_CLOSE_TO(StackTens);
    PROVIDE_CLOSE_TO(Field);
    PROVIDE_CLOSE_TO(DynamicTens);
    
#undef PROVIDE_CLOSE_TO
    
    /// Closes to the best possible thing
    INLINE_FUNCTION constexpr
    auto close() const
    {
      return _closeTo(std::integral_constant<ClosingType,getClosingType()>());
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Assert whether can run on device
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    static bool canRunOnDevice()
    {
      return T::execSpace.isCompatibleWith(execOnGPU);
    }
    
#define PROVIDE_ASSIGN_VARIATION(SYMBOL,OP)				\
    									\
    /*! Assign from another expression */				\
      template <DerivedFromNode Rhs>					\
	requires(not std::is_same_v<T,Rhs>)				\
	constexpr INLINE_FUNCTION					\
	T& operator SYMBOL(const Rhs& u)				\
      {									\
	(~*this).template assign<OP>(u);				\
									\
	return ~*this;							\
      }									\
    									\
    /*! Assign from a callable function */				\
      template <typename Rhs>						\
	requires(not DerivedFromNode<Rhs> and				\
		 std::is_invocable_v<Rhs,Ci...>)			\
	constexpr INLINE_FUNCTION					\
	T& operator SYMBOL(const Rhs& u)				\
      {									\
	(~*this).template assign<OP>(funcNodeWrapper<CompsList<Ci...>>(u,(~*this).getDynamicSizes())); \
									\
	return ~*this;							\
      }									\
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
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
    Node& operator=(const Node& oth)
    {
      (~*this).assign(~oth);
      
      return ~*this;
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
    
    /// Forbids to take a readable accessor if not l-value
    auto getReadable() && = delete;
    
    /// Returns the size of the component
    template <typename Comp>
      constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
      Comp getCompSize() const
    {
      if constexpr(Comp::sizeIsKnownAtCompileTime)
	return Comp::sizeAtCompileTime;
      else
	return std::get<Comp>((~*this).getDynamicSizes());
    }
    
    template <typename MCL>
      constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
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
      DynamicTens<typename T::Comps,typename T::Fund,T::exacSpace> res((~*this).getDynamicSizes());
      
      res=~*this;
      
      return res;
    }
    
    /// Close the expression into an appropriate tensor or field
    template <typename C>
    INLINE_FUNCTION
    auto closeAs(C&& c) const
    {
      /// Helper to get the closed type
      auto getStorage=
	[](auto&&...args)
	{
	  return typename std::decay_t<C>::ClosingType(args...);
	};
      
      /// Result
      auto res=
	std::apply(getStorage,c.getEquivalentStoragePars());
      
      res=~*this;
      
      return res;
    }
    
#define PROVIDE_CALL(ATTRIB,WHAT_TO_PASS)				\
    template <DerivedFromComp...C>					\
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB			\
      decltype(auto) operator()(const C&...cs) ATTRIB			\
      requires(SatisfiesNodeRequirements<T>)				\
    {									\
      using Comps=typename T::Comps;					\
									\
      using SubsComps=std::tuple<C...>;					\
									\
      static_assert(tupleHaveTypes<Comps,SubsComps>,			\
		    "Missing type");					\
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
      constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB			\
      decltype(auto) operator[](const C& c) ATTRIB			\
    {									\
      return (WHAT_TO_PASS)(c);						\
    }
    
    PROVIDE_SUBSCRIBE(const &,~*this);
    
    PROVIDE_SUBSCRIBE(/* non const */&,~*this);
    
    PROVIDE_SUBSCRIBE(/* non const */&&,std::move(~*this));
    
#undef PROVIDE_SUBSCRIBE
    
    /// Computes the squared norm
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    auto norm2() const
    {
      std::decay_t<typename T::Fund> s2=0.0;
      
      compsLoop<typename T::Comps>([t=this->getReadable(),
				    &s2] // DEVICE_ATTRIB // seems to be making conflict...
				   (const auto&...c) INLINE_ATTRIBUTE
      {
	s2+=sqr(t(c...));
      },(~*this).getDynamicSizes());
      
      return s2;
    }
  };
  
#define PROVIDE_CATCH_OP_WITH_ARITHMETIC(NAME,OP)			\
  									\
  namespace impl							\
  {									\
    /*! Not a node but can be cast to a node Fund */			\
    template <typename T,						\
	      typename F>						\
    concept _NotNodeButCan ## NAME ## With=				\
      requires(T&& x)							\
    {									\
      requires not isNode<T>;						\
      std::declval<F>() OP std::declval<T>();				\
      std::declval<T>() OP std::declval<F>();				\
    };									\
  };									\
  									\
  /* Catch operator OP between a node and an arithmetic type */		\
  template <DerivedFromNode T,						\
	    impl::_NotNodeButCan ## NAME ## With<typename T::Fund> Oth>	\
  constexpr INLINE_FUNCTION						\
  auto operator OP(const T& node,					\
		   const Oth& value)					\
  {									\
    return node OP scalar(value);					\
  }									\
  									\
  /* Catch operator OP between an arithmetic type and a node */		\
  template <DerivedFromNode T,						\
	    impl::_NotNodeButCan ## NAME ## With<typename T::Fund> Oth>	\
  constexpr INLINE_FUNCTION						\
  auto operator OP(const Oth& value,					\
		   const T& node)					\
  {									\
    return scalar(value) OP node;					\
  }
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(Sum,+);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(Sub,-);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(Prod,*);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(Div,/);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(Mod,%);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(Comp,==);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(Diff,!=);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(And,and);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(Or,or);
  
  PROVIDE_CATCH_OP_WITH_ARITHMETIC(Xor,xor);
  
#undef PROVIDE_CATCH_OP_WITH_ARITHMETIC
}

#endif
