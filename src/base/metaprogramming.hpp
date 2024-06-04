#ifndef _METAPROGRAMMING_HPP
#define _METAPROGRAMMING_HPP

#include <tuple>
#include <type_traits>
#include <utility>

#include <routines/math_routines.hpp>

namespace nissa
{
  
#if __cplusplus > 201402L
  
  using std::as_const;
  
#else
  
  template <typename T>
  constexpr const T& as_const(T& t) noexcept
  {
    return t;
  }
  
#endif
  
  /// Provides a SFINAE to be used in template par list
  ///
  /// This follows
  /// https://stackoverflow.com/questions/32636275/sfinae-with-variadic-templates
  /// as in this example
  /// \code
  /// template <typename D,
  ///           SFINAE_ON_TEMPLATE_ARG(IsSame<D,int>)>
  /// void foo(D i) {} // fails if D is not int
  /// \endcode
#define SFINAE_ON_TEMPLATE_ARG(...)	\
  std::enable_if_t<(__VA_ARGS__),void*> =nullptr
  
  /// Returns true if T is a const lvalue reference
  template <typename T>
  constexpr bool is_const_lvalue_reference_v=std::is_lvalue_reference<T>::value and std::is_const<std::remove_reference_t<T>>::value;
  
  /// Returns the type without "const" attribute if it is a reference
  template <typename T>
  CUDA_HOST_AND_DEVICE
  decltype(auto) remove_const_if_ref(T&& t)
  {
    using Tv=std::remove_const_t<std::remove_reference_t<T>>;
    
    return (std::conditional_t<is_const_lvalue_reference_v<T>,Tv&,Tv>)t;
  }
  
  /// If the type is an l-value reference, provide the type T&, otherwise wih T
  template <typename T>
  using ref_or_val_t=std::conditional_t<std::is_lvalue_reference<T>::value,T&,T>;
  
  /// Provides also a non-const version of the method \c NAME
  ///
  /// See
  /// https://stackoverflow.com/questions/123758/how-do-i-remove-code-duplication-between-similar-const-and-non-const-member-func
  /// A const method NAME must be already present Example
  ///
  /// \code
  // class ciccio
  /// {
  ///   double e{0};
  ///
  /// public:
  ///
  ///   const double& get() const
  ///   {
  ///     return e;
  ///   }
  ///
  ///   PROVIDE_ALSO_NON_CONST_METHOD(get);
  /// };
  /// \endcode
#define PROVIDE_ALSO_NON_CONST_METHOD(NAME)				\
  /*! Overload the \c NAME const method passing all args             */ \
  template <typename...Ts> /* Type of all arguments                  */	\
  decltype(auto) NAME(Ts&&...ts) /*!< Arguments                      */ \
  {									\
    return remove_const_if_ref(as_const(*this).NAME(std::forward<Ts>(ts)...)); \
  }
  
#define PROVIDE_ALSO_NON_CONST_METHOD_GPU(NAME)				\
  /*! Overload the \c NAME const method passing all args             */ \
  template <typename...Ts> /* Type of all arguments                  */	\
  CUDA_HOST_AND_DEVICE decltype(auto) NAME(Ts&&...ts) /*!< Arguments     */ \
  {									\
    return remove_const_if_ref(as_const(*this).NAME(std::forward<Ts>(ts)...)); \
  }
  
  /// Implements the CRTP pattern
  template <typename T>
  struct Crtp
  {
    /// Crtp access the type
    const T& crtp() const
    {
      return *static_cast<const T*>(this);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD(crtp);
  };
  
  /// Filter a tuple on the basis of a predicate on the type
  ///
  /// Internal implementation working out a single type, forward
  /// declaration
  template <bool,
	    typename>
  struct _TupleFilter;
  
  /// Filter a tuple on the basis of a predicate
  ///
  /// True case, in which the type is filtered
  template <typename T>
  struct _TupleFilter<true,T>
  {
    /// Helper type, used to cat the results
    using type=std::tuple<T>;
    
    /// Filtered value
    const type value;
    
    /// Construct, taking a tuple type and filtering the valid casis
    template <typename Tp>
    _TupleFilter(Tp&& t) : ///< Tuple to filter
      value{std::get<T>(t)}
    {
    }
  };
  
  /// Filter a tuple on the basis of a predicate
  ///
  /// True case, in which the type is filtered out
  template <typename T>
  struct _TupleFilter<false,T>
  {
    /// Helper empty type, used to cat the results
    using type=std::tuple<>;
    
    /// Empty value
    const type value{};
    
    /// Construct without getting the type
    template <typename Tp>
    _TupleFilter(Tp&& t) ///< Tuple to filter
    {
    }
  };
  
  /// Returns a tuple, filtering out the non needed types
  template <template <class> class F,          // Predicate to be applied on the types
	    typename...T>                      // Types contained in the tuple to be filtered
  auto tupleFilter(const std::tuple<T...>& t) ///< Tuple to filter
  {
    return std::tuple_cat(_TupleFilter<F<T>::value,T>{t}.value...);
  }
  
  /// Type obtained applying the predicate filter F on the tuple T
  template <template <class> class F,
	    typename T>
  using TupleFilter=decltype(tupleFilter<F>(T{}));
  
  /////////////////////////////////////////////////////////////////
  
  /// Directly provides the result of filtering out from a tuple
  template <typename F,
	    typename Tp>
  struct _TupleFilterOut;
  
  /// Cannot use directly the TupleFilter, because of some template template limitation
  template <typename...Fs,
	    typename...Tps>
  struct _TupleFilterOut<std::tuple<Fs...>,std::tuple<Tps...>>
  {
    /// Predicate to filter out
    template <typename T>
    struct filter
    {
      /// Predicate result, counting whether the type match
      static constexpr bool value=(sumAll<int>(std::is_same<T,Fs>::value...)==0);
    };
    
    /// Returned type
    typedef TupleFilter<filter,std::tuple<Tps...>> type;
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Directly provides the result of filtering out the types of the tuple F from Tuple Tp
  template <typename F,
	    typename Tp>
  using TupleFilterOut=typename _TupleFilterOut<F,Tp>::type;
  
  /// Predicate returning whether the type is present in the list
  ///
  /// Forward definition
  template <int N,
	    typename Tp>
  struct TypeIsInList;
  
  /// Predicate returning whether the type is present in the list
  template <int N,
	    typename...Tp>
  struct TypeIsInList<N,std::tuple<Tp...>>
  {
    /// Internal implementation
    template <typename T>
    struct t
    {
      /// Predicate result
      static constexpr bool value=(sumAll<int>(std::is_same<T,Tp>::value...)==N);
    };
  };
  
  /////////////////////////////////////////////////////////////////
  
    /// Define a member detecter named hasMember_TAG
  ///
  /// Example:
  ///
  /// \code
  /// DEFINE_HAS_MEMBER(mace);
  ///
  /// struct fuffa
  /// {
  ///    int mace();
  /// };
  ///
  /// int main()
  /// {
  ///   bool has=hasMember_mace(fuffa);
  ///   return 0;
  /// }
  /// \endcode
#define PROVIDE_HAS_MEMBER(TAG)						\
  namespace impl							\
  {									\
    /*! Detect if \c Type has member (variable or method) TAG */	\
    /*!                                                       */	\
    /*! Forward definition                                    */	\
    template <typename Type,						\
	      bool IsClass>						\
      struct HasMember_ ## TAG;						\
    									\
    /*! Detect if \c Type has member (variable or method) TAG */	\
    /*!                                                       */	\
    /*! Internal implementation for class                     */	\
    template <typename Type>						\
      struct HasMember_ ## TAG<Type,true>				\
    {									\
      /*! Internal class of size 1, used if Type has the method */	\
      using Yes=char[1];						\
      									\
      /*! Internal class of size 2 used if Type has not the method */	\
      using No=char[2];							\
      									\
      /*! Internal class which does implement the method TAG */		\
      struct Fallback							\
      {									\
	/*! Member actually implemented */				\
	int TAG;							\
      };								\
      									\
      /*! This class inherits from both Type and Fallback, so it will  */ \
      /*! certainly have the method TAG, possibly two                  */ \
      struct Derived :							\
	      public Type,						\
	      public Fallback						\
	{								\
	};								\
      									\
      /*! This type can be instantiated only if the U type is          */ \
      /*! unambiguosly understood.*/					\
      template <typename U,						\
		U>							\
      struct Check;							\
      									\
      /*! Forward definition of test function, taking a pointer to the */ \
      /*! type of TAG as argument, returning No. The instantiation     */ \
      /*! will fail if Base have two member TAG implemented, which     */ \
      /*! means that Type has the member                               */ \
      template <typename U>						\
	static No& test(Check<int Fallback::*,&U::TAG>*);		\
      									\
      /*! Forward definition of test function, taking a pointer to the */ \
      /*! type of TAG as argument, returning Yes. The instantiation   */ \
      /*! will work when the other fails, which means that Type does  */ \
      /*! have the member                                             */ \
      template <typename U>						\
	static Yes& test(...);						\
      									\
    public:								\
      /*! Result of the check, comparing the size of return type of */	\
      /*! test with the size of yes */					\
      static constexpr bool value=					\
	sizeof(test<Derived>(nullptr))==sizeof(Yes);			\
    };									\
  									\
    /*! Detect if \c Type has member (variable or method) TAG */	\
    /*!                                                       */	\
    /*! Internal implementation for not class                 */	\
    template <typename Type>						\
      struct HasMember_ ## TAG<Type,false>				\
    {									\
    public:								\
      /*! Result of the check, always false */				\
      static constexpr bool value=					\
	false;								\
    };									\
  }									\
									\
  /*! Detect if \c Type has member (variable or method) TAG          */	\
  template <typename Type>						\
  [[ maybe_unused ]]							\
  constexpr bool hasMember_ ## TAG=					\
    impl::HasMember_ ## TAG<std::decay_t<Type>,				\
			    std::is_class<std::decay_t<Type>>::value>::value
  
  /////////////////////////////////////////////////////////////////
  
  /// Returns a tuple containing all types common to the two tuples
  template <typename TupleToSearch,
	    typename TupleBeingSearched>
  using TupleCommonTypes=TupleFilter<TypeIsInList<1,TupleToSearch>::template t,TupleBeingSearched>;
  
  namespace Robbery
  {
    /// Returns the shadow pointer
    template <typename T,
	      const int Tag>
    struct Shadower
    {
      typedef T* T::*type;
      
#if (defined(__GNUC__) && !defined(__clang__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-template-friend"
#endif
      friend auto get(Shadower<T,Tag>);
#if (defined(__GNUC__) && !defined(__clang__))
#pragma GCC diagnostic pop
#endif
    };
    
    /// Performs robbery
    template<const int TagName,
	     typename T,
	     auto M>
    struct Rob
    {
      friend auto get(Shadower<T,TagName>)
      {
	return M;
      }
    };
    
    /// Performs the robbery on the passed pointer
    template <auto Tag,
	      typename P>
    decltype(auto) rob(P* p)
    {
      return p->*get(Shadower<P,Tag>());
    }
  }
}

#endif
