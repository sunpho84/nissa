#ifndef _HASMETHOD_HPP
#define _HASMETHOD_HPP

/// Define a member detecter named hasMember_TAG
///
/// Example:
///
/// \code
/// DECLARE_HAS_MEMBER(ciccio);
///
/// struct fuffa
/// {
///    int ciccio();
/// };
///
/// int main()
/// {
///   bool has=hasMember_ciccio(fuffa);
///   return 0;
/// }
/// \endcode
#define DECLARE_HAS_MEMBER(TAG)						\
  namespace internal							\
  {									\
    /*! Detect if \c Type has member (variable or method) TAG */	\
    /*!                                                       */	\
    /*! Internal implementation                               */	\
    template <typename Type>						\
    struct HasMember_ ## TAG						\
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
	public Type,							\
	public Fallback							\
      {									\
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
      static No& test(Check<int Fallback::*,&U::TAG>*);			\
      									\
      /*! Forward definition of test function, taking a pointer to the */ \
      /*! type of TAG as argument, returning Yes. The instantiation    */ \
      /*! will work when the other fails, which means that Type does   */ \
      /*! have the member                                              */ \
      template <typename U>						\
      static Yes& test(...);						\
      									\
    public:								\
    									\
    /*! Result of the check, comparing the size of return type of */	\
    /*! test with the size of yes */					\
    static constexpr bool result=					\
      sizeof(test<Derived>(nullptr))==sizeof(Yes);			\
    };									\
    									\
    /*! Intemediate function to distinguish the non-class case */	\
    template <typename Type>						\
    [[ maybe_unused ]]							\
    constexpr bool hasMember_ ## TAG ## Helper()			\
    {									\
      if constexpr(std::is_class_v<Type>)				\
	return						                \
	  HasMember_ ## TAG<Type>::result;		                \
									\
      return								\
	  false;							\
    }									\
  }									\
  									\
  /*! Detect if \c Type has member (variable or method) TAG          */ \
  /*!                                                                */ \
  /*! Uses SFINAE to induce ambiguity in the detection of the member */ \
  template <typename Type>                                              \
  [[ maybe_unused ]]                                                    \
  constexpr bool hasMember_ ## TAG=                                     \
    internal::hasMember_ ## TAG ## Helper<std::remove_reference_t<Type>>()

#endif
