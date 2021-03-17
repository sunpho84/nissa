#ifndef _VALUE_WITH_EXTREME_HPP
#define _VALUE_WITH_EXTREME_HPP

#include <limits>

namespace nissa
{
  /// Possible extreme types
  enum Extreme{MINIMUM,MAXIMUM};
  
  /// Class which keeps track of extreme values of a given type
  template <typename T,
	    Extreme E>
  class ValWithExtreme
  {
    /// Stored value
    T val;
    
    /// Extreme value
    T extr;
    
    /// Update the extreme
    ValWithExtreme& updateExtreme()
    {
      /// Result of whether it's extreme or not
      bool is;
      
      switch(E)
	{
	case MINIMUM:
	  is=(val<extr);
	  break;
	case MAXIMUM:
	  is=(val>extr);
	  break;
	}
      
      if(is)
	extr=val;
      
      return *this;
    }
    
  public:
    
    /// Retrurn extreme value
    const T& extreme() const
    {
      return extr;
    }
    
    /// Reset to standard value
    template <typename V=T>
    void reset(const V& init)
    {
      switch(E)
	{
	case MINIMUM:
	  extr=val=init;
	  break;
	case MAXIMUM:
	  extr=val=init;
	  break;
	}
    }
    
    /// Constructor
    template <typename V=T>
    ValWithExtreme(const V& init)
    {
      reset(init);
    }
    
    /// Implicit cast to const value reference
    operator const T&() const
    {
      return val;
    }
    
    /// Provide an unary operator \c OP
#define PROVIDE_UNARY_OPERATOR(OP)		\
    /*! Unary operator \c OP with update */	\
    template <typename V>			\
    ValWithExtreme& operator OP (const V& oth)	\
    {						\
      val OP oth;				\
      						\
      return updateExtreme();			\
    }
    
    PROVIDE_UNARY_OPERATOR(+=);
    
    PROVIDE_UNARY_OPERATOR(-=);
    
    PROVIDE_UNARY_OPERATOR(=);
    
#undef PROVIDE_UNARY_OPERATOR
  };
  
  /// class to keep a value and its maximum
  template <typename T>
  using ValWithMax=ValWithExtreme<T,MAXIMUM>;
}

#endif
