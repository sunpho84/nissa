#ifndef _MULTIUNSIGNEDINT_HPP
#define _MULTIUNSIGNEDINT_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file src/new_types/multiUnsignedInt.hpp

#include <array>
#include <type_traits>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  /// Provides a simple extended unsigned
  template <typename U,
	    int N>
  struct MultiUint
  {
    static_assert(std::is_unsigned_v<U>,"Type needs to be unsigned");
    
    /// Holds the state
    std::array<U,N> val;
    
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)			  \
    /*! Subscribe operator */					  \
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB		  \
    CONST U& operator[](const int& i) CONST			  \
    {								  \
      return val[i];						  \
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* not const */);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    /// Construct from value
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
    MultiUint(const U& u=0) :
      val{u}
    {
    }
    
    /// Increment of a certain amount
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
    MultiUint operator+(const U& z) const
    {
      /// Result
      MultiUint out=(*this);
      
      // Increment
      out[0]=val[0]+z;
      
      // Overflow check
      if(out[0]<=val[0])
	{
	  /// Digit to check for overflow
	  int iDigit=1;
	  
	  // Carry over
	  do out[iDigit]++;
	  while(out[iDigit++]==0 and iDigit<N);
	}
      
      return out;
    }
    
    /// Self increment
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    MultiUint& operator+=(const U& z)
    {
      return(*this)=(*this)+z;
    }
    
    /// Unitary self-increment
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    MultiUint& operator++()
    {
      return (*this)+=1;
    }
    
    /// Unitary self-increment
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    MultiUint operator++(int)
    {
      return (*this)+=1;
    }
  };
  
#define PROVIDE_BEGIN_OR_END(BE,ATTRIB)		\
  template <typename U,				\
	    int N>				\
  decltype(auto) BE(ATTRIB MultiUint<U,N>& m)	\
  {						\
    return BE(m.val);				\
  }

#define PROVIDE_BEGIN_AND_END(ATTRIB)		\
  PROVIDE_BEGIN_OR_END(begin,ATTRIB);		\
  PROVIDE_BEGIN_OR_END(end,ATTRIB)
  
  PROVIDE_BEGIN_AND_END(const);
  PROVIDE_BEGIN_AND_END(/* non const*/);
  
#undef PROVIDE_BEGIN_AND_END
#undef PROVIDE_BEGIN_OR_END
}

#endif
