#ifndef _CUSTOM_REAL_NUMB_HPP
#define _CUSTOM_REAL_NUMB_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <tuple>

namespace nissa
{
  /// A real number with custom mantissa and exponent
  template <size_t n,
	    size_t k>
  struct CustomRealNumb;
  
  /// Double number in the custom representation
  using CustomDouble=
    CustomRealNumb<52,11>;
  
  /// Float number in the custom representation
  using CustomFloat=
    CustomRealNumb<23,8>;
  
  /// Hald number in the custom representation
  using CustomHalf=
    CustomRealNumb<10,5>;
  
  template <size_t n,
	    size_t k>
  struct CustomRealNumb
  {
    // Total length of the custom float
    static constexpr size_t t=
      n+k+1;
    
    static_assert(t<=64,"Total size exceeds 64 bits");
    
    /// Number of bytes needed to host the data
    static constexpr size_t nBytes=
      (t+7)/8;
    
    /// Internal storage type
    using StorageType=
      std::tuple_element_t<nBytes-1,std::tuple<uint8_t,uint16_t,uint32_t,uint32_t,uint64_t,uint64_t,uint64_t,uint64_t>>;
    
    /// Storage for the data
    StorageType data;
    
    /// Set the sign
    constexpr void setSign(const uint8_t& sign)
    {
      data=(data&~(StorageType(1)<<(n+k)))|(StorageType(sign)<<(n+k));
    }
    
    /// Set the exponent
    constexpr void setExponent(const uint64_t& exponent)
    {
      data=(data&~(((StorageType(1)<<k)-1)<<n))|(StorageType(exponent)<<n);
    }
    
    /// Set the mantissa
    constexpr void setMantissa(const uint64_t& mantissa)
    {
      data=(data&~((StorageType(1)<<n)-1))|StorageType(mantissa);
    }
    
    /// Get the sign
    constexpr uint8_t getSign() const
    {
      return (data>>(n+k))&1;
    }
    
    /// Get the exponent
    constexpr uint64_t getExponent() const
    {
      return (data>>n)&((StorageType(1)<<k)-1);
    }
    
    // Get the mantissa
    constexpr uint64_t getMantissa() const
    {
      return data&((StorageType(1)<<n)-1);
    }
    
    /// Initialize from a double
    constexpr CustomRealNumb(const double& value=0)
    {
      if constexpr(std::is_same_v<CustomRealNumb,CustomDouble>)
	memcpy(&data,&value,sizeof(double));
      else
	{
	  /// Temporarily holds the quantity in a CustomDouble
	  CustomDouble tmp(value);
	  
	  /// Calculate new exponent
	  int newExponent=
	    tmp.getExponent()-1023+((1<<(k-1))-1);
	  
	  /// Shift mantissa to fit in the custom mantissa size
	  uint64_t newMantissa=
	    tmp.getMantissa()>>(52-n);
	  
	  // Handle underflow and overflow
	  if(newExponent<0)
	    {
	      newExponent=0;
	      newMantissa=0;
	    }
	  else
	    if(newExponent>=(1<<k)-1)
	      {
		newExponent=(1<<k)-1;
		newMantissa=(1ULL<<n)-1;
	      }
	  
	  // Set the custom float components
	  setMantissa(newMantissa);
	  setExponent(newExponent);
	  setSign(tmp.getSign());
	}
    }
    
    /// Cast to double
    explicit constexpr operator double() const
    {
      // Combine the parts to form the double value
      const uint64_t bits=
	((uint64_t)getSign()<<63)|((uint64_t)(getExponent()-((1<<(k-1))-1)+1023)<<52)|(getMantissa()<<(52-n));
      
      return *reinterpret_cast<const double*>(&bits);
    }
  };
}

#endif
