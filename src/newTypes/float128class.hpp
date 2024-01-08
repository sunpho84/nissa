#ifndef _FLOAT128CLASS_HPP
#define _FLOAT128CLASS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file src/new_types/float128class.hpp

#include <cmath>

#include <base/debug.hpp>
#include <metaprogramming/inline.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  // Quadruple precision float
  struct Float128
  {
    /// Low prec and high prec parts
    double data[2];
    
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)		\
    /*! Takes the high prec or low prec part */		\
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION	\
    CONST double& operator[](const int i) CONST		\
    {							\
      return data[i];					\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* non const */);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    /// Init knowing the two doubles
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
    Float128(const double& a=0.0,
	     const double& b=0.0) :
      data{a,b}
    {
    }
    
    /// Sum another Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    Float128 operator+(const Float128& b) const
    {
      const Float128& a=*this;
      
      const double t1=a[0]+b[0];
      const double e=t1-a[0];
      const double t2=((b[0]-e)+(a[0]-(t1-e)))+a[1]+b[1];
      
      const double t3=t1+t2;
      
      return {t3,t2-(t3-t1)};
    }
    
    /// Subtract a Float128 from a double
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    friend Float128 operator+(const double& a,
			      const Float128& b)
    {
      return b+a;
    }
    
    /// Sumassign another Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    Float128& operator+=(const Float128& b)
    {
      *this=*this+b;
      
      return *this;
    }
    
    /// Product with another Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    Float128 operator*(const Float128& b) const
    {
      const Float128& a=*this;
      const double split=134217729;
      
      const double cona=a[0]*split;
      const double conb=b[0]*split;
      
      const double a1=cona-(cona-a[0]);
      const double b1=conb-(conb-b[0]);
      const double a2=a[0]-a1;
      const double b2=b[0]-b1;
      
      const double c11=a[0]*b[0];
      const double c21=a2*b2+(a2*b1+(a1*b2+(a1*b1-c11)));
      
      const double c2=a[0]*b[1]+a[1]*b[0];
      
      const double t1=c11+c2;
      const double e=t1-c11;
      const double t2=a[1]*b[1]+((c2-e)+(c11-(t1-e)))+c21;
      
      const double t3=t1+t2;
      
      return {t3,t2-(t3-t1)};
    }
    
    /// Multiplies a double for a Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    friend Float128 operator*(const double& a,
			      const Float128& b)
    {
      return b*a;
    }
    
    /// Assigns the product with another Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    Float128& operator*=(const Float128& b)
    {
      *this=*this*b;
      
      return *this;
    }
    
    /// Unary minus
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    Float128 operator-() const
    {
      const Float128& a=*this;
      
      return {-a[0],-a[1]};
    }
    
    /// Subtract another Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    Float128 operator-(const Float128& b) const
    {
      return *this+(-b);
    }
    
    /// Subtract a Float128 from a double
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    friend Float128 operator-(const double& a,
			      const Float128& b)
    {
      return -b+a;
    }
    
    /// Assigns the difference with another Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    Float128& operator-=(const Float128& b)
    {
      *this=*this-b;
      
      return *this;
    }
    
    /// Divide two Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    Float128 operator/(const Float128& b) const
    {
      const double c=1.0/(b[0]+b[1]);
      const Float128& a=*this;
      const Float128 div1=a*c;
      const Float128 rem1=a-div1*b;
      const Float128 div2=rem1*c;
      const Float128 div12=div1+div2;
      const Float128 rem2=a-div12*b;
      const Float128 div3=rem2*c;
      
      return div12+div3;
    }
    
    /// Assigns the ratio with another Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    Float128& operator/=(const Float128& b)
    {
      *this=*this/b;
      
      return *this;
    }
    
    /// Compare two Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    bool operator==(const Float128& b) const
    {
      return (*this)[0]==b[0] and (*this)[1]==b[1];
    }
    
#define PROVIDE_COMPARE(OP)				\
    /*/ Compare two Float128 */				\
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION	\
    bool operator OP(const Float128& b) const		\
    {							\
      return						\
	((*this)[0]==b[0])?				\
	((*this)[1] OP b[1]):				\
	((*this)[0] OP b[0]);				\
    }
    
    PROVIDE_COMPARE(>);
    PROVIDE_COMPARE(>=);
    PROVIDE_COMPARE(<);
    PROVIDE_COMPARE(<=);
    
#undef PROVIDE_COMPARE
    
    /// Returns false if differ from b
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    bool operator!=(const Float128& b) const
    {
      return not ((*this)==b);
    }
    
    /// Sums two double returning a Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    static Float128 sum(const double& a,
			const double& b)
    {
      const double t1=a+b;
      const double e=t1-a;
      const double t2=((b-e)+(a-(t1-e)));
      
      const double t3=t1+t2;
      
      return {t3,t2-(t3-t1)};
    }
    
    /// Multiplies two double returning a Float128
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    static Float128 prod(const double& a,
			 const double& b)
    {
      constexpr double split=134217729;
      
      const double cona=a*split;
      const double conb=b*split;
      
      const double a1=cona-(cona-a);
      const double b1=conb-(conb-b);
      const double a2=a-a1;
      const double b2=b-b1;
      
      const double c11=a*b;
      const double c21=a2*b2+(a2*b1+(a1*b2+(a1*b1-c11)));
      
      const double t3=c11+c21;
      
      return {t3,c21-(t3-c11)};
    }
    
    /// Rounds down
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    double roundDown() const
    {
      if(data[1]>=0)
	return data[0];
      else
	return std::nextafter(data[0],-1e-300);
    }
    
    /// Rounds up
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    double roundUp() const
    {
      if(data[1]<=0)
	return data[0];
      else
	return std::nextafter(data[0],+1e-300);
    }
  };
  
  /// Perform a simple check on 128 bit precision
  inline void check128BitPrec()
  {
    Float128 a=1;
    a+=1e-20;
    a+=-1;
    
    double res=a.roundDown();
    if(fabs(res-1e-20)>1e-30)
      CRASH("float_128, 1+1e-20-1=%lg, difference with 1e-20: %lg",res,res-1e-20);
    VERBOSITY_LV2_MASTER_PRINTF("128 bit precision is working, 1+1e-20-1=%lg where %lg expected in double prec\n",res,1+1e-20-1);
  }

}

#endif
