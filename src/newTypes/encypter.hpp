#ifndef _ENCYPTER_HPP
#define _ENCYPTER_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file src/new_types/encrypter.hpp

#include <array>
#include <cstdint>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  /// Encrypts a number
  struct Encrypter
  {
    /// Type of the key
    using Key=
      std::array<uint64_t,5>;
    
    /// Encrypted word
    using Word=
      std::array<uint64_t,4>;
    
    /// Key used to encrypt
    Key key;
    
    /// Build a key from a word
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    static Key buildKey(const Word& word)
    {
      /// Output
      Key res;
      
      // Copy the first 4
      for(int i=0;i<4;i++)
	res[i]=word[i];
      
      // Set the fifth
      res[4]=0x1BD11BDAA9FC1A22^res[0]^res[1]^res[2]^res[3];
      
      return res;
    }
    
    /// Encrypts the input
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
    Word encrypt(Word x) const         ///< Input to encrypt
    {
      for(int j=0;j<5;j++)
	{
	  for(int i=0;i<2;i++)
	    {
	      constexpr uint64_t mk[2][2]=
		{{14,16},{25,33}};
	      
	      uint64_t& x0=x[2*i];
	      uint64_t& x1=x[(2*i+1)%4];
	      const uint64_t& rx=mk[j%2][i];
	      
	      x1+=key[(2*i+1+j)%5]+((i==1)?j:0);
	      x0+=x1+key[(2*i+j)%5];
	      x1=(x1<<rx)|(x1>>(64-rx));
	      x1^=x0;
	    }
	  
	  for(int l=0;l<3;l++)
	    for(int i=0;i<2;i++)
	      {
		constexpr uint64_t m[2][3][2]=
		  {{{52,57},{23,40},{5,37}},
		   {{46,12},{58,22},{32,32}}};
		
		uint64_t& x0=x[2*i];
		uint64_t& x1=x[(2*i+((l%2==0)?3:1))%4];
		const uint64_t& rx=m[j%2][l][i];
		
		x0+=x1;
		x1=(x1<<rx)|(x1>>(64-rx));
		x1^=x0;
	      }
	}
      
      // Increment entry 3
      x[3]+=5;
      
      return x;
    }
    
    /// Constructor
    Encrypter(const Word& partKey) :
      key(buildKey(partKey))
    {
    }
    
    Encrypter(const uint64_t& x) :
      Encrypter(Word{x})
    {
    }
  };
}

#endif
