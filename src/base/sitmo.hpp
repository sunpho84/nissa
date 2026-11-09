#ifndef _SITMO_HPP
#define _SITMO_HPP

// Single engine implementation shared by existing consumers and NICE
// See doc/sitmo.md for numerical compatibility and ongoing runs

#include <cstdint>

#include <base/debug.hpp>
#include <new_types/arr.hpp>
#include <routines/math_routines.hpp>

namespace Sitmo
{
  /// Type of the key
  using Key=
    nissa::Arr<uint64_t,5>;
  
  /// Four-word counter with the least significant word first, also used for cipher output
  struct Counter :
    nissa::Arr<uint64_t,4>
  {
    /// Add a full-width increment and reject carry beyond the fourth word
    CUDA_HOST_AND_DEVICE
    Counter operator+(const Counter& increment) const
    {
      Counter out{};
      uint64_t carry=0;
      
      for(int i=0;i<4;i++)
        {
          // Retain both the limb sum's carry and the carry from adding the prior limb
          const auto [sum,sumCarry]=nissa::addWithCarry((*this)[i],increment[i]);
          const auto [result,incomingCarry]=nissa::addWithCarry(sum,carry);
          out[i]=result;
          carry=sumCarry|incomingCarry;
        }
      
      // The complete counter cannot wrap into an already used position
      if(carry)
        CRASH("RNG counter overflow");
      
      return out;
    }
    
    /// Add a single-word increment using the same carry and overflow rules
    CUDA_HOST_AND_DEVICE
    Counter operator+(const uint64_t increment) const
    {
      return (*this)+Counter{{{increment,0,0,0}}};
    }
    
    /// Multiply by a single-word factor, retaining high words and checking full overflow
    CUDA_HOST_AND_DEVICE
    Counter operator*(const uint64_t factor) const
    {
      Counter out{};
      uint64_t carry=0;
      
      for(int i=0;i<4;i++)
        {
          const auto [low,high]=nissa::multiplyWithCarry((*this)[i],factor);
          
          // The product's high word is at most UINT64_MAX-1, leaving room for this carry
          const auto [result,sumCarry]=nissa::addWithCarry(low,carry);
          out[i]=result;
          carry=high+sumCarry;
        }
      
      // Reject a product beyond the complete four-word counter
      if(carry!=0)
        CRASH("RNG counter overflow");
      
      return out;
    }
    
    /// Scale in place without narrowing the product to one word
    CUDA_HOST_AND_DEVICE
    Counter& operator*=(const uint64_t factor)
    {
      return (*this)=(*this)*factor;
    }
    
    /// Advance in place by a full-width increment
    CUDA_HOST_AND_DEVICE
    Counter& operator+=(const Counter& increment)
    {
      return (*this)=(*this)+increment;
    }
    
    /// Advance in place by a single-word increment
    CUDA_HOST_AND_DEVICE
    Counter& operator+=(const uint64_t increment)
    {
      return (*this)=(*this)+increment;
    }
    
    /// Advance by one and return the updated counter
    CUDA_HOST_AND_DEVICE
    Counter& operator++()
    {
      return (*this)+=uint64_t{1};
    }
  };
  
  /// Encrypts the input
  CUDA_HOST_AND_DEVICE
  inline Counter encrypt(const Key& key,  ///< Key to encrypt
			 Counter x)       ///< Input to encrypt
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
    
    // Final key injection
    for(int i=0;i<4;i++)
      x[i]+=key[i];
    
    // Increment entry 3
    x[3]+=5;
    
    return x;
  }
  
  /// Build a key from a four-word counter
  inline Key buildKey(const Counter& word)
  {
    /// Output
    Key key;
    
    // Copy the first 4
    for(int i=0;i<4;i++)
      key[i]=word[i];
    
    // Set the fifth
    key[4]=0x1BD11BDAA9FC1A22^key[0]^key[1]^key[2]^key[3];
    
    return key;
  }
  
  /// Class encapsulating the Sitmo random generator
  struct Rng
  {
    /// Default seed
    static constexpr uint32_t DEFAULT_SEED()
    {
      return 3472291050;
    }
    
    /// Counter position retained under the existing state type name
    using State=Counter;
    
    /// Seed, used to encrypt
    Key key;
    
    /// State of the random number generator
    State state{};
    
    /// Type returned
    using ResultType=uint32_t;
    
    static constexpr ResultType max_value=(~(ResultType)0);
    // /// Draw a seed from a random generator
    // template <typename F>
    // void seed(F& f)
    // {
    //   /// Result type of the callabale generator
    //   using FRes=typename F::ResultType;
      
    //   /// Number of calls to be issued to fill the state
    //   constexpr int nCall=sizeof(Counter)/sizeof(FRes);
      
    //   union
    //   {
    // 	/// Partial key
    // 	Counter partKey;
	
    // 	/// Access to the key in a way which allows to fill with the generator
    // 	FRes rawKey[nCall];
    //   };
      
    //   /// Fill the key
    //   for(int i=0;i<nCall;i++)
    // 	rawKey[i]=f();
      
    //   key=buildKey(partKey);
    // }
    
    /// Set the seed and reset the counter to its initial position
    void seed(const uint32_t& s)
    {
      /// Partial word to be used
      const Counter& partKey={s};
      
      key=buildKey(partKey);
      state={};
    }
    
    /// Generate at a full-width offset without changing the retained state
    CUDA_HOST_AND_DEVICE
    ResultType generateNth(const Counter& offset) const
    {
      /// Shifted state
      const State shiftedState=state+offset;
      
      /// Cipher output whose first word supplies the low 32 bits
      const Counter encrypted=encrypt(key,shiftedState);
      
      return static_cast<ResultType>(encrypted[0]);
    }
    
    /// Generate at a single-word offset using the full-width positioning rules
    CUDA_HOST_AND_DEVICE
    ResultType generateNth(const uint64_t offset) const
    {
      return generateNth(Counter{{{offset,0,0,0}}});
    }
    
    /// Default constructor
    explicit Rng(const uint32_t& s=DEFAULT_SEED())
    {
      seed(s);
    }
  };
}

namespace nissa
{
  /// Cursor over a copied Sitmo engine, independent of lattice or distribution mapping
  struct RngView
  {
    /// Copied engine snapshot
    Sitmo::Rng ref;
    
    /// Next integer draw offset, retained at its maximum after exhaustion
    uint64_t irng;
    
    /// Whether the last representable view position has already been consumed
    bool exhausted{};
    
    //Minimal value
    constexpr uint32_t min()
    {
      return 0;
    }
    
    //Maximal value
    constexpr uint32_t max()
    {
      return ~0;
    }
    
    //Constructor
    CUDA_HOST_AND_DEVICE
    RngView(Sitmo::Rng& ref,
	    const uint64_t& irng) :
      ref(ref),
      irng(irng)
    {
    }
    
    /// Position a view at a full-width offset without changing its parent
    CUDA_HOST_AND_DEVICE
    RngView(const Sitmo::Rng& ref,
            const Sitmo::Counter& offset) :
      ref(ref),
      irng(0)
    {
      this->ref.state+=offset;
    }
    
    /// Returns an uint32
    CUDA_HOST_AND_DEVICE
    uint32_t operator()()
    {
      // Reject another draw after the maximum view offset has been consumed
      if(exhausted)
        CRASH("RNG view position overflow");
      
      const uint32_t result=ref.generateNth(irng);
      
      // Consume the final valid offset once without wrapping the next position
      if(irng==UINT64_MAX)
        exhausted=true;
      else
        irng++;
      
      return result;
    }
  };
}

#endif
