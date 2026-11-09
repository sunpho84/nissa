#ifndef _SITMO_HPP
#define _SITMO_HPP

// Single implementation shared by existing consumers and NICE
// See doc/sitmo.md for numerical compatibility and ongoing runs

#include <cstdint>

#include <base/field.hpp>
#include <geometry/geometry_lx.hpp>
#include <threads/threads.hpp>
#include <new_types/arr.hpp>
#include <routines/math_routines.hpp>

/*

The SitmoRng produces a stream of uint32, which must be rectified into
a stream of double. The stream is associated to lattice sites, such
that the i-th random number of the j-th site is the k-th double
number, with k=j+glb_vol*i. In this way, to draw n random numbers r_i
for each lattice site, one needs to move with glb_vol stride along the
stream. It this convenient to define a structure FieldRngStream which
wraps the SitmoRng, providing a view on the stream, in which the
number of random numbers drawn is kept homoegeneous across the
lattice. The user must draw N numbers in turns, requiring the kind of
local structure T to be be filled. The request produces a proxy,
FieldRngOf<T>, which keeps track of the request. When the request is
made, the status of the FieldRngStream is advanced. The user can pick
up the requested numbers from the FieldRngOf<T>, requiring a specific
site. Only at this moment the actual Sitmo is accessed, and the
numbers generated.

SitmoRng
--------
-initialized on a given seed
-given an uint64, it produces an uint32 in the range 0-max

FieldRngStream
--------------
-embeds the SitmoRng
-keeps track of the number of real random real generated across the
lattice, which in turns correspond to a certain number of calls to the
SitmoRng
-accessed requiring to draw a certain number of double
-this call returns a proxy and advances the number of random reals

FieldRngOf<T>
--------------------
-keeps reference to the Sitmo of the FieldRngStream
-when passed a T, and the local lattice site, it fill T with the required elements
 */

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
  //Embeds the Sitmo rng at a certain point in the stream
  struct RngView
  {
    /// Field random number generator
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
  
  /// Consume two 32-bit draws, low significance first, and clamp rounded unity
  template <typename G>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  double uniformDouble(G&& generator)
  {
    constexpr double F=4294967296.0;
    
    const uint32_t low=generator();
    const uint32_t high=generator();
    
    if(const double result=(low/F+high)/F;result>=1)
      return 0.999999999999999889;
    else
      return result;
  }
  
  template <typename T>
  struct FieldRngOf
  {
    //Mark if the generator can be used or not
    bool used;
    
    //Reference to the Sitmo
    Sitmo::Rng rng;
    
    //Number of reals to be shifted, in units of global volume
    const Sitmo::Counter offsetReal;
    
    //Number of reals for each site
    static constexpr int nRealsPerSite=
      sizeof(T)/sizeof(double);
    
    //Constructor
    FieldRngOf(const Sitmo::Rng& rng,
	       const uint64_t& offsetReal) :
      FieldRngOf(rng,Sitmo::Counter{{{offsetReal,0,0,0}}})
    {
    }
    
    /// Retain the full-width per-site double offset for this drawer
    FieldRngOf(const Sitmo::Rng& rng,
               const Sitmo::Counter& offsetReal) :
      used(false),
      rng(rng),
      offsetReal(offsetReal)
    {
    }
    
    /// Returns a view on a specific site and real number
    CUDA_HOST_AND_DEVICE
    RngView getRngViewOnGlbSiteIRndReal(const uint64_t glblx,
                                        const unsigned int iRndRealPerSite)
    {
      // Require a valid site in a nonempty geometry before constructing its offset
      if(nissa::glbVol<=0 or glblx>=static_cast<uint64_t>(nissa::glbVol))
        CRASH("RNG global site index out of range");
      
      const uint64_t volume=static_cast<uint64_t>(nissa::glbVol);
      
      // Preserve every word when flattening the site-fastest field position
      const Sitmo::Counter slot=offsetReal+uint64_t{iRndRealPerSite};
      const Sitmo::Counter irndDouble=slot*volume+glblx;
      
      // Each uniform double consumes two consecutive 32-bit generator outputs
      const Sitmo::Counter irndUint32=irndDouble*uint64_t{2};
      
      return RngView(rng,irndUint32);
    }
    
    
    /// Gets a number between [0;1)
    ///
    /// Consume two 32-bit words, low significance first, and clamp rounded unity
    template <typename G>
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    double getUniformDistributed(G&& g)
    {
      return uniformDouble(g);
    }
    
    /// Fill a specific site
    CUDA_HOST_AND_DEVICE
    void fillGlbSite(T& out,
		     const uint64_t& glblx)
    {
      for(int irnd_real=0;irnd_real<nRealsPerSite;irnd_real++)
	{
	  auto view=
	    getRngViewOnGlbSiteIRndReal(glblx,irnd_real);
	  
	  ((double*)out)[irnd_real]=getUniformDistributed(view);
	}
    }
    
    /// Fill a specific site given its local index
    CUDA_HOST_AND_DEVICE
    void fillLocSite(T& out,
		     const uint64_t& loclx)
    {
      fillGlbSite(out,glblxOfLoclx[loclx]);
    }
    
    /// Claims that the generator must be used only once
    void enforce_single_usage()
    {
      if(used)
	CRASH("cannot use a rng field twice");
      used=true;
    }
    
    /// Fill with origin
    void fillWithOrigin(T& out)
    {
      enforce_single_usage();
      
      fillGlbSite(out,0);
    }
    
    /// Fill all sites
    void fillField(LxField<T>& out)
    {
      enforce_single_usage();
      
      PAR(0,
	  locVol,
	  CAPTURE(self=*this,
		  TO_WRITE(out)),
	  site,
	{
	  for(int iDeg=0;
              iDeg<std::remove_reference_t<decltype(out)>::nInternalDegs;
	      iDeg++)
	    {
              auto view=
		self.getRngViewOnGlbSiteIRndReal(glblxOfLoclx[site],iDeg);
              out(site,iDeg)=self.getUniformDistributed(view);
            }
	});
    }
  };
  
  /// Stream of numbers for all lattice sites
  ///
  /// This is the object that is accessed from external routines
  struct FieldRngStream
  {
    /// Take note of initialization
    bool inited{};
    
    /// Embeds the random number generator
    Sitmo::Rng rng{};
    
    //Number of double precision numbers generated per site
    Sitmo::Counter nGeneratedDouble{};
    
    /// Skip n drawers
    template <typename T>
    void skipDrawers(const uint64_t& nDrawers)
    {
      constexpr uint64_t nReals=FieldRngOf<T>::nRealsPerSite;
      
      if constexpr(nReals>0)
        {
          // Promote before multiplying so large skips retain their high words
          nGeneratedDouble+=Sitmo::Counter{{{nReals,0,0,0}}}*nDrawers;
        }
    }
    
    /// Return a proxy from which to fetch the random numbers
    template <typename T>
    auto getDrawer()
    {
      static_assert(sizeof(T)%sizeof(double)==0,"Type not multiple in size of double");
      
      /// Result
      FieldRngOf<T> res(rng,nGeneratedDouble);
      
      skipDrawers<T>(1);
      
      return res;
    }
    
    /// Draw a single instance of T
    template <typename T>
    void drawScalar(T& out)
    {
      static_assert(sizeof(T)%sizeof(double)==0,"Type not multiple in size of double");
      
      auto drawer=getDrawer<T>();
      
      drawer.fillWithOrigin(out);
    }
    
    /// Initializes with a seed
    void init(const uint32_t& seed)
    {
      if(inited)
	CRASH("already inited");
      inited=true;
      
      rng.seed(seed);
    }
    
    /// Returns the commands to exec hen starting/ending the benchmark
    bool maybeBackupForBenchmark()
    {
#ifdef USE_CUDA
      if(insideParallelFor)
	{
	  benchmarkEndActions.emplace_back([this,
					    n=nGeneratedDouble]()
	  {
	    VERBOSITY_LV3_MASTER_PRINTF("Restoring FieldRngStream\n");
	    this->nGeneratedDouble=n;
	  });
	
	VERBOSITY_LV3_MASTER_PRINTF("Added backup action for FieldRngStream\n");
	
	return true;
      }
    else
      {
	VERBOSITY_LV3_MASTER_PRINTF("Skipping backup actions for FieldRngStream as we are not inside a parallel for\n");
	
	return false;
      }
#else
      return false;
#endif
    }
  };
  
  static constexpr complex zero_complex={0.0,0.0};
  static constexpr complex sqrt_2_half_complex{M_SQRT1_2,M_SQRT1_2};
  
  template <typename C>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  void BoxMullerTransform(C&& out,
				 const complex& ave=zero_complex,
				 const complex& sig=sqrt_2_half_complex)
  {
    const double r=sqrt(-2*log(1-out[RE]));
    const double q=2*M_PI*out[IM];
    
    out[RE]=r*cos(q)*sig[RE]+ave[RE];
    out[IM]=r*sin(q)*sig[IM]+ave[IM];
  }
  
  /// Remap a number in the range [0;1) to {-1;+1}
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  double remapUniformToSign(const double& in)
  {
    return int(in*2.0)*2-1;
  }
  
  template <typename C>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  void z2Transform(C&& out)
  {
    out[RE]=remapUniformToSign(out[RE]);
    out[IM]=0.0;
  }
  
  template <typename C>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  void z4Transform(C&& out)
  {
    UNROLL_FOR_RI(ri)
      out[ri]=remapUniformToSign(out[ri])*M_SQRT1_2;
  }
}

#endif
