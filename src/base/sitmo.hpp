#ifndef _SITMO_HPP
#define _SITMO_HPP

#include <cstdint>
#include <random>

#include <base/vectors.hpp>
#include <geometry/geometry_lx.hpp>
#include <threads/threads.hpp>

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

namespace nissa
{
  /// Simple array to overcome limitations of std::array
  template <typename T,
	    int N>
  struct Arr
  {
    /// Internal data
    T data[N];
    
    /// Access to internal data, constant attribute
    CUDA_HOST_AND_DEVICE
    const T& operator[](const int& i) const
    {
      return data[i];
    }
  
    /// Access to internal data
    CUDA_HOST_AND_DEVICE
    T& operator[](const int& i)
    {
      return data[i];
    }
  };
}

namespace Sitmo
{
  /// Type of the key
  using Key=
    nissa::Arr<uint64_t,5>;
  
  /// Encrypted word
  using Word=
    nissa::Arr<uint64_t,4>;
  
  /// Encrypts the input
  CUDA_HOST_AND_DEVICE
  inline Word encrypt(const Key& key,  ///< Key to encrypt
		      Word x)          ///< Input to encrypt
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
  
  /// Build a key from a word
  inline Key buildKey(const Word& word)
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
    
    /// Holds the state of the generator
    struct State : public Sitmo::Word
    {
      /// Increment of a certain amount
      CUDA_HOST_AND_DEVICE
      State operator+(const uint64_t& z) const
      {
	/// Result
	State out=(*this);
	
	// Increment
	out[0]=(*this)[0]+z;
	
	// Overflow check
	if(out[0]<=(*this)[0])
	  {
	    /// Digit to check for overflow
	    int iDigit=1;
	    
	    // Carry over
	    do out[iDigit]++;
	    while(out[iDigit++]==0 and iDigit<4);
	  }
	
	return out;
      }
      
      /// Self increment
      State operator+=(const uint64_t& z)
      {
	return(*this)=(*this)+z;
      }
      
      /// Unitary self-increment
      State operator++()
      {
	return (*this)+=1;
      }
    };
    
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
    //   constexpr int nCall=sizeof(Word)/sizeof(FRes);
      
    //   union
    //   {
    // 	/// Partial key
    // 	Word partKey;
	
    // 	/// Access to the key in a way which allows to fill with the generator
    // 	FRes rawKey[nCall];
    //   };
      
    //   /// Fill the key
    //   for(int i=0;i<nCall;i++)
    // 	rawKey[i]=f();
      
    //   key=buildKey(partKey);
    // }
    
    /// Seed from a seed
    void seed(const uint32_t& s)
    {
      /// Partial word to be used
      const Word& partKey={s};
      
      key=buildKey(partKey);
    }
    
    /// Generate a number with a given offset w.r.t current state
    CUDA_HOST_AND_DEVICE
    ResultType generateNth(const uint64_t& offset)
    {
      union
      {
	/// Temporary result of the encpter
	Word temp;
	
	/// Output to be returned
	const ResultType out{};
      };
      
      /// Shifted state
      State shiftedState=state+offset;
      
      temp=encrypt(key,static_cast<Word>(shiftedState));
      
      return out;
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
    Sitmo::Rng& ref;
    
    //Number of uint32_t generated so far
    uint64_t irng;
    
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
    RngView(Sitmo::Rng& ref,const uint64_t& irng) :
      ref(ref),irng(irng)
    {
    }
    
    /// Returns an uint32
    CUDA_HOST_AND_DEVICE
    uint32_t operator()()
    {
      return ref.generateNth(irng++);
    }
  };
  
  template <typename T>
  struct FieldRngOf
  {
    //Mark if the generator can be used or not
    bool used;
    
    //Reference to the Sitmo
    Sitmo::Rng& rng;
    
    //Distribution [0,1)
    std::uniform_real_distribution<> distr;
    
    //Number of reals to be shifted, in units of global volume
    const uint64_t offsetReal;
    
    //Number of reals for each site
    static constexpr int nRealsPerSite=
      sizeof(T)/sizeof(double);
    
    //Constructor
    FieldRngOf(Sitmo::Rng& rng,const uint64_t& offsetReal) :
      used(false),rng(rng),offsetReal(offsetReal)
    {
      //master_printf("All entries of FieldRng initialized\n");
    }
    
    /// Returns a view on a specific site and real number
    CUDA_HOST_AND_DEVICE
    RngView getRngViewOnGlbSiteIRndReal(const int& glblx,const int& irnd_real_per_site)
    {
      //Computes the number in the stream of reals
      const uint64_t irnd_double=offsetReal+glblx+nissa::glbVol*irnd_real_per_site;
      
      //Computes the offset in the rng stream
      const uint64_t irnd_uint32_t=2*irnd_double;
      
      return RngView(rng,irnd_uint32_t);
    }
    
    /// Fill a specific site
    CUDA_HOST_AND_DEVICE
    void fillGlbSite(T& out,const uint64_t glblx)
    {
      for(int irnd_real=0;irnd_real<nRealsPerSite;irnd_real++)
	{
	  auto view=getRngViewOnGlbSiteIRndReal(glblx,irnd_real);
	  
	  ((double*)out)[irnd_real]=distr(view);
	}
    }
    
    /// Fill a specific site given its local index
    CUDA_HOST_AND_DEVICE
    void fillLocSite(T& out,const uint64_t loclx)
    {
      //Finds the global site of local one
      const int& glblx=glblxOfLoclx[loclx];
      fillGlbSite(out,glblx);
    }
    
    void enforce_single_usage()
    {
      if(used)
	crash("cannot use a rng field twice");
      used=true;
    }
    
    /// Fill with origin
    void fillWithOrigin(T& out)
    {
      enforce_single_usage();
      
      fillGlbSite(out,0);
    }
    
    /// Fill all sites
    void fillField(T* out)
    {
      enforce_single_usage();
      
      NISSA_PARALLEL_LOOP(loclx,0,locVol)
	{
	  const int& glblx=glblxOfLoclx[loclx];
	  fillGlbSite(out[glblx],glblx);
	}
      NISSA_PARALLEL_LOOP_END;
      
      set_borders_invalid(out);
    }
  };
  
  struct FieldRngStream
  {
    //Embeds the random number generator
    Sitmo::Rng rng;
    
    //Number of double precision numbers generated per site
    uint64_t ndouble_gen;
    
    //Skip n drawers
    template <typename T>
    void skipDrawers(const uint64_t& nDrawers)
    {
      ndouble_gen+=FieldRngOf<T>::nRealsPerSite*nDrawers;
    }
    
    //Return a proxy from which to fetch the random numbers
    template <typename T>
    auto getDrawer()
    {
      static_assert(sizeof(T)%sizeof(double)==0,"Type not multiple in size of double");
      
      //Result
      FieldRngOf<T> res(rng,ndouble_gen);
      
      skipDrawers<T>(1);
      
      return res;
    }
    
    //Draw a single instance of T
    template <typename T>
    void drawScalar(T& out)
    {
      static_assert(sizeof(T)%sizeof(double)==0,"Type not multiple in size of double");
      
      auto drawer=getDrawer<T>();
      
      drawer.fillWithOrigin(out);
    }
    
    //Initializes with a seed
    void init(const uint32_t& seed)
    {
      rng.seed(seed);
      ndouble_gen=0;
    }
  };
  
  static constexpr complex zero_complex={0.0,0.0};
  static constexpr complex sqrt_2_half_complex{M_SQRT1_2,M_SQRT1_2};
  
  inline void BoxMullerTransform(complex out,const complex ave=zero_complex,const complex sig=sqrt_2_half_complex)
  {
    const double r=sqrt(-2*log(1-out[RE]));
    const double q=2*M_PI*out[IM];
    
    out[RE]=r*cos(q)*sig[RE]+ave[RE];
    out[IM]=r*sin(q)*sig[IM]+ave[IM];
  }
  
  CUDA_HOST_AND_DEVICE
  inline void z2Transform(double& out)
  {
    out=(out>0.5)?M_SQRT1_2:-M_SQRT1_2;
  }
  
  CUDA_HOST_AND_DEVICE
  inline void z4Transform(complex out)
  {
    for(int ri=0;ri<2;ri++)
      z2Transform(out[ri]);
  }
}

#endif
