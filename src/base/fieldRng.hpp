#ifndef _FIELD_RNG_HPP
#define _FIELD_RNG_HPP

#include <base/field.hpp>
#include <base/randomDistributions.hpp>
#include <base/sitmo.hpp>
#include <geometry/geometry_lx.hpp>
#include <threads/threads.hpp>

// Lattice drawers assign double index (offsetReal+component)*volume+globalSite
// The stream reserves components per site when a drawer is requested; actual
// values are generated only when that drawer is accessed

namespace nissa
{
  /// Lazy site-indexed uniform drawer backed by the selected Sitmo engine
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
	
	VERBOSITY_LV3_MASTER_PRINTF("Added backup callback for FieldRngStream\n");
	
	return true;
      }
    else
      {
	VERBOSITY_LV3_MASTER_PRINTF("Skipping backup callbacks for FieldRngStream as we are not inside a parallel for\n");
	
	return false;
      }
#else
      return false;
#endif
    }
  };
}

#endif
