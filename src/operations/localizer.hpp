#ifndef _LOCALIZER_HPP
#define _LOCALIZER_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <expr/comp.hpp>
#include <expr/comps.hpp>
#include <expr/conj.hpp>
#include <expr/field.hpp>
#include <expr/mergedComps.hpp>
#include <operations/allToAll.hpp>

#ifndef EXTERN_LOCALIZER
# define EXTERN_LOCALIZER extern
# define INITIALIZE_LOCALIZER_TO(ARGS...)
#else
# define INITIALIZE_LOCALIZER_TO(ARGS...) ARGS
#endif

namespace nissa::localizer
{
  DECLARE_DYNAMIC_COMP(OrthoSpaceTime);
  
  DECLARE_DYNAMIC_COMP(FullLocDirCoord);
  
  /// Merged component for the temporary storage
  using MC=
    MergedComp<CompsList<OrthoSpaceTime,FullLocDirCoord>>;
  
  /// Type of the communicator which makes local a given direction
  using LocDirMaker=
    AllToAllComm<MC,LocLxSite>;
  
  /// Type for the communicator which changes the current local direction to another one
  using LocDirChanger=
    AllToAllComm<MC,MC>;
  
  /// Type for the communicator which makes ordinary a local direction
  using LocDirUnmaker=
    AllToAllComm<LocLxSite,MC>;
  
  /// Communicator which makes local the first direction
  EXTERN_LOCALIZER LocDirMaker* firstLocDirMaker;
  
  /// Communicators which transform the local direction to the next in list
  EXTERN_LOCALIZER std::vector<LocDirChanger> locDirChanger;
  
  /// Communicator which brings back the local direction to lexicograèhic
  EXTERN_LOCALIZER LocDirUnmaker* lastLocDirUnmaker;
  
  /// Store whether the cycler has been initialized
  EXTERN_LOCALIZER bool initialized INITIALIZE_LOCALIZER_TO({false});
  
  /// Dimensions to be used for the temporary storage
  EXTERN_LOCALIZER StackTens<CompsList<Dir>,std::tuple<FullLocDirCoord,OrthoSpaceTime>> dimensions;
  
  /// Computes the dimension for each direction
  std::tuple<FullLocDirCoord,OrthoSpaceTime> computeDimensions(const Dir& dir);
  
  /// Creates the communicator which makes local a given direction
  LocDirMaker getLocDirMaker(const Dir& dir);
  
  /// Initializes the communicators
  void init();
  
  /// Release the communicators
  void dealloc();
  
  /// Internal implementation making a dimension local in turn
  template <DerivedFromNode Out,
	    DerivedFromNode In,
	    typename CompsBef,
	    typename CompsAft>
  struct _Transform
  {
    /// Buffer components
    using BufComps=
      TupleCat<CompsBef,CompsList<OrthoSpaceTime,FullLocDirCoord>,CompsAft>;
    
    /// Fundamental type of the expression to cycle
    using Fund=
      typename In::Fund;
    
    /// Buffer to be used
    using Buf=
      DynamicTens<BufComps,Fund,getMemoryType<In::execSpace>()>;
    
    /// Types needed to store extra dynamic sizes
    using ExtraDynamicComps=
      In::DynamicComps;
    
    /// Extra dynamic sizes
    const ExtraDynamicComps extraDynamicSizes;
    
    /// Construct taking extra dynamic sizes
    _Transform(const ExtraDynamicComps& extraDynamicSizes) :
      extraDynamicSizes(extraDynamicSizes)
    {
    }
    
    /// Prepare for the first iteration
    template <typename F>
    auto prepare(const In& in,
		 const F& f) const
    {
      /// Result
      auto res=
	getBuf(nDim-1);
      
      firstLocDirMaker->communicate(res,in);
      
      f(res,std::integral_constant<Dir,nDim-1>());
      
      return res;
    }
    
    template <DerivedFromComp...Dc>
    auto getBuf(const Dir& D) const
    {
      /// Dynamical components of the result
      const auto dynamicSizes=
	std::tuple_cat(extraDynamicSizes,dimensions[D]);
      
      return
	Buf(dynamicSizes).template mergeComps<CompsList<OrthoSpaceTime,FullLocDirCoord>>();
    }
    
    /// Iteration
    template <Dir D,
	      DerivedFromNode N,
	      typename F>
    requires(D>=0 and D<nDim)
    decltype(auto) iter(N&& n,
			const F& f) const
    {
      /// This direction in a constant form
      constexpr auto DIC=
	std::integral_constant<Dir,D>();
      
      if constexpr(D==nDim-1)
	return prepare(std::forward<N>(n),f);
      else
	{
	  /// Inner part to be used
	  const auto inner=
	    iter<D+1>(n,f);
	  
	  /// Returned result
	  auto res=
	    getBuf(D);
	  
	  locDirChanger[D()].communicate(res,inner);
	  
	  f(res,DIC);
	  
	  return res;
	}
    }
    
    /// Cycle on all directions localization
    template <typename F>
    void cycle(Out& out,
	       const In& in,
	       const F& f)
    {
      lastLocDirUnmaker->communicate(out,iter<Dir(0)>(in,f));
    }
  };
  
  /// Computes the dimensions for the localizer
  inline std::tuple<FullLocDirCoord,OrthoSpaceTime> computeDimensions(const Dir& dir)
  {
    /// Fully local size in a given direction is the global size
    const FullLocDirCoord fcSize=
      lat->getGlbSizes(dir)();
    
    /// Perpendicular size across the whole lattice
    const OrthoSpaceTime glbOsdSize=
      lat->getGlbVol()()/lat->getGlbSizes(dir)();
    
    /// Portion of the perpendicular size relative to each lattice
    const OrthoSpaceTime locOsdSize=
      (glbOsdSize+nRanks()-1)/nRanks();
    
    return {fcSize,locOsdSize};
  }
  
  /// Returns a localizer for a given direction
  inline LocDirMaker getLocDirMaker(const Dir& dir)
  {
    return {lat->getLocVol(),
	    [dir](const LocLxSite& locLxSite)
	    {
	      /// Dimensions of the current direction
	      const auto& [fcSize,locOsdSize]=
		dimensions[dir];
	      
	      /// Coordinate in the current direction of the requires site
	      const FullLocDirCoord fc=
		lat->getGlbCoordsOfLocLx(locLxSite)[dir]();
	      
	      /// Global index in the space perpendicular to the current direction
	      int64_t glbOsd=0;
	      for(PerpDir pDir=0;pDir<NDIM-1;pDir++)
		{
		  const Dir jDir=perpDirOf[dir][pDir];
		  glbOsd=lat->getGlbSizes(jDir)()*glbOsd+lat->getGlbCoordsOfLocLx(locLxSite)(jDir)();
		}
	      
	      /// Rank hosting global site
	      const MpiRank orthoRank=
		(int)(glbOsd/locOsdSize());
	      
	      /// Local index in the rank
	      const OrthoSpaceTime locOsd=
		glbOsd-orthoRank()*locOsdSize();
	      
	      /// Merges the components
	      const auto mc=
		MC::merge(std::make_tuple(fcSize,locOsdSize),fc,OrthoSpaceTime(locOsd));
	      
	      return std::make_tuple(orthoRank,mc);
	    }};
  }
  
  /// Initialize the localizers
  inline void init()
  {
    if(initialized)
      crash("already initialized");
    
    for(Dir dir=0;dir<NDIM;dir++)
      dimensions[dir]=computeDimensions(dir);
    
    /// Communicators which make a given direction local
    std::vector<LocDirMaker> locDirMaker;
    for(Dir dir=0;dir<NDIM;dir++)
      locDirMaker.push_back(getLocDirMaker(dir));
    
    /// Inverse localizer
    std::vector<LocDirUnmaker> locDirUnmaker;
    for(Dir dir=0;dir<NDIM;dir++)
      locDirUnmaker.push_back(locDirMaker[dir()].inverse());
    
    for(Dir dir=0;dir<nDim-1;dir++)
      locDirChanger.push_back(locDirMaker[dir()]*locDirUnmaker[dir()+1]);
    
    firstLocDirMaker=new LocDirMaker(std::move(locDirMaker.back()));
    
    lastLocDirUnmaker=new LocDirUnmaker(std::move(locDirUnmaker.front()));
    
    initialized=true;
  }
  
  /// Deallocate the localizers
  inline void dealloc()
  {
    if(initialized)
      {
	delete firstLocDirMaker;
	locDirChanger.clear();
	delete lastLocDirUnmaker;
      }
    else
      crash("localizer not initialized");
    
    initialized=false;
  }
}

namespace nissa
{
  /// Cycle making each direction local in turn, executing the function f
  template <typename BufCompsBef,
	    typename BufCompsAft,
	    DerivedFromNode Out,
	    typename F,
	    DerivedFromNode In,
	    DerivedFromComp...D>
  void cycleOnAllLocalDirections(Out&& out,
				 const In& in,
				 const F& f,
				 const CompsList<D...>& dynamicSizes)
  {
    localizer::_Transform<Out,In,BufCompsBef,BufCompsAft> cycler(dynamicSizes);
    cycler.cycle(out,in,f);
  }
}

#undef EXTERN_LOCALIZER
#undef INITIALIZE_LOCALIZER_TO

#endif
