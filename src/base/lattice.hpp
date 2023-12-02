#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <memory>

#include <expr/cWiseCombine.hpp>
#include <expr/field.hpp>
#include <operations/allToAll.hpp>

#ifndef EXTERN_LATTICE
# define EXTERN_LATTICE extern
#endif

namespace nissa
{
  /// Holds all the info on the lattice
  template <bool IsRef=false>
  struct Lattice
  {
    /// Global volume
    GlbLxSite glbVol;
    
    /// Local volume
    LocLxSite locVol;
    
    /// Gloabl extension in each direction
    GlbCoords glbSizes;
    
    /// Local extension in each direction
    LocCoords locSizes;
    
    /// Global sites of local sites
    MirroredTens<OfComps<LocLxSite>,ConstIf<IsRef,GlbLxSite>,IsRef> glbLxOfLocLx;
    
    /// Global coordinates of local sites
    MirroredTens<OfComps<LocLxSite,Dir>,ConstIf<IsRef,GlbCoord>,IsRef> glbCoordsOfLocLx;
    
    /// Initializes
    void init()
      requires(not IsRef)
    {
      locVol=nissa::locVol;
      glbVol=nissa::glbVol;
      locSizes=[](const Dir& dir){return nissa::locSize[dir()];};
      glbSizes=[](const Dir& dir){return nissa::glbSizes[dir()];};
      glbCoordsOfLocLx.allocate(std::make_tuple(locVol));
      glbCoordsOfLocLx.getFillable()=
	[g=nissa::glbCoordOfLoclx](const LocLxSite& site,
				   const Dir& dir)
	{
	  return g[site()][dir()];
	};
      
      glbLxOfLocLx.allocate(std::make_tuple(locVol));
      glbLxOfLocLx.getFillable()=
	[g=nissa::glblxOfLoclx](const LocLxSite& site)
	{
	  return g[site()];
	};
    }
    
    /// Gets a global coordinate
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    auto glbCoord(const Dir& dir) const
    {
      return glbCoordsOfLocLx(dir);
    }
    
    /// Functor used to mak the spatial origin
    struct SpatOriginMaskFunctor
    {
      /// Can run on both GPU and CPU as it is trivially copyable
      static constexpr ExecSpace execSpace=
	execOnCPUAndGPU;
      
      Lattice lat;
      
      constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
      SpatOriginMaskFunctor(const Lattice& lat)
	requires(IsRef) :
	lat(lat)
      {
      }
      
      SpatOriginMaskFunctor(const SpatOriginMaskFunctor&) = default;
      
      /// Evaluate
      constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
      auto operator()(const LocLxSite& site) const
      {
	bool isSpatOrigin=true;
	
	for(Dir nu=1;nu<nDim;nu++)
	  isSpatOrigin&=(lat.glbCoordsOfLocLx(site,nu)==0);
	
	return isSpatOrigin;
      }
    };
    
    /// Returns a function which evaluates to true on spatial origins
    auto spatialOriginsMask() const
    {
      return
	funcNodeWrapper<CompsList<LocLxSite>>(SpatOriginMaskFunctor{this->getRef()},std::make_tuple(locVol));
    }
    
    /// Default constructor
    Lattice() = default;
    
#define COPY_CONSTRUCTOR_BODY				\
    glbVol(oth.glbVol),					\
      locVol(oth.locVol),				\
      glbSizes(oth.glbSizes),				\
      locSizes(oth.locSizes),				\
      glbCoordsOfLocLx(oth.glbCoordsOfLocLx.getRef())	\
    {							\
    }
    
    /// Copy construct from reference
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    Lattice(const Lattice<false>& oth)
      requires(IsRef) :
      COPY_CONSTRUCTOR_BODY;
    
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    Lattice(const Lattice& oth) :
      COPY_CONSTRUCTOR_BODY;
    
#undef COPY_CONSTRUCTOR_BODY
    
    /// Returns a reference
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    Lattice<true> getRef() const
    {
      return *this;
    }
    
    /// Returns a versor in the direction extDir
    template <typename Fund=bool>
    constexpr static auto getVersor(const Dir& extDir)
    {
      StackTens<CompsList<Dir>,Fund> res;
      for(Dir dir=0;dir<nDim;dir++)
	res[dir]=dir==extDir;
      
      return res;
    }
    
    /// Hypercube diagonal
    static constexpr StackTens<CompsList<Dir>> hCubeDiag=1.0;
    
    /// Versor in a given direction
    static constexpr StackTens<CompsList<Dir>,StackTens<CompsList<Dir>,bool>> versors=
      getVersor<bool>;
    
    /// List of perpendicular directions
    static constexpr StackTens<CompsList<Dir>,StackTens<CompsList<Dir>,bool>> perpDirs=
      hCubeDiag-versors;
  };
  
  /// Stores the actual lattice
  EXTERN_LATTICE Lattice<>* _lat;
  
  /// Reference to the lattcice
  EXTERN_LATTICE std::unique_ptr<Lattice<true>> lat;
}

#undef EXTERN_LATTICE

#endif
