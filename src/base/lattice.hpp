#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <memory>
#include <sstream>

#include <expr/cWiseCombine.hpp>
#include <expr/compReduce.hpp>
#include <expr/fieldDeclaration.hpp>
#include <expr/prod.hpp>
#include <operations/allToAll.hpp>
#include <routines/mpiRoutines.hpp>

namespace nissa
{
  /// Global coordinates
  using GlbCoords=
    Coords<GlbCoord>;
  
  /// Local coordinates
  using LocCoords=
    Coords<LocCoord>;
  
  /// Mpi Rank coordinates
  using MpiRankCoords=
    Coords<MpiRankCoord>;
  
  /////////////////////////////////////////////////////////////////
  
  /// Holds all the info on the lattice
  template <bool IsRef=false>
  struct Lattice
  {
    /// Index of time direction
    static constexpr Dir timeDir=0;
    
    /// Mapping of ILDG directons
    static constexpr Coords<Dir> scidacDirOfNissaDir=
      [](const Dir& in)
      {
	if(in==timeDir)
	  return in;
	else
	  return nDim-in;
      };
    
    /// Backward substitution
    static constexpr Coords<Dir> nissaDirOfScidacDir=
      scidacDirOfNissaDir;
    
#define PROVIDE_MEMBER_WITH_ACCESSOR(TYPE,NAME,CAP_NAME)	\
    TYPE _ ## NAME;						\
    								\
    template <typename...Args>					\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION		\
    decltype(auto) get ## CAP_NAME(Args&&...args) const		\
    {								\
      if constexpr(sizeof...(args))				\
	return _ ## NAME(std::forward<Args>(args)...);		\
      else							\
	return _ ## NAME;					\
    }
    
    PROVIDE_MEMBER_WITH_ACCESSOR(GlbLxSite,glbVol,GlbVol);
    
    PROVIDE_MEMBER_WITH_ACCESSOR(LocLxSite,locVol,LocVol);
    
    PROVIDE_MEMBER_WITH_ACCESSOR(LocLxSite,surfSize,SurfSize);
    
    PROVIDE_MEMBER_WITH_ACCESSOR(LocCoords,surfSizePerDir,SurfSizePerDir);
    
    PROVIDE_MEMBER_WITH_ACCESSOR(LocCoords,surfOffsetOfDir,SurfOffsetOfDir);
    
    PROVIDE_MEMBER_WITH_ACCESSOR(GlbCoords,glbSizes,GlbSizes);
    
    PROVIDE_MEMBER_WITH_ACCESSOR(LocCoords,locSizes,LocSizes);
    
    using LocLxToLocLxMap=
      MirroredTens<OfComps<LocLxSite>,LocLxSite,IsRef>;
    
    PROVIDE_MEMBER_WITH_ACCESSOR(LocLxToLocLxMap,surfSiteOfHaloSite,SurfSiteOfHaloSite);
    
    using LocLxToLocNeighsMap=
      MirroredTens<OfComps<Ori,LocLxSite,Dir>,LocLxSite,IsRef>;
    
    PROVIDE_MEMBER_WITH_ACCESSOR(LocLxToLocNeighsMap,locLxNeigh,LocLxNeigh);
    
    using LocLxToGlbLxMap=
      MirroredTens<OfComps<LocLxSite>,GlbLxSite,IsRef>;
    
    PROVIDE_MEMBER_WITH_ACCESSOR(LocLxToGlbLxMap,glbLxOfLocLx,GlbLxOfLocLx);
    
    using LocLxToGlbCoordsMap=
      MirroredTens<OfComps<LocLxSite,Dir>,GlbCoord,IsRef>;
    
    PROVIDE_MEMBER_WITH_ACCESSOR(LocLxToGlbCoordsMap,glbCoordsOfLocLx,GlbCoordsOfLocLx);
    
#undef PROVIDE_MEMBER_WITH_ACCESSOR
    
    /////////////////////////////////////////////////////////////////
    
    /// Compute internal volume
    template <DerivedFromComp C>
    constexpr INLINE_FUNCTION
    static C bulkVolume(const Coords<C>& L)
    {
      C out=1;
      
      for(Dir mu=0;mu<nDim;mu++)
	{
	  if(L(mu)>2) out*=L(mu)-2;
	  else out=0;
	}
      
      return out;
    }
    
    /// Compute the variance of the surface
    template <DerivedFromComp C>
    constexpr INLINE_FUNCTION
    static double computeSurfaceVariance(const Coords<C>& L)
    {
      const Coords<C> s=
	compProd<Dir>(L)/L;
      
      return
	compSum<Dir>(sqr(s))()/(double)nDim-
	sqr(compSum<Dir>(s)()/(double)nDim);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Returns the optimal paritioning {r}, which divides {l}={g}/{r}
    template <DerivedFromComp C,
	      DerivedFromComp G>
    static Coords<C> findOptimalPartitioning(const int64_t& r,
					     const Coords<G>& g,
					     const Coords<C>& fixR)
    {
      constexpr bool partitionDebug=false;
      
      /// Result
      Coords<C> res;
      res=0;
      
      [[maybe_unused]]
      auto coordsToStr=
	[]<typename T>(const Coords<T>& c)
	{
	  std::ostringstream os;
	  os<<"Coords["<<demangle(typeid(T).name()).c_str()<<"] : {"<<c[timeDir]();
	  for(Dir mu=1;mu<nDim;mu++)
	    os<<","<<c[mu]();
	  os<<"}";
	  
	  return os.str();
	};
      
      [[maybe_unused]]
      auto vectToStr=
	[](const std::vector<int64_t>& c)
	{
	  std::ostringstream os;
	  os<<"{"<<c[timeDir()];
	  for(Dir mu=1;mu<nDim;mu++)
	    os<<","<<c[mu()];
	  os<<"}";
	  
	  return os.str();
	};
      
      const int64_t l=
	compProd<Dir>(g)()()/r;
      
      if constexpr(partitionDebug)
	master_printf("Going to factorize l=%ld, r=%ld, g=%s with constraints: %s\n",l,r,coordsToStr(g).c_str(),coordsToStr(fixR).c_str());
      
      /// Factors of l
      const std::vector<int64_t> lFactors=
	factorize(l);
      
      if constexpr(partitionDebug)
	master_printf("l factors: %s\n",vectToStr(lFactors).c_str());
      
      /// Factors of r
      const std::vector<int64_t> rFactors=
	factorize(r);
      
      if constexpr(partitionDebug)
	master_printf("r factors: %s\n",vectToStr(rFactors).c_str());
      
      /// Chooses if to factorize r, otherwise l
      const bool factorizeR=
	lFactors.size()>=rFactors.size();
      
      master_printf("factorize r: %d\n",factorizeR);
      
      /// Copy the facts to be used
      const std::vector<int64_t>& facts=
	factorizeR?
	rFactors:
	lFactors;
      
       const int64_t nFacts=
	 facts.size();
      
      if constexpr(partitionDebug)
	master_printf("nFacts: %ld\n",nFacts);
      
      /// We will decompose the index putting each factor to the correct direction
      const int64_t baseMask=
	pow(nDim(),nFacts-1);
      
      if constexpr(partitionDebug)
	master_printf("baseMask: %ld\n",baseMask);
      
      /// Total number of combo
      const int64_t nCombo=
	baseMask*nDim();
      
      master_printf("nCombo: %ld\n",nCombo);
      
      int64_t minLocSurf=l;
      double minSurfVariance=-1;
      int64_t nextSupposedValid=0;
      for(int64_t iCombo=0;iCombo<nCombo;)
	{
	  if constexpr(partitionDebug)
	    master_printf("iCombo: %ld/%ld\n",iCombo,nCombo);
	  
	  //number of ranks in each direction for current partitioning
	  Coords<C> R;
	  R=1;
	  
	  //find the partioning corresponding to icombo
	  int64_t mask=baseMask;
	  int64_t iFact=facts.size()-1;
	  bool validPartitioning=true;
	  do
	    {
	      if constexpr(partitionDebug)
		{
		  master_printf(" iFact: %ld/%ld\n",iFact,nFacts);
		  master_printf(" mask %ld\n",mask);
		}
	      
	      /// Direction: this is given by the ifact digit of icombo wrote in base nDim
	      const Dir mu=
		(int)((iCombo/mask)%nDim());
	      
	      if constexpr(partitionDebug)
		master_printf(" mu: %d\n",mu());
	      
	      const int64_t f=
		facts[iFact];
	      R[mu]*=f;
	      
	      if constexpr(partitionDebug)
		master_printf(" fact: %ld\n",f);
	      
	      //check that the total volume L is a multiple and it is larger than the number of proc
	      validPartitioning&=
		(g[mu]()%R[mu]()==0 and g[mu]()>=R[mu]());
	      
	      if constexpr(partitionDebug)
		master_printf("g[mu] mod R[mu] = %ld mod  %ld = %ld\n",g[mu](),R[mu](),g[mu]()%R[mu]());
	      
	      if(validPartitioning)
		{
		  iFact--;
		  mask/=nDim();
		}
	      
	      if constexpr(partitionDebug)
		{
		  master_printf("valid: %d \n",validPartitioning);
		  printf("\n");
		}
	    }
	  while(validPartitioning and iFact>=0);
	  
	  if(not factorizeR)
	    for(Dir mu=0;mu<nDim;mu++)
	      R[mu]=g[mu]()/R[mu];
	  
	  for(Dir mu=0;mu<nDim;mu++)
	    if(fixR[mu]())
	      validPartitioning&=(fixR[mu]==R[mu]);
	  
	  if constexpr(partitionDebug)
	    master_printf("Valid: %d (should not before: %ld)\n",validPartitioning,nextSupposedValid);
	  
	  //validity coulde have changed
	  if(validPartitioning)
	    {
	      const Coords<C> tryLSizes=
		g.template reinterpretFund<C>()/R;
	      
	      const int64_t locSurf=
		l-bulkVolume(tryLSizes)();
	      
	      const double surfVariance=
		computeSurfaceVariance(tryLSizes);
		
	      if(locSurf<minLocSurf or
		 (locSurf==minLocSurf and
		  surfVariance<minSurfVariance))
		{
		  minLocSurf=locSurf;
		  minSurfVariance=surfVariance;
		  
		  res=R;
		}
	      
	      if constexpr(partitionDebug)
		if(iCombo<nextSupposedValid)
		  master_printf("hey! it was not supposed to be valid!\n");
	      
	      iCombo++;
	    }
	  //skip all remaining factorization using the same structure
	  else
	    {
	      const int64_t skip=
		pow(nDim(),std::max(0l,iFact-1));
	      
	      // iCombo+=skip;
	      if(iCombo>=nextSupposedValid)
		{
		  nextSupposedValid=iCombo+skip;
		  if constexpr(partitionDebug)
		    master_printf("Would need to skip up to: %ld\n",iCombo+skip);
		}
	      if constexpr(partitionDebug)
		iCombo++;
	      else
		iCombo=nextSupposedValid;
	    }
	  if constexpr(partitionDebug)
	    printf("\n");
	}
      
      return res;
    }
    
    /// Sets the global sizes
    void setGlbSizes(const GlbCoords& extGlbSizes)
    {
      _glbSizes=extGlbSizes;
      _glbVol=compProd<Dir>(getGlbSizes()).close()();
      
      master_printf("Global lattice:\t%ld",getGlbSizes().dirRow(0));
      for(Dir mu=1;mu<NDIM;mu++)
	master_printf("x%ld",getGlbSizes()[mu]());
      master_printf(" = %ld\n",getGlbVol());
    }
    
    /// Sets the mpi ranks
    void setMpiRanks()
    {
      using namespace resources;
      
      if(const MpiRankCoords temp=
	 Lattice<>::findOptimalPartitioning<MpiRankCoord>(nRanks(),getGlbSizes(),MpiRankCoords{});compProd<Dir>(temp).close()==0)
	crash("Unable to partition");
      else
	_nRanksPerDir=temp;
      
      for(Dir dir=0;dir<nDim;dir++)
	_isDirParallel[dir]=(nRanksPerDir[dir]()>1);
      
      _thisRankCoords=decomposeLxToCoords(thisRank(),nRanksPerDir);
      
      for(Ori ori=0;ori<2;ori++)
	for(Dir dir=0;dir<nDim;dir++)
	  {
	    MpiRankCoords neighRankCoords=thisRankCoords;
	    neighRankCoords[dir]=(neighRankCoords[dir]+nRanksPerDir[dir]+(2*ori()-1))%nRanksPerDir[dir];
	    _neighRanks(ori,dir)=lxOfCoords<MpiRank>(neighRankCoords,nRanksPerDir);
	  }
    }
    
    /// Sets the local sizes
    void setLocSizes()
    {
      _locVol=getGlbVol()()/nRanks();
      
      for(Dir dir=0;dir<nDim;dir++)
	_locSizes()[dir]=
	  getGlbSizes()[dir]()/nRanksPerDir[dir]();
    }
    
    /// Sets the surface sizes
    void setSurfSizes()
    {
      _surfSizePerDir=
	getLocVol()/getLocSizes();
      
      _surfSize=compSum<Dir>(getSurfSizePerDir());
      
      _surfOffsetOfDir[timeDir]=0;
      for(Dir dir=1;dir<nDim;dir++)
	_surfOffsetOfDir[dir]=_surfOffsetOfDir[dir-1]+getSurfSizePerDir()[dir-1];
    }
    
    /// Initializes
    void init(const GlbCoords& extGlbSizes)
      requires(not IsRef)
    {
      setGlbSizes(extGlbSizes);
      
      setMpiRanks();
      
      setLocSizes();

      //glbCoordsOfLocLx;
    }
    
    /// Initializes from T and L
    void init(const GlbCoord& T,
	      const GlbCoord& L)
      requires(not IsRef)
    {
      GlbCoords _glbSizes=L;
      _glbSizes[timeDir]=T;
      
      this->init(_glbSizes);
    }
    
    /// Gets a global coordinate
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    auto glbCoord(const Dir& dir) const
    {
      return getGlbCoordsOfLocLx(dir);
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
	  isSpatOrigin&=(lat.getGlbCoordsOfLocLx(site)(nu)==0);
	
	return isSpatOrigin;
      }
    };
    
    /// Returns a function which evaluates to true on spatial origins
    auto spatialOriginsMask() const
    {
      return
	funcNodeWrapper<CompsList<LocLxSite>>(SpatOriginMaskFunctor{this->getRef()},std::make_tuple(getLocVol()));
    }
    
    /// Default constructor
    Lattice() = default;
    
#define COPY_CONSTRUCTOR_BODY					\
    _glbVol(oth._glbVol),					\
      _locVol(oth._locVol),					\
      _glbSizes(oth._glbSizes),					\
      _locSizes(oth._locSizes)/*,					\
				_glbCoordsOfLocLx(oth._glbCoordsOfLocLx.getRef())*/ \
    {								\
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
  inline Lattice<>* _lat;
  
  /// Reference to the lattice
  inline std::unique_ptr<Lattice<true>> lat;
}

#endif
