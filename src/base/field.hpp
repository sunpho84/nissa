#ifndef _FIELD_HPP
#define _FIELD_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#ifdef USE_CUDA
# include <thrust/complex.h>
# include <thrust/execution_policy.h>
# include <thrust/inner_product.h>
# include <thrust/reduce.h>
#endif

#include <cstddef>
#include <mpi.h>
#include <optional>
#include <type_traits>

#include <base/bench.hpp>
#include <base/memory_manager.hpp>
#include <base/vectors.hpp>
#include <communicate/communicate.hpp>
#include <geometry/geometry_eo.hpp>
#include <geometry/geometry_lx.hpp>
#include <linalgs/reduce.hpp>
#include <metaprogramming/constnessChanger.hpp>
#include <metaprogramming/extent.hpp>
#include <metaprogramming/feature.hpp>
#include <metaprogramming/unroll.hpp>
#include <routines/ios.hpp>
#include <threads/benchmarks.hpp>

namespace nissa
{
  /// Start the communications of buffer interpreted as halo
  std::vector<MPI_Request> startBufHaloNeighExchange(const int& divCoeff,
						     const size_t& bps,
						     const std::optional<std::vector<std::pair<int,int>>>& dirs=std::nullopt);
  
  /// Start the communications of buffer interpreted as halo
  template <typename T>
  std::vector<MPI_Request> startBufHaloNeighExchange(const int& divCoeff,
						     const std::optional<std::vector<std::pair<int,int>>>& dirs=std::nullopt)
  {
    return startBufHaloNeighExchange(divCoeff,sizeof(T),dirs);
  }
  
  /// Start the communications of buffer interpreted as edges
  std::vector<MPI_Request> startBufEdgesNeighExchange(const int& divCoeff,
						      const size_t& bps);
  
  /// Start the communications of buffer interpreted as edges
  template <typename T>
  std::vector<MPI_Request> startBufEdgesNeighExchange(const int& divCoeff)
  {
    return startBufEdgesNeighExchange(divCoeff,sizeof(T));
  }
  
  /// Wait for communications to finish
  inline void waitAsyncCommsFinish(std::vector<MPI_Request> requests)
  {
    VERBOSITY_LV3_MASTER_PRINTF("Entering MPI comm wait\n");
    
    MPI_Waitall(requests.size(),&requests[0],MPI_STATUS_IGNORE);
  }
  
  /// Communicates the buffers as halo
  template <typename T>
  void exchangeHaloNeighBuf(const int& divCoeff)
  {
    waitAsyncCommsFinish(startBufHaloNeighExchange<T>(divCoeff));
  }
  
  /// Communicates the buffers as edges
  template <typename T>
  void exchangeEdgesNeighBuf(const int& divCoeff)
  {
    waitAsyncCommsFinish(startBufEdgesNeighExchange<T>(divCoeff));
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Subscription of the field
  ///
  /// Forward declaration
  template <typename F,
	    typename P>
  struct SubscribedField;
  
  /////////////////////////////////////////////////////////////////
  
  /// Arrangement of internal DOF w.r.t spacetime: AoS or SoA
  enum class SpaceTimeLayout{CPU,GPU};
  
  /// Crashes if not running on the given memory space
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  void assertRunningOnMemorySpace(const MemorySpace& ms)
  {
    if(ms!=currentMemorySpace)
      CRASH("should be running on memory space %s, is running on %s",
	    memorySpaceName(ms),memorySpaceName(currentMemorySpace));
  }
  
  /// Coverage of the field
  enum FieldCoverage{EVEN_SITES,ODD_SITES,FULL_SPACE,EVEN_OR_ODD_SITES};
  
  /// Has or not the halo and the edges
  enum HaloEdgesPresence{WITHOUT_HALO,WITH_HALO,WITH_HALO_EDGES};
  
  /// Predefinite arrangement of internal DOF
  constexpr SpaceTimeLayout defaultSpaceTimeLayout=
	      SpaceTimeLayout::DEFAULT_MEMORY_LAYOUT;
  
  /// Returns the name of the spacetime layout sl
  INLINE_FUNCTION
  constexpr const char* spaceTimeLayoutName(const SpaceTimeLayout& sl)
  {
    switch(sl)
      {
      case SpaceTimeLayout::CPU:
	return "CPU";
	break;
      case SpaceTimeLayout::GPU:
	return "GPU";
	break;
      }
    
    return "";
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Number of sites contained in the field
  template <FieldCoverage fieldCoverage>
  struct FieldSizes
  {
    /// Assert that the coverage is definite
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static void assertHasDefinedCoverage()
    {
      static_assert(fieldCoverage==FULL_SPACE or
		    fieldCoverage==EVEN_SITES or
		    fieldCoverage==ODD_SITES,
		    "Trying to probe some feature of a field with unknwon coverage. If you are accessing an EvnOrOddField subfield, do it without subscribing them, or cast to the specific subspace");
    }
    
    /// Number of sites covered by the field
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr const int64_t& nSites()
    {
      if constexpr(fieldCoverage==FULL_SPACE)
	return locVol;
      else
	return locVolh;
    }
    
    /// Number of sites in the halo of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr const int64_t& nHaloSites()
    {
      if constexpr(fieldCoverage==FULL_SPACE)
	return bordVol;
      else
	return bordVolh;
    }
    
    /// Number of sites in the edges of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr const int64_t& nEdgesSites()
    {
      if constexpr(fieldCoverage==FULL_SPACE)
	return edgeVol;
      else
	return edgeVolh;
    }
    
    /// Computes the size to allocate
    static int64_t nSitesToAllocate(const HaloEdgesPresence& haloEdgesPresence)
    {
      int64_t res=nSites();
      
      if(haloEdgesPresence>=WITH_HALO)
	res+=nHaloSites();
      
      if(haloEdgesPresence>=WITH_HALO_EDGES)
	res+=nEdgesSites();
      
      return res;
    }
    
    /// Surface site of a site in the halo
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static auto surfSiteOfHaloSite(const int64_t& iHalo)
    {
      assertHasDefinedCoverage();
      
      if constexpr(fieldCoverage==FULL_SPACE)
	return surflxOfBordlx[iHalo];
      else
	if constexpr(fieldCoverage==EVEN_SITES or fieldCoverage==ODD_SITES)
	  return surfeo_of_bordeo[fieldCoverage][iHalo];
    }
    
    /// Surface site of a site in the e
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static auto surfSiteOfEdgeSite(const int64_t& iEdge)
    {
      assertHasDefinedCoverage();
      
      if constexpr(fieldCoverage==FULL_SPACE)
	return surflxOfEdgelx[iEdge];
      else
	if constexpr(fieldCoverage==EVEN_SITES or fieldCoverage==ODD_SITES)
	  return surfeo_of_edgeo[fieldCoverage][iEdge];
    }
    
#define PROVIDE_NEIGH(UD)						\
									\
    /* Neighbor in the UD orientation */				\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION				\
    static auto locNeigh ## UD(const int64_t& site,			\
			       const int& mu)				\
    {									\
      assertHasDefinedCoverage();					\
      									\
      if constexpr(fieldCoverage==FULL_SPACE)				\
	return loclxNeigh ## UD[site][mu];				\
      else								\
	if constexpr(fieldCoverage==EVEN_SITES or			\
		     fieldCoverage==ODD_SITES)				\
	  return loceo_neigh ## UD[fieldCoverage][site][mu];		\
    }
    
    PROVIDE_NEIGH(dw);
    
    PROVIDE_NEIGH(up);
    
#undef PROVIDE_NEIGH
    
    /////////////////////////////////////////////////////////////////
    
    /// Global coordinates
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static auto glbCoord(const int64_t& site,
			 const int& mu)
    {
      assertHasDefinedCoverage();
      
      if constexpr(fieldCoverage==FULL_SPACE)
	return glbCoordOfLoclx[site][mu];
      else
	if constexpr(fieldCoverage==EVEN_SITES or fieldCoverage==ODD_SITES)
	  return glbCoordOfLoclx[loclx_of_loceo[fieldCoverage][site]][mu];
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Reference to field
  template <typename F>
  struct FieldRef;
  
  /// Hosts all the data and provides low-level access
  template <typename D,
	    typename Comps,
	    MemorySpace MS>
  struct FieldData
  {
    /// Fundamental type
    using Fund=
      std::remove_all_extents_t<Comps>;
    
    /// Number of degrees of freedom
    static constexpr int nInternalDegs=
      sizeof(Comps)/sizeof(Fund);
    
    /// Size of the external nominal index
    int64_t externalSize;
    
    /// Pointer to the data
    Fund* _data;
    
    /// Calls the asserter to prove that we are running on the proper memory space
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    void assertMemorySpaceIs(const MemorySpace& DMS) const
    {
      if(DMS!=MS)
	CRASH("asked a pointer to be valid on memory space %s, but field is defined on %s",
	      memorySpaceName(DMS),memorySpaceName(MS));
    }
    
    /// Computes the index of the data
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    static int64_t index(const int64_t& site,
			 const int& internalDeg,
			 const int& externalSize)
    {
      if constexpr(D::spaceTimeLayout==SpaceTimeLayout::CPU)
	return internalDeg+nInternalDegs*site;
      else
	return site+externalSize*internalDeg;
    }
    
    /// Gets the pointer to the given combo of site/internalDeg
#define PROVIDE_GET_PTR(CONST)				\
    template <MemorySpace DMS>				\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION	\
    Fund* CONST& getPtr() CONST				\
    {							\
      assertMemorySpaceIs(DMS);				\
      							\
      return _data;					\
    }
    
    PROVIDE_GET_PTR(const);
    
    PROVIDE_GET_PTR(/* not const */);
    
#undef PROVIDE_GET_PTR
    
    /// Gets the pointer to the given combo of site/internalDeg
#define PROVIDE_GET_PTR_TO(CONST)			\
    template <MemorySpace DMS>				\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr	\
    CONST Fund* getPtrTo(const int64_t& site,		\
			 const int& internalDeg) CONST	\
    {							\
      return getPtr<DMS>()+index(site,internalDeg);	\
    }
    
    PROVIDE_GET_PTR_TO(const);
    
    PROVIDE_GET_PTR_TO(/* not const */);
    
    /// Nonstatic index
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    int64_t index(const int64_t& site,
		  const int& internalDeg) const
    {
      return index(site,internalDeg,externalSize);
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_FLATTENED_CALLER(CONST)					\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    CONST Fund& operator()(const int64_t& site,				\
			   const int& internalDeg) CONST		\
    {									\
      return *getPtrTo<currentMemorySpace>(site,internalDeg);		\
    }
    
    PROVIDE_FLATTENED_CALLER(const);
    
    PROVIDE_FLATTENED_CALLER(/* not const */);
    
#undef PROVIDE_FLATTENED_CALLER
    
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    decltype(auto) operator[](const int64_t& site) CONST		\
    {									\
      if constexpr(not std::is_array_v<Comps>)				\
	return								\
	  (*this)(site,0);						\
      else								\
	{								\
	  CONST D& self=						\
	    *static_cast<CONST D*>(this);				\
	  								\
	  if constexpr(D::spaceTimeLayout==SpaceTimeLayout::CPU)	\
	    return							\
	      ((CONST Comps*)self.template getPtr<currentMemorySpace>())[site]; \
	  else								\
	    return							\
	      SubscribedField<CONST D,					\
			      std::remove_extent_t<Comps>>(self,site,nullptr); \
	}								\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* not const */);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    /// Default Constructor
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    FieldData() :
      externalSize(0),
      _data(nullptr)
    {
    }
  };
  
  /// Field
  template <typename F>
  struct FieldRef :
    FieldData<FieldRef<F>,typename F::Comps,F::memorySpace>
  {
    /// Type describing the components
    using Comps=
      typename F::Comps;
    
    /// Type of the fundamental data
    using Fund=
      ConstIf<std::is_const_v<F>,typename F::Fund>;
    
    /// Memory space of the field
    static constexpr MemorySpace memorySpace=
      F::memorySpace;
    
    /// Spacetime layout of the field
    static constexpr SpaceTimeLayout spaceTimeLayout=
      F::spaceTimeLayout;
    
    /////////////////////////////////////////////////////////////////
    
    /// Pointer to the original field
    F* fptr;
    
    /// Set from f
    void setFrom(F& f)
    {
      if(fptr!=nullptr)
	CRASH("Unable to set an already set ref");
      
      fptr=&f;
      this->externalSize=f.externalSize;
      this->template getPtr<F::memorySpace>()=f.template getPtr<F::memorySpace>();
      fptr->nRef++;
    }
    
    /// Default constructor
    FieldRef()=default;
    
    /// Construct from field
    FieldRef(F& f) :
      fptr{nullptr}
    {
      setFrom(f);
    }
    
    /// Copy constructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    FieldRef(const FieldRef& oth) :
      fptr(oth.fptr)
    {
      this->externalSize=oth.externalSize;
      this->_data=oth._data;
      
#ifndef COMPILING_FOR_DEVICE
      fptr->nRef++;
#endif
    }
    
    /// Removes the reference
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    void unset()
    {
      if(fptr==nullptr)
	CRASH("Unable to unset a non set ref");
      
#ifndef COMPILING_FOR_DEVICE
      fptr->nRef--;
      if(fptr->nRef==0)
	if constexpr(not std::is_const_v<F>)
	  fptr->invalidateHalo();
#endif
      fptr=nullptr;
      this->_data=nullptr;
      this->externalSize=0;
    }
    
    /// Destructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    ~FieldRef()
    {
      unset();
    }
    
    /// Return itself as a readable
    constexpr INLINE_FUNCTION
    const FieldRef& getReadable() const
    {
      return *this;
    }
    
    /// Return itself as a writable
    constexpr INLINE_FUNCTION
    FieldRef& getWritable()
    {
      return *this;
    }
    
    /// Returns the same or a copy if not with the given spacetime layout
    template <SpaceTimeLayout Dest>
    decltype(auto) getSurelyWithSpaceTimeLayout() const
    {
      return fptr->template getSurelyWithSpaceTimeLayout<Dest>();
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  PROVIDE_FEATURE(Field);
  
  /// Field
  template <typename T,
	    FieldCoverage FC,
	    SpaceTimeLayout STL=defaultSpaceTimeLayout,
	    MemorySpace MS=defaultMemorySpace>
  struct Field :
    FieldData<Field<T,FC,STL,MS>,T,MS>,
    FieldFeat<Field<T,FC,STL,MS>>,
    FieldSizes<FC>
  {
    using FD=FieldData<Field<T,FC,STL,MS>,T,MS>;
    
    /// Coefficent which divides the space time, if the field is covering only half the space
    static constexpr const int divCoeff=
      (FC==FULL_SPACE)?1:2;
    
    mutable int nRef;
    
    /// Name of the field
    const char* name;
    
    /// Components
    using Comps=T;
    
    /// Coverage of sites
    static constexpr FieldCoverage fieldCoverage=FC;
    
    /// Arrangement of internal DOF w.r.t spacetime
    static constexpr SpaceTimeLayout spaceTimeLayout=STL;
    
    /// Memory space where the data is stored
    static constexpr MemorySpace memorySpace=MS;
    
    /// Presence of halo and edges
    const HaloEdgesPresence haloEdgesPresence;
    
    /// Import from FieldData
    static constexpr int nInternalDegs=
      FD::nInternalDegs;
    
    /// Import Fund from FieldData
    using typename FD::Fund;
    
    /// Store the data in case of a backup
    Field* backup;
    
    /// States whether the halo is updated
    mutable bool haloIsValid;
    
    /// States whether the edges are updated
    mutable bool edgesAreValid;
    
    // /// Exec the operation f on each site and degree of freedom
    // template <typename F>
    // Field& forEachSiteDeg(const F& f)
    // {
    //   PAR(0,this->nSites(),
    // 	  CAPTURE(f,t=this->getWritable()),site,
    // 	  {
    // 	    for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
    // 	      f(t(site,internalDeg),site,internalDeg);
    // 	  });
      
    //   invalidateHalo();
      
    //   return *this;
    // }
    
#define FOR_EACH_SITE_DEG_OF_FIELD(F,					\
				   CAPTURES,				\
				   SITE,				\
				   IDEG,				\
				   CODE...)				\
    PAR(0,								\
	(F).nSites(),							\
	CAPTURE(CAPTURES),						\
	SITE,								\
	{								\
	  constexpr int nInternalDegs=					\
	   std::remove_reference_t<decltype(F)>::nInternalDegs;		\
	 UNROLL_FOR(IDEG,0,nInternalDegs)				\
	   CODE								\
	})
    
    // /// Exec the operation f on each site
    // template <typename F>
    // Field& forEachSite(F&& f)
    // {
    //   PAR(site,0,this->nSites())
    // 	f((*this)[site]);
    //   NISSA_PARALLEL_LOOP_END;
      
    //   invalidateHalo();
      
    //   return *this;
    // }
    
#define PROVIDE_SELFOP(OP)						\
    Field& operator OP ## =(const Field& oth)				\
    {									\
      PAR(0,this->nSites(),						\
	  CAPTURE(t=this->getWritable(),				\
		  TO_READ(oth)),					\
	  site,								\
	  {								\
	    UNROLL_FOR(internalDeg,0,nInternalDegs)			\
	      t(site,internalDeg) OP ## =oth(site,internalDeg);		\
	  });								\
									\
      return *this;							\
    }
    
    PROVIDE_SELFOP(+);
    PROVIDE_SELFOP(-);
    
#undef PROVIDE_SELFOP
    
#define PROVIDE_SELF_SCALOP(OP)						\
    Field& operator OP ## =(const Fund& oth)				\
    {									\
      PAR(0,this->nSites(),						\
	  CAPTURE(oth,							\
		  t=this->getWritable()),				\
	  site,								\
	  {								\
	    UNROLL_FOR(internalDeg,0,nInternalDegs)			\
	      t(site,internalDeg) OP ## =oth;				\
	  });								\
									\
      return *this;							\
    }
    
    PROVIDE_SELF_SCALOP(*);
    PROVIDE_SELF_SCALOP(/);
    
#undef PROVIDE_SELF_SCALOP
    
    /// Reset to 0
    void reset()
    {
      PAR(0,this->nSites(),
	  CAPTURE(t=this->getWritable()),
	  site,
	  {
	    UNROLL_FOR(internalDeg,0,nInternalDegs)
	      t(site,internalDeg)=0.0;
	  });
    }
    
    /// Squared norm
    Fund norm2() const
    {
      return realPartOfScalarProdWith(*this);
    }
    
    /// Put the field to the new norm, returning the reciprocal of normalizing factor
    Fund normalize(const Fund& newNorm=1.0)
    {
      const Fund f=
	newNorm/sqrt(this->norm2());
      
      (*this)*=f;
      
      return 1/f;
    }
    
    /// Put the field to the given norm
    Fund normalize(const Field& oth,
		   const Fund& newNorm=1.0)
    {
      //compute current norm
      const Fund oldNorm2=
	sqrt(oth.norm2());
      
      //compute normalizing factor
      const Fund fact=
	newNorm/sqrt(oldNorm2);
      
      PAR(0,this->nSites(),
	  CAPTURE(t=this->getWritable(),
		  TO_READ(oth),
		  fact),
	  site,
	  {
	    UNROLL_FOR(internalDeg,0,nInternalDegs)
	      t(site,internalDeg)=oth(site,internalDeg)*fact;
	  });
      
      return 1/fact;
    }
    
    /// Reduce the field summing all elements and storing the result inside 'out'
    template <typename F=GlbReduceSumFunctor>
    void reduceNoPreserve(T& out,
			  const F& f=GlbReduceSumFunctor(),
			  const MPI_Op& mpiOp=MPI_Op_sum_for_type<T>())
    {
      Field& self=*this;
      
      int64_t n=
	this->nSites();
      
      // const double init_time=take_time();
      while(n>1)
	{
	  const int64_t stride=(n+1)/2;
	  const int64_t nReductions=n/2;
	  
	  PAR(0,nReductions,
	      CAPTURE(stride,
		      TO_WRITE(self),
		      f),
		      ireduction,
	      {
		const int64_t first=ireduction;
		const int64_t second=first+stride;
		
		UNROLL_FOR(iDeg,0,nInternalDegs)
		  f(self(first,iDeg),self(second,iDeg));
	      });
	  n=stride;
	};
      
      UNROLL_FOR(iDeg,0,nInternalDegs)
#if defined(USE_CUDA)
 	if constexpr(memorySpace==MemorySpace::GPU)
	  cudaMemcpy((Fund*)&out+iDeg,
		     self.template getPtrTo<MemorySpace::GPU>(0,iDeg),
		     sizeof(Fund),
		     cudaMemcpyDeviceToHost);
	else
#endif
	  ((Fund*)&out)[iDeg]=self(0,iDeg);
      
      MPI_Allreduce(MPI_IN_PLACE,&out,nInternalDegs,MPI_Datatype_of<Fund>(),mpiOp,MPI_COMM_WORLD);
    }
    
    /// Reduce the field summing all elements and storing the result inside 'out'
    template <typename F=GlbReduceSumFunctor>
    void reduce(T& out,
		const F& f=GlbReduceSumFunctor(),
		const MPI_Op& mpiOp=MPI_Op_sum_for_type<T>()) const
    {
      /// Copy to a temporary to avoid destroying the field
      Field tmp("tmp");
      tmp=*this;
      
      tmp.reduceNoPreserve(out,f,mpiOp);
    }
    
    /// Reduce with accuracy
    template <typename F=GlbReduceSumFunctor>
    void preciseReduce(T& out,
		       const F& f=GlbReduceSumFunctor(),
		       const MPI_Op& mpiOp=MPI_Op_sum_for_type<Float128>()) const
    {
      using T128=
	DuplicateExtents<T,Float128>;
      
      /// Copy to a high prec field
      Field<T128,FC,STL,MS> tmp("tmp");
      FOR_EACH_SITE_DEG_OF_FIELD(tmp,
				 CAPTURE(TO_WRITE(tmp),
					 self=this->getReadable()),
				 site,
				 iDeg,
				 {
				   tmp(site,iDeg)=self(site,iDeg);
				 });
      
      T128 tmpOut;
      tmp.reduceNoPreserve(tmpOut,f,mpiOp);
      
      for(int i=0;i<nInternalDegs;i++)
	((Fund*)&out)[i]=((Float128*)&tmpOut)[i].roundDown();
    }
    
    /// (*this,out)
    void scalarProdWith(complex& out,
			const Field& oth) const
    {
      Field<complex,FC> buf("buf");
      
      using NT=Fund[nInternalDegs][2];
      const auto& l=castComponents<NT>();
      const auto& r=oth.castComponents<NT>();
      
      PAR(0,this->nSites(),
	  CAPTURE(TO_WRITE(buf),
		  TO_READ(r),
		  TO_READ(l)),
	  site,
	  {
	    complex c;
	    complex_put_to_zero(c);
	    UNROLL_FOR(internalDeg,0,nInternalDegs)
	      complex_summ_the_conj1_prod(c,l[site][internalDeg],r[site][internalDeg]);
	    complex_copy(buf[site],c);
	  });
      
      complex res;
      glb_reduce(&res,buf,this->nSites());
    }
    
    /// Re(*this,out)
    Fund realPartOfScalarProdWith(const Field& oth) const
    {
      Fund res;
      
// #ifdef USE_CUDA
//       if(spaceTimeLayout==SpaceTimeLayout::CPU or bord_vol==0 or haloEdgesPresence==WITHOUT_HALO)
// 	{
// 	  res=
// 	    thrust::inner_product(thrust::device,
// 				  _data,
// 				  _data+this->nSites()*nInternalDegs,
// 				  oth._data,
// 				  Fund{});
// 	  non_loc_reduce(&res);
// 	}
//       else
// 	{
// #endif
	  Field<Fund,FC> buf("buf");
	  
	  /// Parsing is a bit failing in some version of cuda
	  const auto& self=
	    *this;
	  
	  PAR(0,oth.nSites(),
	      CAPTURE(TO_WRITE(buf),
		      TO_READ(self),
		      TO_READ(oth)),
	      site,
	      {
		Fund r{};
		UNROLL_FOR(internalDeg,0,nInternalDegs)
		  r+=self(site,internalDeg)*oth(site,internalDeg);
		buf[site]=r;
	      });
	  
	  buf.reduce(res);
	  
// #ifdef USE_CUDA
// 	}
// #endif
      
      return res;
    }
    
#define PROVIDE_CASTS(CONST)						\
    /* Cast to a different fieldCoverage */				\
    template <FieldCoverage NFC,					\
	      bool Force=false>						\
    CONST Field<T,NFC,STL>& castFieldCoverage() CONST			\
    {									\
      if constexpr(not (Force or FC== EVEN_OR_ODD_SITES))		\
	static_assert(NFC==EVEN_SITES or NFC==ODD_SITES,		\
		      "incompatible fieldCoverage! Force the change if needed"); \
      									\
      return *(CONST Field<T,NFC,STL>*)this;				\
    }									\
									\
    /* Cast to a different comps */					\
    template <typename NT,						\
	      bool Force=false>						\
    CONST Field<NT,FC,STL>& castComponents() CONST			\
    {									\
      static_assert(Force or sizeof(T)==sizeof(NT),			\
		    "incompatible components! Force the change if needed"); \
      									\
      return *(CONST Field<NT,FC,STL>*)this;				\
    }
    
    PROVIDE_CASTS(const);
    
    PROVIDE_CASTS(/* not const */);
    
#undef PROVIDE_CASTS
    
    /// Returns a copy cast to the given fund
    template <typename NewFund>
    auto fundCast() const
    {
      /// Type to be used for the new underlying struct
      using NewT=
	DuplicateExtents<T,NewFund>;
      
      Field<NewT,FC,STL> res(name,haloEdgesPresence);
      res=*this;
      
      return res;
    }
    
    /// Constructor
    Field(const char *name,
	  const HaloEdgesPresence& haloEdgesPresence=WITHOUT_HALO) :
      nRef(0),
      name(name),
      haloEdgesPresence(haloEdgesPresence),
      haloIsValid(false),
      edgesAreValid(false)
    {
      this->externalSize=FieldSizes<fieldCoverage>::nSitesToAllocate(haloEdgesPresence);
      VERBOSITY_LV3_MASTER_PRINTF("Allocating field %s\n",name);
      this->_data=memoryManager<MS>()->template provide<Fund>(this->externalSize*nInternalDegs);
    }
    
    /// Construct from other layout
    template <SpaceTimeLayout OFL,
	      ENABLE_THIS_TEMPLATE_IF(OFL!=STL)>
    Field(const Field<T,FC,OFL>& oth) :
      Field(oth.name,oth.haloEdgesPresence)
    {
      *this=oth;
    }
    
    /// Move constructor
    Field(Field&& oth) :
      nRef(oth.nRef),
      name(oth.name),
      haloEdgesPresence(oth.haloEdgesPresence),
      haloIsValid(oth.haloIsValid),
      edgesAreValid(oth.edgesAreValid)
    {
      this->externalSize=oth.externalSize;
      this->_data=oth._data;
      oth.nRef=0;
      oth._data=nullptr;
    }
    
    /// Copy constructor
    Field(const Field& oth) :
      Field(oth.name,oth.haloEdgesPresence)
    {
      *this=oth;
    }
    
    /// Destructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    ~Field()
    {
#ifndef COMPILING_FOR_DEVICE
      VERBOSITY_LV3_MASTER_PRINTF("Deallocating field %s\n",name);
      if(this->_data)
	{
	  if(nRef==0)
	    memoryManager<MS>()->release(this->_data);
	  else
	    CRASH("Trying to destroying field %s with dangling references",name);
	}
#endif
    }
    
    /// Register the actions to backup and restore the field
    constexpr INLINE_FUNCTION
    FieldRef<Field> getWritable()
    {
#ifdef USE_CUDA
      if(insideParallelFor)
	{
	  benchmarkBeginActions.emplace_back([this]()
	  {
	    VERBOSITY_LV3_MASTER_PRINTF("Backing up field %s\n",name);
	    backup=new Field("backup",haloEdgesPresence);
	    *backup=*this;
	  });
	  
	  benchmarkEndActions.emplace_back([this]()
	  {
	    *this=*backup;
	    delete backup;
	    VERBOSITY_LV3_MASTER_PRINTF("Restored field %s\n",name);
	  });
	  VERBOSITY_LV3_MASTER_PRINTF("Added backup actions for field %s\n",name);
	}
      else
	VERBOSITY_LV3_MASTER_PRINTF("Skipping backup actions for field %s as we are not inside a parallel for\n",name);
#endif
      
      return *this;
    }
    
    /// Gets a readable view
    constexpr INLINE_FUNCTION
    FieldRef<const Field> getReadable() const
    {
      return *this;
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Set edges as invalid
    void invalidateEdges()
    {
      edgesAreValid=false;
    }
    
    /// Set halo as invalid
    void invalidateHalo()
    {
      haloIsValid=false;
      
      invalidateEdges();
    }
    
    /// Fill the sending buf using the data with a given function
    template <typename F>
    void fillSendingBufWith(const F& f,
			    const int& n) const
    {
      PAR(0,n,
	  CAPTURE(f,
		  sb=sendBuf,
		  t=this->getReadable()),
	  i,
	  {
	    UNROLL_FOR(internalDeg,0,nInternalDegs)
	      ((Fund*)sb)[internalDeg+nInternalDegs*i]=
		t(f(i),internalDeg);
	  });
    }
    
    /// Fill the sending buf using the data on the surface of a field
    void fillSendingBufWithSurface() const
    {
      fillSendingBufWith(LAMBDA_FUNCTION_FOR_DEFAULT_THREADS(INLINE_ATTRIBUTE,
							     CAPTURE(),
							     i,
							     {
							       return Field::surfSiteOfHaloSite(i);
							     }),
			 bordVol/divCoeff);
    }
    
    /// Fill the sending buf using the data on the surface edge
    void fillSendingBufWithEdgesSurface() const
    {
      fillSendingBufWith(LAMBDA_FUNCTION_FOR_DEFAULT_THREADS(INLINE_ATTRIBUTE,
							     CAPTURE(),
							     i,
							     {
							       return Field::surfSiteOfEdgeSite(i);
							     }),edgeVol/divCoeff);
    }
    
    /// Fill the surface using the data from the buffer
    template <typename B,
	      typename F>
    void fillSurfaceWithReceivingBuf(const F& f)
    {
      for(int bf=0;bf<2;bf++)
	for(int mu=0;mu<NDIM;mu++)
	  PAR(0,bordDirVol[mu]/divCoeff,
	      CAPTURE(f,bf,mu,
		      t=this->getWritable()),
	      iHaloOriDir,
	      {
		const int iHalo=
		  bf*bordVolh/divCoeff+
		  iHaloOriDir+bordOffset[mu]/divCoeff;
		
		const int iSurf=
		  Field::surfSiteOfHaloSite(iHalo);
		
		f(t[iSurf],
		  ((B*)recvBuf)[iHalo],
		  bf,
		  mu);
	      });
    }
    
    /// Fill the sending buf using the data with a given function
    void fillHaloOrEdgesWithReceivingBuf(const int& offset,
					 const int& n) const
    {
      PAR(0,n,
	  CAPTURE(data=this->template getPtr<defaultMemorySpace>(),
		  offset,
		  rb=recvBuf,
		  defaultMemorySpace=defaultMemorySpace,
		  externalSize=this->externalSize),
	  i,
	  {
	    assertRunningOnMemorySpace(defaultMemorySpace);
	    
	    for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	      data[FD::index(offset+i,internalDeg,externalSize)]=
		((Fund*)rb)[internalDeg+nInternalDegs*i];
	  });
    }
    
    /// Fills the halo with the received buffer
    void fillHaloWithReceivingBuf() const
    {
      assertHasHalo();
      
      fillHaloOrEdgesWithReceivingBuf(locVol/divCoeff,bordVol/divCoeff);
    }
    
    /// Fills the halo with the received buffer
    void fillEdgesWithReceivingBuf() const
    {
      assertHasEdges();
      
      fillHaloOrEdgesWithReceivingBuf((locVol+bordVol)/divCoeff,edgeVol/divCoeff);
    }
    
    /// Fills the sending buffer with the halo, compressing into elements of B using f
    template <typename B,
	      typename F>
    void fillSendingBufWithHalo(const F& f) const
    {
      for(int bf=0;bf<2;bf++)
	for(int mu=0;mu<NDIM;mu++)
	  PAR(0,bordDirVol[mu]/divCoeff,
	      CAPTURE(bf,mu,f,
		      t=this->getReadable()),
	      iHaloOriDir,
	      {
		const int iHalo=
		  bf*bordVolh/divCoeff+
		  iHaloOriDir+bordOffset[mu]/divCoeff;
		
		f(((B*)sendBuf)[iHalo],
		  t[locVol/divCoeff+iHalo],
		  bf,
		  mu);
	      });
    }
    
    /// Start the communications of halo
    static std::vector<MPI_Request> startHaloAsyincComm()
    {
      return startBufHaloNeighExchange<T>(divCoeff);
    }
    
    /// Start the communications of edges
    static std::vector<MPI_Request> startEdgesAsyincComm()
    {
      return startBufEdgesNeighExchange<T>(divCoeff);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Assert that can communicate n sites
    static void assertCanCommunicate(const int& n)
    {
      /// Needed size in the buffer
      const size_t neededBufSize=
	sizeof(T)*n;
      
      const size_t maxBufSize=
	std::min(sendBufSize,recvBufSize);
      
      if(neededBufSize>maxBufSize)
	CRASH("asking to create a communicator that needs %ld large buffer (%ld allocated)",
		  neededBufSize,maxBufSize);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Start communication using halo
    std::vector<MPI_Request> startCommunicatingHalo() const
    {
      /// Pending requests
      std::vector<MPI_Request> requests;
      
      assertCanCommunicate(Field::nHaloSites());
      
      //take time and write some debug output
      START_TIMING(tot_comm_time,ntot_comm);
      VERBOSITY_LV3_MASTER_PRINTF("Start communication of borders of %s\n",name);
      
      //fill the communicator buffer, start the communication and take time
      fillSendingBufWithSurface();
      requests=startHaloAsyincComm();
      STOP_TIMING(tot_comm_time);
      
      return requests;
    }
    
    /// Finalize communications
    void finishCommunicatingHalo(std::vector<MPI_Request> requests) const
    {
      //take note of passed time and write some debug info
      START_TIMING(tot_comm_time,ntot_comm);
      VERBOSITY_LV3_MASTER_PRINTF("Finish communication of halo of %s\n",name);
      
      //wait communication to finish, fill back the vector and take time
      waitAsyncCommsFinish(requests);
      fillHaloWithReceivingBuf();
      STOP_TIMING(tot_comm_time);
      
      haloIsValid=true;
    }
    
    /// Crash if the halo is not allocated
    void assertHasHalo() const
    {
      if(not (haloEdgesPresence>=WITH_HALO))
	CRASH("needs halo allocated for %s!",name);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Start communication using edges
    std::vector<MPI_Request> startCommunicatingEdges() const
    {
      /// Pending requests
      std::vector<MPI_Request> requests;
      
      assertCanCommunicate(Field::nEdgesSites());
      
      //take time and write some debug output
      START_TIMING(tot_comm_time,ntot_comm);
      VERBOSITY_LV3_MASTER_PRINTF("Start communication of edges of %s\n",name);
      
      //fill the communicator buffer, start the communication and take time
      fillSendingBufWithEdgesSurface();
      requests=startEdgesAsyincComm();
      STOP_TIMING(tot_comm_time);
      
      return requests;
    }
    
    /// Finalize communications
    void finishCommunicatingEdges(std::vector<MPI_Request> requests) const
    {
      //take note of passed time and write some debug info
      START_TIMING(tot_comm_time,ntot_comm);
      VERBOSITY_LV3_MASTER_PRINTF("Finish communication of edges of %s\n",name);
      
      //wait communication to finish, fill back the vector and take time
      waitAsyncCommsFinish(requests);
      fillEdgesWithReceivingBuf();
      STOP_TIMING(tot_comm_time);
      
      haloIsValid=true;
    }
    
    /// Crash if the edges are not allocated
    void assertHasEdges() const
    {
      if(not (haloEdgesPresence>=WITH_HALO_EDGES))
	CRASH("needs edges allocated on field %s!",name);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Communicate the halo
    void updateHalo(const bool& force=false) const
    {
      if(force or not haloIsValid)
	{
	  VERBOSITY_LV3_MASTER_PRINTF("Sync communication of halo of %s\n",name);
	  
	  const std::vector<MPI_Request> requests=
	    startCommunicatingHalo();
	  finishCommunicatingHalo(requests);
      }
    }
    
    /// Communicate the edges
    void updateEdges(const bool& force=false) const
    {
      updateHalo(force);
      
      if(force or not edgesAreValid)
	{
	  VERBOSITY_LV3_MASTER_PRINTF("Sync communication of edges of %s\n",name);
	  
	  const std::vector<MPI_Request> requests=
	    startCommunicatingEdges();
	  finishCommunicatingEdges(requests);
      }
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Compare
    INLINE_FUNCTION
    bool operator==(const Field& oth) const
    {
      return this->_data==oth._data;
    }
    
    /// Negate comparison
    INLINE_FUNCTION
    bool operator!=(const Field& oth) const
    {
      return not ((*this)==oth);
    }
    
    /// Assigns with no check
    template <typename O>
    INLINE_FUNCTION
    void assign(const O& oth)
    {
      const bool b=doNotBackupDuringBenchmark;
      doNotBackupDuringBenchmark=true;
      
      /// Parsing is a bit failing in some version of cuda
      auto& self=
	*this;
      
      FOR_EACH_SITE_DEG_OF_FIELD(self,
				 CAPTURE(TO_READ(oth),
					 TO_WRITE(self)),
				 site,
				 iDeg,
				 {
				   self(site,iDeg)=oth(site,iDeg);
				 }
      );
      
      doNotBackupDuringBenchmark=b;
    }
    
    /// Assigns from a different layout
    template <typename OT,
	      SpaceTimeLayout OFl,
	      ENABLE_THIS_TEMPLATE_IF(std::is_same_v<T,OT> or
				      std::is_same_v<T,DuplicateExtents<OT,Fund>>)>
    INLINE_FUNCTION
    Field& operator=(const Field<OT,FC,OFl,MS>& oth)
    {
      assign(oth);
      
      return *this;
    }
    
    /// Assigns from the same layout but a different memory space
    template <MemorySpace OMS>
    INLINE_FUNCTION
    Field& operator=(const Field<T,FC,STL,OMS>& oth)
    {
      if(oth.haloEdgesPresence!=haloEdgesPresence)
	CRASH("can copy across memory spaces only if the edges/halo are present equally on source and destination");
#ifdef USE_CUDA
      const double initTime=take_time();
      const size_t size=this->externalSize*nInternalDegs*sizeof(Fund);
      cudaMemcpy(this->template getPtr<MS>(),
		 oth.template getPtr<OMS>(),
		 size,
		 memcpyKindForCopy<MS,OMS>);
      const double usedTime=take_time()-initTime;
      MASTER_PRINTF("Bare time for cudaMemcpy of %zu bytes: %lg s, %lg GB/s\n",size,usedTime,size/usedTime*1e-9);
#else
      CRASH("Not compiled with cuda");
#endif
      
      return *this;
    }
    
    /// Returns the same or a copy if not on the given memory space
    template <MemorySpace Dest>
    decltype(auto) getSurelyReadableOn() const
    {
      if constexpr(Dest!=MS)
	{
	  const double tin=take_time();
	  Field<T,FC,STL,Dest> out("out",haloEdgesPresence);
	  out=*this;
	  
	  MASTER_PRINTF("Time to copy %zu bytes: %lg s\n",sizeof(Fund)*this->externalSize*nInternalDegs,take_time()-tin);
	  return out;
	}
      else
	return this->getReadable();
    }
    
    /// Make it possible to write to the field, and pass it to the function f
    template <MemorySpace Dest,
	      typename F>
    void passSurelyWritableOn(F f,
			      const bool& uninited)
    {
      if constexpr(Dest!=MS)
	{
	  Field<T,FC,STL,Dest> tmp("tmp",haloEdgesPresence);
	  
	  if(not uninited)
	    tmp=*this;
	  
	  f(tmp);
	  
	  *this=tmp;;
	}
      else
	return f(*this);
    }
    
    /// Make it possible to write to the field, and pass it to the function f, ignoring original status
    template <MemorySpace Dest,
	      typename F>
    void initOn(F f)
    {
      passSurelyWritableOn<Dest>(f,true);
    }
    
    /// Returns the same or a copy if not with the given spacetime layout
    template <SpaceTimeLayout Dest>
    decltype(auto) getSurelyWithSpaceTimeLayout() const
    {
      if constexpr(Dest!=STL)
	{
	  Field<T,FC,Dest,MS> out("out",haloEdgesPresence);
	  out=*this;
	  return out;
	}
      else
	return *this;
    }
    
    /// Assign the same layout
    INLINE_FUNCTION
    Field& operator=(const Field& oth)
    {
      if(this!=&oth)
	assign(oth);
      
      return *this;
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Hack
  template <typename F>
  void set_borders_invalid(FieldFeat<F>& field)
  {
    static_cast<F*>(&field)->invalidateHalo();
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Lexicographic field
  template <typename T,
	    SpaceTimeLayout STL=defaultSpaceTimeLayout,
	    MemorySpace MS=defaultMemorySpace>
  using LxField=Field<T,FULL_SPACE,STL,MS>;
  
  /// Field over even sites
  template <typename T,
	    SpaceTimeLayout STL=defaultSpaceTimeLayout,
	    MemorySpace MS=defaultMemorySpace>
  using EvnField=Field<T,EVEN_SITES,STL,MS>;
  
  /// Field over odd sites
  template <typename T,
	    SpaceTimeLayout STL=defaultSpaceTimeLayout,
	    MemorySpace MS=defaultMemorySpace>
  using OddField=Field<T,ODD_SITES,STL,MS>;
  
  /// Field over even or odd sites
  template <typename T,
	    SpaceTimeLayout STL=defaultSpaceTimeLayout,
	    MemorySpace MS=defaultMemorySpace>
  using EvenOrOddField=Field<T,EVEN_OR_ODD_SITES,STL,MS>;
  
  /////////////////////////////////////////////////////////////////
  
  template <int EO>
  struct Par
  {
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    operator int() const
    {
      return EO;
    }
    
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    int operator()() const
    {
      return EO;
    }
    
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Par<1-EO> operator!() const
    {
      return {};
    }
  };
  
  /// Structure to hold an even/old field
  template <typename T,
	    SpaceTimeLayout STL=defaultSpaceTimeLayout,
	    MemorySpace MS=defaultMemorySpace,
	    typename Fevn=Field<T,EVEN_SITES,STL,MS>,
	    typename Fodd=Field<T,ODD_SITES,STL,MS>,
	    typename FevnOrOdd=Field<T,EVEN_OR_ODD_SITES,STL,MS>>
  struct EoField
  {
    /// Type representing a pointer to type T
    template <FieldCoverage EO>
    using F=
      std::conditional_t<EO==EVN,Fevn,Fodd>;
    
    /// Fundamental type
    using Fund=
      typename FevnOrOdd::Fund;
    
    Fevn evenPart;
    
    Fodd oddPart;
    
    /////////////////////////////////////////////////////////////////
    
    /// Gets the even or odd part
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    template <int EO>							\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr			\
    CONST auto& operator[](Par<EO>) CONST				\
    {									\
      if constexpr(EO==EVN)						\
	return evenPart;						\
      else								\
	return oddPart;							\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* const*/ );
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    /////////////////////////////////////////////////////////////////
    
    /// Gets the even or odd part, non-constexpr version
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION				\
    CONST auto& operator[](const int eo) CONST				\
    {									\
      using EOOF=							\
	CONST FevnOrOdd;						\
      									\
      EOOF* t[2]={(EOOF*)&evenPart,(EOOF*)&oddPart};			\
      									\
      return *t[eo];							\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* const*/ );
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    /////////////////////////////////////////////////////////////////
    
    constexpr INLINE_FUNCTION
    EoField<T,STL,MS,
	    FieldRef<Field<T,EVEN_SITES,STL,MS>>,
	    FieldRef<Field<T,ODD_SITES,STL,MS>>,
	    FieldRef<Field<T,EVEN_OR_ODD_SITES,STL,MS>>>
    getWritable()
    {
      return {evenPart,oddPart};
    }
    
    constexpr INLINE_FUNCTION
    EoField<T,STL,MS,
	    FieldRef<const Field<T,EVEN_SITES,STL,MS>>,
	    FieldRef<const Field<T,ODD_SITES,STL,MS>>,
	    FieldRef<const Field<T,EVEN_OR_ODD_SITES,STL,MS>>>
    getReadable() const
    {
      return {evenPart,oddPart};
    }
    
    /// Constructor
    EoField(Fevn&& ev,
	    Fodd&& od) :
      evenPart(ev),
      oddPart(od)
    {
    }
    
    /// Constructor
    EoField(const char* name,
	    const HaloEdgesPresence& haloEdgesPresence=WITHOUT_HALO) :
      evenPart(name,haloEdgesPresence),
      oddPart(name,haloEdgesPresence)
    {
    }
    
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    ~EoField()
    {
    }
    
    /// Reset
    INLINE_FUNCTION
    void reset()
    {
      forBothParities([this](const auto& par)
      {
	(*this)[par].reset();
      });
    }
    
    /// Update the halo of both parities
    INLINE_FUNCTION
    void updateHalo(const bool& force=0) const
    {
      forBothParities([force,
		       this](const auto& par)
      {
	(*this)[par].updateHalo(force);
      });
    }
    
    /// Update the edges of both parities
    INLINE_FUNCTION
    void updateEdges() const
    {
      forBothParities([this](const auto& par)
      {
	(*this)[par].updateEdges();
      });
    }
    
    /// Invalidate the halo of both parities
    INLINE_FUNCTION
    void invalidateHalo()
    {
      forBothParities([this](const auto& par)
      {
	(*this)[par].invalidateHalo();
      });
    }
    
    /// Invalidate the edges of both parities
    INLINE_FUNCTION
    void invalidateEdges()
    {
      forBothParities([this](const auto& par)
      {
	(*this)[par].invalidateEdges();
      });
    }
    
#define PROVIDE_SELFOP(OP)						\
    EoField& operator OP ## =(const EoField& oth)			\
    {									\
      evenPart OP ## =oth.evenPart;					\
      oddPart OP ## =oth.oddPart;					\
									\
      return *this;							\
    }
    
    PROVIDE_SELFOP(+);
    PROVIDE_SELFOP(-);
    
#undef PROVIDE_SELFOP
    
#define PROVIDE_SELF_SCALOP(OP)						\
    EoField& operator OP ## =(const Fund& oth)				\
    {									\
      evenPart OP ## =oth.evenPart;					\
      oddPart OP ## =oth.oddPart;					\
    }
    
    PROVIDE_SELF_SCALOP(*);
    PROVIDE_SELF_SCALOP(/);
    
#undef PROVIDE_SELF_SCALOP
    
    /// Returns a copy cast to the given fund
    template <typename NewFund>
    auto fundCast() const
    {
      using NewT=
	DuplicateExtents<T,NewFund>;
      
      return
	EoField<NewT,STL>(this->evenPart.template fundCast<NewFund>(),
			  this->oddPart.template fundCast<NewFund>());
    }
    
    /// Compare
    INLINE_FUNCTION
    bool operator==(const EoField& oth) const
    {
      return evenPart==oth.evenPart and oddPart==oth.oddPart;
    }
    
    /// Negate comparison
    INLINE_FUNCTION
    bool operator!=(const EoField& oth) const
    {
      return not ((*this)==oth);
    }
    
    /// Assign
    INLINE_FUNCTION
    EoField& operator=(const EoField& oth)
    {
      evenPart=oth.evenPart;
      oddPart=oth.oddPart;
      
      return *this;
    }
    
    /// Returns the squared norm
    Fund norm2() const
    {
      return
	evenPart.norm2()+
	oddPart.norm2();
    }
  };
  
  /// Loop on both parities
  template <typename F>
  INLINE_FUNCTION void forBothParities(F&& f)
  {
    f(Par<0>{});
    f(Par<1>{});
  }

#define FOR_PARITY(PAR,			       \
		   ID,			       \
		   CODE...)		       \
  {						\
    Par<ID> PAR;				\
  CODE						\
    }
  
#define FOR_BOTH_PARITIES(PAR,	       \
			  CODE...)     \
  FOR_PARITY(PAR,0,CODE);\
  FOR_PARITY(PAR,1,CODE)
  
  /////////////////////////////////////////////////////////////////
  
  template <typename F,
	    typename P>
  struct SubscribedField
  {
    F& f;
    
    using Fund=
      typename F::Fund;
    
    using Comps=
      typename F::Comps;
    
    const int64_t site;
    
    const P* ptr;
    
#define PROVIDE_EVAL(CONST)						\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr			\
    CONST ConstIf<std::is_const_v<std::remove_reference_t<F>>,Fund>& eval(const int& i) CONST \
    {									\
      const int internalDeg=						\
	(int)(size_t)(&ptr[i])/sizeof(Fund);				\
      									\
      return f(site,internalDeg);					\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/* not const */);
    
#undef PROVIDE_EVAL
    
#define PROVIDE_SELFOP(OP)						\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr			\
    SubscribedField& operator OP ## =(const SubscribedField& oth)	\
      {									\
	UNROLL_FOR(internalDeg,0,F::nInternalDegs)			\
	eval(internalDeg) OP ## =oth.eval(internalDeg);			\
									\
	return *this;							\
      }
    
    PROVIDE_SELFOP(+);
    PROVIDE_SELFOP(-);
    
#undef PROVIDE_SELFOP
    
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    decltype(auto) operator[](const int& i) CONST			\
    {									\
      if constexpr(std::is_array_v<P>)					\
	return								\
	  SubscribedField<F,std::remove_extent_t<P>>(f,site,ptr[i]);	\
      else								\
	return eval(i);							\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* not const */);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    SubscribedField(F& f,
		    const int64_t& site,
		    const P* ptr) :
      f(f),
      site(site),
      ptr(ptr)
    {
    }
  };
}

#endif
