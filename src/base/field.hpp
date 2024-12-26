#ifndef _FIELD_HPP
#define _FIELD_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <cstddef>
#include <mpi.h>
#include <optional>
#include <type_traits>

#include <base/bench.hpp>
#include <base/vectors.hpp>
#include <communicate/communicate.hpp>
#include <geometry/geometry_eo.hpp>
#include <geometry/geometry_lx.hpp>
#include <linalgs/reduce.hpp>
#include <metaprogramming/constnessChanger.hpp>
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
    verbosity_lv3_master_printf("Entering MPI comm wait\n");
    
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
  
  /// Memory layout
  enum class FieldLayout{CPU,GPU};
  
  /// Coverage of the field
  enum FieldCoverage{EVEN_SITES,ODD_SITES,FULL_SPACE,EVEN_OR_ODD_SITES};
  
  /// Has or not the halo and the edges
  enum HaloEdgesPresence{WITHOUT_HALO,WITH_HALO,WITH_HALO_EDGES};
  
  /// Predefinite memory layout
  constexpr FieldLayout defaultFieldLayout=
	      FieldLayout::
#ifdef USE_CUDA
	      GPU
#else
	      CPU
#endif
	      ;
  
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
    static constexpr int64_t& nSites()
    {
      if constexpr(fieldCoverage==FULL_SPACE)
	return locVol;
      else
	return locVolh;
    }
    
    /// Number of sites in the halo of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr int64_t& nHaloSites()
    {
      if constexpr(fieldCoverage==FULL_SPACE)
	return bord_vol;
      else
	return bord_volh;
    }
    
    /// Number of sites in the edges of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr int64_t& nEdgesSites()
    {
      if constexpr(fieldCoverage==FULL_SPACE)
	return edge_vol;
      else
	return edge_volh;
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
  
  /// Field
  template <typename F>
  struct FieldRef
  {
    using Comps=
      typename F::Comps;
    
    using Fund=
      ConstIf<std::is_const_v<F>,typename F::Fund>;
    
    static constexpr int nInternalDegs=
      F::nInternalDegs;
    
    F& f;
    
    const int64_t externalSize;
    
    Fund* _data;
    
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    decltype(auto) operator[](const int64_t& site) CONST		\
    {									\
      if constexpr(not std::is_array_v<Comps>)				\
	return								\
	  _data[site];							\
      else								\
	return								\
	  SubscribedField<CONST FieldRef,				\
	  std::remove_extent_t<Comps>>(*this,site,nullptr);		\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* not const */);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_FLATTENED_CALLER(CONST)					\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    CONST Fund& operator()(const int64_t& site,				\
			   const int& internalDeg) CONST		\
    {									\
      return _data[index(site,internalDeg)];				\
    }
    
    PROVIDE_FLATTENED_CALLER(const);
    
    PROVIDE_FLATTENED_CALLER(/* not const */);
    
#undef PROVIDE_FLATTENED_CALLER
    
    /////////////////////////////////////////////////////////////////
    
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    int64_t index(const int64_t& site,
		  const int& internalDeg) const
    {
      return F::index(site,internalDeg,externalSize);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Construct from field
    FieldRef(F& f) :
      f(f),
      externalSize(f.externalSize),
      _data(f._data)
    {
      f.nRef++;
    }
    
    /// Copy constructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    FieldRef(const FieldRef& oth) :
      f(oth.f),
      externalSize(oth.externalSize),
      _data(oth._data)
    {
#ifndef COMPILING_FOR_DEVICE
      f.nRef++;
#endif
    }
    
    /// Destructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    ~FieldRef()
    {
#ifndef COMPILING_FOR_DEVICE
      f.nRef--;
      if(f.nRef==0)
	if constexpr(not std::is_const_v<F>)
	  f.invalidateHalo();
#endif
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  PROVIDE_FEATURE(Field);
  
  /// Field
  template <typename T,
	    FieldCoverage FC,
	    FieldLayout FL=defaultFieldLayout>
  struct Field :
    FieldFeat<Field<T,FC,FL>>,
    FieldSizes<FC>
  {
    /// Coefficient which divides the space time, if the field is covering only half the space
    static constexpr const int divCoeff=
      (FC==FULL_SPACE)?1:2;
    
    mutable int nRef;
    
    /// Name of the field
    const char* name;
    
    /// Fundamental type
    using Fund=
      std::remove_all_extents_t<T>;
    
    /// Components
    using Comps=T;
    
    /// Coverage of sites
    static constexpr FieldCoverage fieldCoverage=FC;
    
    /// Memory layout of the field
    static constexpr FieldLayout fieldLayout=FL;
    
    /// Number of degrees of freedom
    static constexpr int nInternalDegs=
      sizeof(Comps)/sizeof(Fund);
    
    /// Presence of halo and edges
    const HaloEdgesPresence haloEdgesPresence;
    
    /// Total allocated sites
    const int64_t externalSize;
    
    /// Container for actual data
    Fund* _data;
    
    /// Store the data in case of a backup
    Field* backup;
    
    /// States whether the halo is updated
    mutable bool haloIsValid;
    
    /// States whether the edges are updated
    mutable bool edgesAreValid;
    
    /// Computes the index of the data
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    static int64_t index(const int64_t& site,
			 const int& internalDeg,
			 const int& externalSize)
    {
      if constexpr(FL==FieldLayout::CPU)
	return internalDeg+nInternalDegs*site;
      else
	return site+externalSize*internalDeg;
    }
    
    /// Nonstatic index
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    int64_t index(const int64_t& site,
		  const int& internalDeg) const
    {
      return Field::index(site,internalDeg,externalSize);
    }
    
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
    
#define FOR_EACH_SITE_DEG_OF_FIELD(F,CAPTURES,SITE,IDEG,CODE...)	\
    PAR(0,(F).nSites(),CAPTURE(CAPTURES),SITE,				\
    {									\
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
    double norm2() const
    {
      Field<Fund,FC> buf("buf");
      
      PAR(0,this->nSites(),
	  CAPTURE(t=this->getReadable(),
		  TO_WRITE(buf)),
	  site,
	  {
	    double s2=0.0;
	    UNROLL_FOR(internalDeg,0,nInternalDegs)
	      s2+=sqr(t(site,internalDeg));
	    buf[site]=s2;
	  });
      
      double res;
      glb_reduce(&res,buf,this->nSites());
      
      return res;
    }
    
    /// Put the field to the new norm, returning the reciprocal of normalizing factor
    double normalize(const double& newNorm=1.0)
    {
      const double f=
	newNorm/sqrt(this->norm2());
      
      (*this)*=f;
      
      return 1/f;
    }
    
    double normalize(const Field& oth,
		     const double& newNorm=1.0)
    {
      //compute current norm
      const double oldNorm2=
	sqrt(oth.norm2());
      
      //compute normalizing factor
      const double fact=
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
    
    void reduce(T& out) const
    { //hack
      Field tmp("tmp");
      tmp=*this;
      
      int64_t n=this->nSites();
      
      // const double init_time=take_time();
      while(n>1)
	{
	  const int64_t stride=(n+1)/2;
	  const int64_t nReductions=n/2;
	  
	  PAR(0,nReductions,
	      CAPTURE(stride,
		      TO_WRITE(tmp)),
		      ireduction,
	      {
		const int64_t first=ireduction;
		const int64_t second=first+stride;
		
		UNROLL_FOR(iDeg,0,nInternalDegs)
		  tmp(first,iDeg)+=tmp(second,iDeg);
	      });
	  n=stride;
	};
      
      UNROLL_FOR(iDeg,0,nInternalDegs)
	((Fund*)out)[iDeg]=tmp(0,iDeg);
      
      MPI_Allreduce(MPI_IN_PLACE,out,nInternalDegs,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
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
    double realPartOfScalarProdWith(const Field& oth) const
    {
      Field<double,FC> buf("buf");
      
      PAR(0,this->nSites(),
	  CAPTURE(TO_WRITE(buf),
		  t=this->getReadable(),
		  TO_READ(oth)),
	  site,
	  {
	    double r=0;
	    UNROLL_FOR(internalDeg,0,nInternalDegs)
	      r+=t(site,internalDeg)*oth(site,internalDeg);
	    buf[site]=r;
	  });
      
      double res;
      glb_reduce(&res,buf,this->nSites());
      
      return res;
    }
    
#define PROVIDE_CASTS(CONST)						\
    /* Cast to a different fieldCoverage */				\
    template <FieldCoverage NFC,					\
	      bool Force=false>						\
    CONST Field<T,NFC,FL>& castFieldCoverage() CONST			\
    {									\
      if constexpr(not (Force or FC== EVEN_OR_ODD_SITES))	\
	static_assert(NFC==EVEN_SITES or NFC==ODD_SITES, \
		      "incompatible fieldCoverage! Force the change if needed"); \
									\
      return *(CONST Field<T,NFC,FL>*)this;				\
    }									\
									\
    /* Cast to a different comps */					\
    template <typename NT,						\
	      bool Force=false>						\
    CONST Field<NT,FC,FL>& castComponents() CONST			\
    {									\
      static_assert(Force or sizeof(T)==sizeof(NT),			\
		    "incompatible components! Force the change if needed"); \
      									\
      return *(CONST Field<NT,FC,FL>*)this;				\
    }
    
    PROVIDE_CASTS(const);
    
    PROVIDE_CASTS(/* not const */);
    
#undef PROVIDE_CASTS
    
    /// Constructor
    Field(const char *name,
	  const HaloEdgesPresence& haloEdgesPresence=WITHOUT_HALO) :
      nRef(0),
      name(name),
      haloEdgesPresence(haloEdgesPresence),
      externalSize(FieldSizes<fieldCoverage>::nSitesToAllocate(haloEdgesPresence)),
      haloIsValid(false),
      edgesAreValid(false)
    {
      verbosity_lv3_master_printf("Allocating field %s\n",name);
      _data=nissa_malloc(name,externalSize*nInternalDegs,Fund);
    }
    
    /// Construct from other layout
    template <FieldLayout OFL,
	      ENABLE_THIS_TEMPLATE_IF(OFL!=FL)>
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
      externalSize(oth.externalSize),
      _data(oth._data),
      haloIsValid(oth.haloIsValid),
      edgesAreValid(oth.edgesAreValid)
    {
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
      verbosity_lv3_master_printf("Deallocating field %s\n",name);
      if(_data)
	{
	  if(nRef==0)
	    nissa_free(_data);
	  else
	    crash("Trying to destroying field %s with dangling references",name);
	}
#endif
    }
    
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    decltype(auto) operator[](const int& site) CONST			\
    {									\
      if constexpr(not std::is_array_v<T>)				\
	return								\
	  _data[site];							\
      else								\
	if constexpr(FL==FieldLayout::CPU)				\
	  return ((CONST T*)_data)[site];				\
	else								\
	  return							\
	    SubscribedField<CONST Field,				\
	    std::remove_extent_t<T>>(*this,site,nullptr);		\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* not const */);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    /// Register the actions to backup and restore the field
    constexpr INLINE_FUNCTION
    FieldRef<Field> getWritable()
    {
#ifdef USE_CUDA
      if(insideParallelFor)
	{
	  benchmarkBeginActions.emplace_back([this]()
	  {
	    verbosity_lv3_master_printf("Backing up field %s\n",name);
	    backup=new Field("backup",haloEdgesPresence);
	    *backup=*this;
	  });
	  
	  benchmarkEndActions.emplace_back([this]()
	  {
	    *this=*backup;
	    delete backup;
	    verbosity_lv3_master_printf("Restored field %s\n",name);
	  });
	  verbosity_lv3_master_printf("Added backup actions for field %s\n",name);
	}
      else
	verbosity_lv3_master_printf("Skipping backup actions for field %s as we are not inside a parallel for\n",name);
#endif
      
      return *this;
    }
    
    constexpr INLINE_FUNCTION
    FieldRef<const Field> getReadable() const
    {
      return *this;
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_FLATTENED_CALLER(CONST)					\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    CONST Fund& operator()(const int64_t& site,				\
			   const int& internalDeg) CONST		\
    {									\
      return _data[index(site,internalDeg,externalSize)];		\
    }
    
    PROVIDE_FLATTENED_CALLER(const);
    
    PROVIDE_FLATTENED_CALLER(/* not const */);
    
#undef PROVIDE_FLATTENED_CALLER
    
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
		  t=this->getReadable()),
	  i,
	  {
	    UNROLL_FOR(internalDeg,0,nInternalDegs)
	      ((Fund*)send_buf)[internalDeg+nInternalDegs*i]=
		t(f(i),internalDeg);
	  });
    }
    
    /// Fill the sending buf using the data on the surface of a field
    void fillSendingBufWithSurface() const
    {
      fillSendingBufWith([] CUDA_DEVICE(const int& i) INLINE_ATTRIBUTE
      {
	return Field::surfSiteOfHaloSite(i);
      },bord_vol/divCoeff);
    }
    
    /// Fill the sending buf using the data on the surface edge
    void fillSendingBufWithEdgesSurface() const
    {
      fillSendingBufWith([] CUDA_DEVICE(const int& i) INLINE_ATTRIBUTE
      {
	return Field::surfSiteOfEdgeSite(i);
      },edge_vol/divCoeff);
    }
    
    /// Fill the surface using the data from the buffer
    template <typename B,
	      typename F>
    void fillSurfaceWithReceivingBuf(const F& f)
    {
      for(int bf=0;bf<2;bf++)
	for(int mu=0;mu<NDIM;mu++)
	  PAR(0,bord_dir_vol[mu]/divCoeff,
	      CAPTURE(f,bf,mu,
		      t=this->getWritable()),
	      iHaloOriDir,
	      {
		const int iHalo=
		  bf*bord_volh/divCoeff+
		  iHaloOriDir+bord_offset[mu]/divCoeff;
		
		const int iSurf=
		  Field::surfSiteOfHaloSite(iHalo);
		
		f(t[iSurf],
		  ((B*)recv_buf)[iHalo],
		  bf,
		  mu);
	      });
    }
    
    /// Fill the sending buf using the data with a given function
    void fillHaloOrEdgesWithReceivingBuf(const int& offset,
					 const int& n) const
    {
      PAR(0,n,
	  CAPTURE(_data=this->_data,offset,externalSize=this->externalSize),
	  i,
	  {
	    for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	      _data[index(offset+i,internalDeg,externalSize)]=
		((Fund*)recv_buf)[internalDeg+nInternalDegs*i];
	  });
    }
    
    /// Fills the halo with the received buffer
    void fillHaloWithReceivingBuf() const
    {
      assertHasHalo();
      
      fillHaloOrEdgesWithReceivingBuf(locVol/divCoeff,bord_vol/divCoeff);
    }
    
    /// Fills the halo with the received buffer
    void fillEdgesWithReceivingBuf() const
    {
      assertHasEdges();
      
      fillHaloOrEdgesWithReceivingBuf((locVol+bord_vol)/divCoeff,edge_vol/divCoeff);
    }
    
    /// Fills the sending buffer with the halo, compressing into elements of B using f
    template <typename B,
	      typename F>
    void fillSendingBufWithHalo(const F& f) const
    {
      for(int bf=0;bf<2;bf++)
	for(int mu=0;mu<NDIM;mu++)
	  PAR(0,bord_dir_vol[mu]/divCoeff,
	      CAPTURE(bf,mu,f,
		      t=this->getReadable()),
	      iHaloOriDir,
	      {
		const int iHalo=
		  bf*bord_volh/divCoeff+
		  iHaloOriDir+bord_offset[mu]/divCoeff;
		
		f(((B*)send_buf)[iHalo],
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
	std::min(send_buf_size,recv_buf_size);
      
      if(neededBufSize>maxBufSize)
	crash("asking to create a communicator that needs %ld large buffer (%ld allocated)",
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
      verbosity_lv3_master_printf("Start communication of borders of %s\n",name);
      
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
      verbosity_lv3_master_printf("Finish communication of halo of %s\n",name);
      
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
	crash("needs halo allocated for %s!",name);
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
      verbosity_lv3_master_printf("Start communication of edges of %s\n",name);
      
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
      verbosity_lv3_master_printf("Finish communication of edges of %s\n",name);
      
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
	crash("needs edges allocated on field %s!",name);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Communicate the halo
    void updateHalo(const bool& force=false) const
    {
      if(force or not haloIsValid)
	{
	  verbosity_lv3_master_printf("Sync communication of halo of %s\n",name);
	  
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
	  verbosity_lv3_master_printf("Sync communication of edges of %s\n",name);
	  
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
      return _data==oth._data;
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
      
      FOR_EACH_SITE_DEG_OF_FIELD(*this,
				 CAPTURE(TO_READ(oth),
					 t=this->getWritable()),
				 site,
				 iDeg,
				 {
				   t(site,iDeg)=oth(site,iDeg);
				 }
      );
      
      doNotBackupDuringBenchmark=b;
    }
    
    /// Assigns from a different layout
    template <FieldLayout OFl>
    INLINE_FUNCTION
    Field& operator=(const Field<T,FC,OFl>& oth)
    {
      assign(oth);
      
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
	    FieldLayout FL=defaultFieldLayout>
  using LxField=Field<T,FULL_SPACE,FL>;
  
  /// Field over even sites
  template <typename T,
	    FieldLayout FL=defaultFieldLayout>
  using EvnField=Field<T,EVEN_SITES,FL>;
  
  /// Field over odd sites
  template <typename T,
	    FieldLayout FL=defaultFieldLayout>
  using OddField=Field<T,ODD_SITES,FL>;
  
  /// Field over even or odd sites
  template <typename T,
	    FieldLayout FL=defaultFieldLayout>
  using EvenOrOddField=Field<T,EVEN_OR_ODD_SITES,FL>;
  
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
	    FieldLayout FL=defaultFieldLayout,
	    typename Fevn=Field<T,EVEN_SITES,FL>,
	    typename Fodd=Field<T,ODD_SITES,FL>,
	    typename FevnOrOdd=Field<T,EVEN_OR_ODD_SITES,FL>>
  struct EoField
  {
    /// Type representing a pointer to type T
    template <FieldCoverage EO>
    using F=std::conditional_t<EO==EVN,Fevn,Fodd>;
    
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
    EoField<T,FL,
	    FieldRef<Field<T,EVEN_SITES,FL>>,
	    FieldRef<Field<T,ODD_SITES,FL>>,
	    FieldRef<Field<T,EVEN_OR_ODD_SITES,FL>>>
    getWritable()
    {
      return {evenPart,oddPart};
    }
    
    constexpr INLINE_FUNCTION
    EoField<T,FL,
	    FieldRef<const Field<T,EVEN_SITES,FL>>,
	    FieldRef<const Field<T,ODD_SITES,FL>>,
	    FieldRef<const Field<T,EVEN_OR_ODD_SITES,FL>>>
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
    void updateHalo() const
    {
      forBothParities([this](const auto& par)
      {
	(*this)[par].updateHalo();
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
      f(f),site(site),ptr(ptr)
    {
    }
  };
}

#endif
