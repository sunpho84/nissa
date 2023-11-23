#ifndef _FIELD_HPP
#define _FIELD_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/field.hpp

#include <base/old_field.hpp>
#include <expr/assignDispatcher.hpp>
#include <expr/comps.hpp>
#include <expr/dynamicTens.hpp>
#include <expr/fieldDeclaration.hpp>
#include <expr/mergedComps.hpp>
#include <expr/stackTens.hpp>
#include <geometry/geometry_eo.hpp>
#include <threads/threads.hpp>

/// \todo: Components must all be provided externally

namespace nissa
{
  /// Start the communications of buffer interpreted as halo
  std::vector<MPI_Request> startBufHaloNeighExchange(const int& divCoeff,
						     const size_t& bps);
  
  DECLARE_UNTRANSPOSABLE_COMP(Parity,int,2,createParity);
  DECLARE_UNTRANSPOSABLE_COMP(Dir,int,NDIM,dir);
  
  DECLARE_UNTRANSPOSABLE_COMP(Ori,int,2,createOri);
  
  /// Backward, see real imag comment
#define bw Ori(0)
  
  /// Forward
#define fw Ori(1)
  
  /// Number of dimensions
#define nDim Dir(NDIM)
  
  DECLARE_PARALLELIZABLE_COMP(LocLxSite,int64_t,locLxSite);
  DECLARE_PARALLELIZABLE_COMP(LocEoSite,int64_t,locEoSite);
  DECLARE_PARALLELIZABLE_COMP(LocEvnSite,int64_t,locEvnSite);
  DECLARE_PARALLELIZABLE_COMP(LocOddSite,int64_t,locOddSite);

  DECLARE_PARALLELIZABLE_COMP(GlbLxSite,int64_t,glbLxSite);
  
  /// Specifies the order of components
  template <typename TP,
	    typename F,
	    FieldCoverage FC,
	    FieldLayout FL>
  struct FieldCompsProvider;
  
#define PROVIDE_FIELD_COMPS_PROVIDER(COVERAGE,LAYOUT,SITE,TYPES...)	\
  									\
  template <typename...C,						\
	    typename F>							\
  struct FieldCompsProvider<CompsList<C...>,F,				\
			    FieldCoverage::COVERAGE,			\
			    FieldLayout::LAYOUT>			\
  {									\
    using Comps=CompsList<TYPES>;					\
    									\
    using Site=SITE;							\
    									\
    using Fund=F;							\
  }
  
  PROVIDE_FIELD_COMPS_PROVIDER(FULL_SPACE,CPU,LocLxSite,LocLxSite,C...);
  PROVIDE_FIELD_COMPS_PROVIDER(EVEN_OR_ODD_SITES,CPU,LocEoSite,Parity,LocEoSite,C...);
  PROVIDE_FIELD_COMPS_PROVIDER(EVEN_SITES,CPU,LocEvnSite,LocEvnSite,C...);
  PROVIDE_FIELD_COMPS_PROVIDER(ODD_SITES,CPU,LocOddSite,LocOddSite,C...);
  
  PROVIDE_FIELD_COMPS_PROVIDER(FULL_SPACE,GPU,LocLxSite,C...,LocLxSite);
  PROVIDE_FIELD_COMPS_PROVIDER(EVEN_OR_ODD_SITES,GPU,LocEoSite,C...,Parity,LocEoSite);
  PROVIDE_FIELD_COMPS_PROVIDER(EVEN_SITES,GPU,LocEvnSite,C...,LocEvnSite);
  PROVIDE_FIELD_COMPS_PROVIDER(ODD_SITES,GPU,LocOddSite,C...,LocOddSite);
  
#undef PROVIDE_FIELD_COMPS_PROVIDER
  
  /////////////////////////////////////////////////////////////////
  
#define FIELD_COMPS_PROVIDER FieldCompsProvider<CompsList<C...>,_Fund,FC,FL>
  
#define FIELD_COMPS typename FIELD_COMPS_PROVIDER::Comps
  
#define THIS						\
  Field<CompsList<C...>,_Fund,FC,FL,MT,IsRef>
  
#define BASE					\
  Node<THIS,FIELD_COMPS>
  
  /// Defines a field
  template <typename...C,
	    typename _Fund,
	    FieldCoverage FC,
	    FieldLayout FL,
	    MemoryType MT,
	    bool IsRef>
  struct THIS :
    FieldFeat,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    /// Importing assignment operator from BaseTens
    using Base::operator=;
    
    ///  Equivalent Field on device
    using DeviceEquivalent=
      Field<CompsList<C...>,_Fund,FC,FieldLayout::GPU,maybeGpuMemoryType,IsRef>;
    
    template <FieldLayout OFL,
	      MemoryType OMT,
	      bool OIR>
    requires(OFL!=FL or OMT!=MT)
    INLINE_FUNCTION
    Field& operator=(const Field<CompsList<C...>,_Fund,FC,OFL,OMT,OIR>& oth)
    {
      if(this->getDynamicSizes()!=oth.getDynamicSizes())
	crash("trying to assign fields on different memory space, having different dynamic sizes");
      
      this->assign(oth.template copyToMemorySpaceIfNeeded<MT>());
      
      invalidateHalo();
      
      return *this;
    }
    
    /// Assign from another expression
    template <typename OP=DirectAssign,
	      typename O>
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    void assign(O&& oth)
    {
      if constexpr(compilingForDevice)
	assert(compilingForDevice);
      else
	{
#define LOOP(LOOP_TYPE)							\
      LOOP_TYPE(0,this->nSites(),					\
		CAPTURE(self=this->getWritable(),			\
			TO_READ(oth)),					\
		site,							\
		{							\
		  using RhsComps=typename std::decay_t<O>::Comps;	\
		  							\
		  /*! we need to take care that the rhs might not have the site (such in the case of a scalar) */ \
		  							\
		  if constexpr(tupleHasType<RhsComps,Site>)		\
		    OP::dispatch(self(site),oth(site));			\
		  else							\
		    OP::dispatch(self(site),oth);			\
		})
      
#ifdef USE_CUDA
      if constexpr(MT==MemoryType::GPU)
	LOOP(DEVICE_PARALLEL_LOOP);
      else
#endif
	if constexpr(MT==MemoryType::CPU)
	  LOOP(HOST_PARALLEL_LOOP);
	else
	  crash("unkwnown condition");
      
#undef LOOP
	}
    }
    
    /// Copy assign
    INLINE_FUNCTION
    Field& operator=(const Field& oth)
    {
      assign(oth);
      
      return *this;
    }
    
    /// Move assign
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    Field& operator=(Field&& oth)
    {
      nTotalAllocatedSites=oth.nTotalAllocatedSites;
      data=std::move(oth.data);
      haloIsValid=oth.haloIsValid;
      edgesAreValid=oth.edgesAreValid;
      haloEdgesPresence=oth.haloEdgesPresence;
      
      return *this;
    }
    
    /// Assigns another node
    template <DerivedFromNode O>
    INLINE_FUNCTION
    Field& operator=(const O& oth)
    {
      assign(oth);
      
      return *this;
    }
    
    /// Sumassigns another node
    template <DerivedFromNode O>
    INLINE_FUNCTION
    Field& operator+=(const O& oth)
    {
      assign<SumAssign>(oth);
      
      return *this;
    }
    
    // /// Assigns a scalar
    // Field& operator=(const _Fund& oth)
    // {
    //   *this=Scalar<_Fund>(oth);
      
    //   return *this;
    // }
    
    static constexpr FieldCoverage fieldCoverage=FC;
    
    static constexpr FieldLayout fieldLayout=FL;
    
    /// Internal components
    using InnerComps=
      CompsList<C...>;
    
    /// Components
    using Comps=
      FIELD_COMPS;
    
    /// Import dynamic comps
    using DynamicComps=
      typename DynamicCompsProvider<FIELD_COMPS>::DynamicComps;
    
    /// Determine if we have dynamic comps
    using DynamicCompsProvider<Comps>::hasDynamicComps;
    
    /// Type used for the site
    using Site=
      typename FIELD_COMPS_PROVIDER::Site;
    
    /// Fundamental tye
    using Fund=
      typename FIELD_COMPS_PROVIDER::Fund;
    
#undef FIELD_COMPS
    
#undef FIELD_COMPS_PROVIDER
    
    /// Coefficient which divides the space time, if the field is covering only half the space
    static constexpr const int divCoeff=
      (FC==FULL_SPACE)?1:2;
    
    /// Internal storage type
    using Data=
      DynamicTens<Comps,Fund,MT,IsRef>;
    
    /// Executes where Data exec
    static constexpr ExecSpace execSpace=
      Data::execSpace;
    
    /////////////////////////////////////////////////////////////////
    
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
    static constexpr const Site nSites()
    {
      if constexpr(fieldCoverage==FULL_SPACE)
	return locVol;
      else
	return locVolh;
    }
    
    /// Number of sites in the halo of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr const Site nHaloSites()
    {
      if constexpr(fieldCoverage==FULL_SPACE)
	return bord_vol;
      else
	return bord_volh;
    }
    
    /// Number of sites in the edges of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr const Site nEdgesSites()
    {
      if constexpr(fieldCoverage==FULL_SPACE)
	return edge_vol;
      else
	return edge_volh;
    }
    
    /// Computes the size to allocate
    static Site nSitesToAllocate(const HaloEdgesPresence& haloEdgesPresence)
    {
      Site res=nSites();
      
      if(haloEdgesPresence>=WITH_HALO)
	res+=nHaloSites();
      
      if(haloEdgesPresence>=WITH_HALO_EDGES)
	res+=nEdgesSites();
      
      return res;
    }
    
    /// Surface site of a site in the halo
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static Site surfSiteOfHaloSite(const Site& iHalo)
    {
      assertHasDefinedCoverage();
      
      if constexpr(fieldCoverage==FULL_SPACE)
	return surflxOfBordlx[iHalo()];
      else
	if constexpr(fieldCoverage==EVEN_SITES or fieldCoverage==ODD_SITES)
	  return surfeo_of_bordeo[fieldCoverage][iHalo()];
    }
    
    /// Surface site of a site in the e
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static Site surfSiteOfEdgeSite(const Site& iEdge)
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
    static Site locNeigh ## UD(const Site& site,			\
			       const Dir& mu)				\
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
    static auto glbCoord(const Site& site,
			 const Dir& mu)
    {
      assertHasDefinedCoverage();
      
      if constexpr(fieldCoverage==FULL_SPACE)
	return glbCoordOfLoclx[site][mu];
      else
	if constexpr(fieldCoverage==EVEN_SITES or fieldCoverage==ODD_SITES)
	  return glbCoordOfLoclx[loclx_of_loceo[fieldCoverage][site]][mu];
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Total allocated sites
    Site nTotalAllocatedSites;
    
    /// Number of internal degrees of freedom - this will be made dynamic
    static constexpr int nInternalDegs=
      indexMaxValue<C...>();
    
    /// Storage data
    mutable Data data;
    
    /// Presence of halo and edges
    HaloEdgesPresence haloEdgesPresence;
    
    /// States whether the halo is updated
    mutable bool haloIsValid;
    
    /// States whether the edges are updated
    mutable bool edgesAreValid;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    auto getDynamicSizes() const
    {
      return std::make_tuple(nSites());
    }
    
#define PROVIDE_EVAL(ATTRIB)					\
    template <typename...U>					\
    CUDA_HOST_AND_DEVICE constexpr INLINE_FUNCTION		\
    ATTRIB Fund& eval(const U&...cs) ATTRIB			\
    {								\
      return data(cs...);					\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/* non const */);
    
#undef PROVIDE_EVAL
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=
      not std::is_const_v<Fund>;
    
    /// Describe a field
    void describe(const std::string& pref="") const
    {
      master_printf("%sField %s address %p\n",pref.c_str(),demangle(typeid(*this).name()).c_str(),this);
      master_printf("%s Components:\n",pref.c_str());
      (master_printf("%s  %s\n",pref.c_str(),demangle(typeid(C).name()).c_str()),...);
      master_printf("%s Site type: %s\n",pref.c_str(),demangle(typeid(decltype(this->nSites())).name()).c_str());
      master_printf("%s Fund: %s\n",pref.c_str(),demangle(typeid(_Fund).name()).c_str());
      master_printf("%s FieldCoverage: %d\n",pref.c_str(),FC);
      master_printf("%s FieldLayout: %d\n",pref.c_str(),FL);
      master_printf("%s MemoryType: %d\n",pref.c_str(),MT);
      master_printf("%s FieldLayout: %d\n",pref.c_str(),FL);
      master_printf("%s IsRef: %d\n",pref.c_str(),IsRef);
      master_printf("%sEnd of Field\n",pref.c_str());
    }
    /////////////////////////////////////////////////////////////////
    
    /// Computes the squared norm, overloading default expression
    Fund norm2() const
    {
      Field<CompsList<>,_Fund,FC,FL,MT> buf;
      
      PAR(0,this->nSites(),
	  CAPTURE(t=this->getReadable(),
		  TO_WRITE(buf)),
	  site,
	  {
	    buf(site)=t(site).norm2();
	  });
      
      return buf.glbReduce();
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Make the origin have the sum over all sites
    void selfReduce()
    {
      Site n=nSites();
      
      while(n>1)
	{
	  const Site stride=(n+1)/2;
	  const Site nReductions=n/2;
	  //verbosity_lv3_master_printf("n: %d, stride: %d, nreductions: %d \n",n(),stride(),nReductions());
	  
	  PAR(0,nReductions,
	      CAPTURE(stride,
		      t=this->getWritable()),
	      iReduction,
	      {
		const Site first=iReduction;
		const Site second=first+stride;
		
		t(first)=t(first)+t(second);
	      });
	  
	  n=stride;
	}
    }
    
    /// Returns the sum over all sites
    StackTens<CompsList<C...>,Fund> locReduce() const
    {
      //verbosity_lv2_master_printf("n: %d, nori: %d\n",n(),nOri());
      
      /// Make spacetime the external component
      Field<CompsList<C...>,Fund,FieldCoverage::FULL_SPACE,FieldLayout::CPU,MT> buf(*this);
      buf.selfReduce();
      
      StackTens<CompsList<C...>,Fund> res;
      memcpy<MemoryType::CPU,MT>(res.storage,buf.data.storage,res.nElements*sizeof(Fund));
      
      return res;
    }
    
    /// Performs a global reduction
    auto glbReduce() const
    {
      auto res=locReduce();
      
      // StackTens<OfComps<C...>,Fund> res=;
      
      MPI_Allreduce(MPI_IN_PLACE,res.storage,res.nElements,
		    MPI_Datatype_of<Fund>(),MPI_SUM,MPI_COMM_WORLD);
      
      return res;
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)						\
    Field<CompsList<C...>,ATTRIB _Fund,FC,FL,MT,true> getRef() ATTRIB	\
    {									\
      return *this;							\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
    using FlattenedInnerComp=MergedComp<InnerComps>;
    
#define PROVIDE_FLATTEN(ATTRIB)						\
    ATTRIB auto& flatten() ATTRIB					\
    {									\
      return *((Field<CompsList<FlattenedInnerComp>,ATTRIB _Fund,FC,FL,MT,true>*)this); \
    }
    
    PROVIDE_FLATTEN(const);
    
    PROVIDE_FLATTEN(/* non const */);
    
#undef PROVIDE_FLATTEN
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_RECREATE_FROM_EXPR(ATTRIB)			\
    /*! Returns itself */					\
    INLINE_FUNCTION						\
    decltype(auto) recreateFromExprs() ATTRIB			\
    {								\
      return *this;						\
    }
    
    PROVIDE_RECREATE_FROM_EXPR(/* non const */);
    
    PROVIDE_RECREATE_FROM_EXPR(const);
    
#undef PROVIDE_RECREATE_FROM_EXPR
    
    /// Type to define the closing of the expression
    using ClosingType=
      Field<CompsList<C...>,std::decay_t<_Fund>,FC,FL,MT>;
    
    /// Parameters to recreate an equivalent storage
    auto getEquivalentStoragePars() const
    {
      return std::make_tuple(haloEdgesPresence);
    }
    
    /// Creates a copy
    ClosingType createEquivalentStorage() const
    {
      return haloEdgesPresence;
    }
    
    /// Creates a copy
    auto createCopy() const
    {
      return createEquivalentStorage()=*this;
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// States whether the field can be simdified
    static constexpr bool canSimdify=
      false;//      FL!=FieldLayout::SIMDIFIED and Data::canSimdify;
    
    /// Simdifying component
    using SimdifyingComp=
      void;//      typename Data::SimdifyingComp;
    
    /// We keep referring to the original object
    static constexpr bool storeByRef=not IsRef;
    
    /// Returns that can assign
    constexpr INLINE_FUNCTION
    bool canAssign()
    {
      return canAssignAtCompileTime;
    }
    
    /// Constructor avoiding allocation
    Field(const DoNotAllocate)
    requires(not IsRef)
    {
      master_printf("avoiding allocation\n");
    }
    
    void allocate(const HaloEdgesPresence& _haloEdgesPresence=WITHOUT_HALO)
    {
      nTotalAllocatedSites=FieldSizes<fieldCoverage>::nSitesToAllocate(haloEdgesPresence);
      data.allocate(std::make_tuple(nTotalAllocatedSites));
      haloEdgesPresence=_haloEdgesPresence;
      
      invalidateHalo();
    }
    
    /// Create a field
    Field(const HaloEdgesPresence& haloEdgesPresence=WITHOUT_HALO) :
      haloEdgesPresence(haloEdgesPresence)
    {
      static_assert(not IsRef,"Can allocate only if not a reference");
      
      allocate(haloEdgesPresence);
    }
    
    /// Assign another expression
    template <DerivedFromNode O>
    INLINE_FUNCTION
    explicit Field(const O& oth) :
      Field()
    {
      assign<DirectAssign>(oth);
    }
    
    /// Assign from fund
    INLINE_FUNCTION
    explicit Field(const Fund& oth) :
      Field()
    {
      assign<DirectAssign>(scalar(oth));
    }
    
    /// Used to dispatch the copy constructor
    struct _CopyConstructInternalDispatcher;
    
    /// Copy constructor, internal implementation
    template <typename O>
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    Field(O&& oth,
	   _CopyConstructInternalDispatcher*) :
      nTotalAllocatedSites(oth.nTotalAllocatedSites),
      data(oth.data),
      haloEdgesPresence(oth.haloEdgesPresence),
      haloIsValid(oth.haloIsValid),
      edgesAreValid(oth.edgesAreValid)
    {
#ifndef COMPILING_FOR_DEVICE
      if constexpr(not IsRef)
	verbosity_lv3_master_printf("Using copy constructor of Field, isRef: %d\n",IsRef);
#endif
    }
    
    /// Return a copy on the given memory space
    template <MemoryType OMT>
    Field<CompsList<C...>,_Fund,FC,FL,OMT> copyToMemorySpace() const
    {
      Field<CompsList<C...>,_Fund,FC,FL,OMT> res(haloEdgesPresence);
      res.data=data;
      
      return res;
    }
    
    /* Return a copy on the given memory space, only if needed */
#define PROVIDE_COPY_TO_MEMORY_SPACE_IF_NEEDED(ATTRIB,REF,RETURNED_IN_SAME_MT_CASE)	\
    template <MemoryType OMT>						\
    std::conditional_t<OMT==MT,						\
		       RETURNED_IN_SAME_MT_CASE,			\
		       Field<CompsList<C...>,Fund,FC,FL,OMT>>		\
    copyToMemorySpaceIfNeeded() ATTRIB REF				\
    {									\
      return *this;							\
    }
    
    PROVIDE_COPY_TO_MEMORY_SPACE_IF_NEEDED(const,&,const Field&);
    
    PROVIDE_COPY_TO_MEMORY_SPACE_IF_NEEDED(/* non const */,&,Field&);
    
    PROVIDE_COPY_TO_MEMORY_SPACE_IF_NEEDED(/* non const */,&&,Field);
    
#undef PROVIDE_COPY_TO_MEMORY_SPACE_IF_NEEDED
    
    /// Copy construct, taking as input a non-reference when this is a reference
    template <DerivedFromField O>
    requires(IsRef)
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    Field(O&& oth) :
      Field(std::forward<O>(oth),(_CopyConstructInternalDispatcher*)nullptr)
    {
    }
    
    /// Copy constructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    Field(const Field& oth) :
      Field(oth,(_CopyConstructInternalDispatcher*)nullptr)
    {
    }
    
    /// Move constructor
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    Field(Field&& oth) :
      nTotalAllocatedSites(oth.nTotalAllocatedSites),
      data(oth.data),
      haloEdgesPresence(oth.haloEdgesPresence),
      haloIsValid(oth.haloIsValid),
      edgesAreValid(oth.edgesAreValid)
    {
#ifndef COMPILING_FOR_DEVICE
      verbosity_lv3_master_printf("Using move constructor of Field\n");
#endif
    }
    
    /// Construct from another exec space and/or field layout
    template <FieldLayout OFL,
	      MemoryType OMT,
	      bool OIR>
    requires(OFL!=FL or OMT!=MT)
    INLINE_FUNCTION
    Field(const Field<CompsList<C...>,_Fund,FC,OFL,OMT,OIR>& oth) :
      Field(oth.haloEdgesPresence)
    {
      if constexpr(OMT!=MT)
	(*this)=oth.template copyToMemorySpace<MT>();
      else
	(*this)=oth;
    }
    
    /// Construct from another exec space and/or field layout
    template <FieldLayout OFL,
	      bool OIR>
    INLINE_FUNCTION
    Field(const Field<CompsList<C...>,_Fund,FC,OFL,MT,OIR>& oth) :
      Field(oth.haloEdgesPresence)
    {
      (*this)=oth;
    }
    
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
    
    /////////////////////////////////////////////////////////////////
    
    /// Fill the sending buf using the data with a given function
    template <typename F>
    void fillSendingBufWith(const F& f,
			    const Site& n) const
    {
      PAR(0,n,
	  CAPTURE(f,
		  t=this->getReadable(),
		  dynamicSizes=this->getDynamicSizes()),
	  i,
	  {
	    compsLoop<InnerComps>([i,
				   t,
				   dynamicSizes,
				   f] CUDA_DEVICE(const auto&...c) INLINE_ATTRIBUTE
	    {
	      const auto internalDeg=index(dynamicSizes,c...);
	      
	      ((std::remove_const_t<Fund>*)send_buf)[internalDeg+nInternalDegs*i()]=
		t(f(i),c...);
	    },dynamicSizes);
	  });
    }
    
    /// Fill the sending buf using the data on the surface of a field
    void fillSendingBufWithSurface() const
    {
      fillSendingBufWith([] CUDA_DEVICE(const Site& i) INLINE_ATTRIBUTE
      {
	return surfSiteOfHaloSite(i);
      },bord_vol/divCoeff);
      
      // for(size_t i=0;i<bord_vol*StackTens<CompsList<C...>,Fund>::nElements;i++)
      // 	master_printf("s %zu %lg\n",i,((Fund*)send_buf)[i]);
    }
    
    /// Fill the sending buf using the data on the surface edge
    void fillSendingBufWithEdgesSurface() const
    {
      fillSendingBufWith([] CUDA_DEVICE(const Site& i) INLINE_ATTRIBUTE
      {
	return surfSiteOfEdgeSite(i);
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
		  surfSiteOfHaloSite(iHalo);
		
		f(t[iSurf],
		  ((B*)recv_buf)[iHalo],
		  bf,
		  mu);
	      });
    }
    
    /// Fill the sending buf using the data with a given function
    INLINE_FUNCTION
    void fillHaloOrEdgesWithReceivingBuf(const Site& offset,
					 const Site& n) const
    {
      PAR(0,n,
	  CAPTURE(data=this->data.getWritable(),
		  offset,
		  dynamicSizes=this->getDynamicSizes()),
	  i,
	  {
	    compsLoop<InnerComps>([offset,
				   i,
				   &data,
				   dynamicSizes] CUDA_DEVICE(const auto&...c) INLINE_ATTRIBUTE
	    {
	      const auto internalDeg=index(dynamicSizes,c...);
	      
	      asMutable(data(offset+i,c...))=
		((Fund*)recv_buf)[internalDeg+nInternalDegs*i()];
	    },dynamicSizes);
	  });
      
      // for(size_t i=0;i<bord_vol*StackTens<CompsList<C...>,Fund>::nElements;i++)
      // 	master_printf("r %zu %lg\n",i,((Fund*)recv_buf)[i]);
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
      return startBufHaloNeighExchange(divCoeff,nInternalDegs*sizeof(Fund));
    }
    
    /// Start the communications of edges
    static std::vector<MPI_Request> startEdgesAsyincComm()
    {
      return startBufEdgesNeighExchange(divCoeff,nInternalDegs*sizeof(Fund));
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Assert that can communicate n sites
    static void assertCanCommunicate(const Site& n)
    {
      /// Needed size in the buffer
      const size_t neededBufSize=
	nInternalDegs*n();
      
      const size_t maxBufSize=
	std::min(send_buf_size,recv_buf_size);
      
      if(neededBufSize>maxBufSize)
	crash("asking to create a communicator that needs %d large buffer (%d allocated)",
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
      verbosity_lv3_master_printf("Start communication of borders\n");
      
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
      verbosity_lv3_master_printf("Finish communication of halo\n");
      
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
	crash("needs halo allocated!");
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
      verbosity_lv3_master_printf("Start communication of edges\n");
      
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
      verbosity_lv3_master_printf("Finish communication of edgess\n");
      
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
	crash("needs edges allocated!");
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Communicate the halo
    void updateHalo(const bool& force=false) const
    {
      if(force or not haloIsValid)
	{
	  verbosity_lv3_master_printf("Sync communication of halo\n");
	  
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
	  verbosity_lv3_master_printf("Sync communication of edges\n");
	  
	  const std::vector<MPI_Request> requests=
	    startCommunicatingEdges();
	  finishCommunicatingEdges(requests);
      }
    }
  };
}

#endif
