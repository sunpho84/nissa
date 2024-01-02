#ifndef _FIELD_HPP
#define _FIELD_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/field.hpp

#include <base/lattice.hpp>
#include <base/universe.hpp>
#include <communicate/communicate.hpp>
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
  /// Used to dispatch the copy constructor
  struct _CopyConstructInternalDispatcher
  {
  };
  
  /// Start the communications of buffer interpreted as halo
  std::vector<MPI_Request> startBufHaloNeighExchange(const size_t& bps);
  
  /// Specifies the order of components
  template <typename TP,
	    typename F,
	    FieldLayout FL>
  struct FieldCompsProvider;
  
#define PROVIDE_FIELD_COMPS_PROVIDER(LAYOUT,SITE,TYPES...)	\
  									\
  template <typename...C,						\
	    typename F>							\
  struct FieldCompsProvider<CompsList<C...>,F,				\
			    FieldLayout::LAYOUT>			\
  {									\
    using Comps=CompsList<TYPES>;					\
    									\
    using Site=SITE;							\
    									\
    using Fund=F;							\
  }
  
  PROVIDE_FIELD_COMPS_PROVIDER(CPU,LocLxSite,LocLxSite,C...);
  
  PROVIDE_FIELD_COMPS_PROVIDER(GPU,LocLxSite,C...,LocLxSite);
  
#undef PROVIDE_FIELD_COMPS_PROVIDER
  
  /////////////////////////////////////////////////////////////////
  
#define FIELD_COMPS_PROVIDER FieldCompsProvider<CompsList<C...>,_Fund,FL>
  
#define FIELD_COMPS typename FIELD_COMPS_PROVIDER::Comps
  
#define THIS						\
  Field<CompsList<C...>,_Fund,FL,MT,IsRef>
  
#define BASE					\
  Node<THIS,FIELD_COMPS>
  
  /// Defines a field
  template <typename...C,
	    typename _Fund,
	    FieldLayout FL,
	    MemoryType MT,
	    bool IsRef>
  struct THIS :
    FieldFeat,
    SingleSubExpr<THIS>,
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
      Field<CompsList<C...>,_Fund,FieldLayout::GPU,maybeGpuMemoryType,IsRef>;
    
    template <FieldLayout OFL,
	      MemoryType OMT,
	      bool OIR>
    requires(OFL!=FL or OMT!=MT)
    INLINE_FUNCTION
    Field& operator=(const Field<CompsList<C...>,_Fund,OFL,OMT,OIR>& oth)
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
    INLINE_FUNCTION
    void assign(O&& oth)
    {
#define LOOP(LOOP_TYPE)							\
      LOOP_TYPE(0,lat->getLocVol(),					\
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
    
    /// Remap the possible Dir components from/to the nissa to/from Ildg order
    INLINE_FUNCTION constexpr
    void scidacNissaCompsDirRemap()
    {
      /// Helper holding the two set of inner components: those referring to the directions, (Valid) and the other (Invalid)
      using T=
	TupleTellApart<InnerComps,CompsList<DirRow,DirCln>>;
      
      /// Holds the possible Dir components
      using MaybeDirs=
	T::Valid;
      
      /// Holds all the other components
      using OtherComps=
	T::Invalid;
      
      if constexpr(std::tuple_size_v<MaybeDirs>)
	PAR_ON_EXEC_SPACE(execSpace,
			  0,
			  lat->getLocVol(),
			  CAPTURE(self=this->getWritable()),
			  site,
			  {
			    compsLoop<OtherComps>([&self,
						   &site](const auto&...rc)
			      {
				/// Holds the temporarily remapped component
				StackTens<MaybeDirs,Fund> tmp;
				
				compsLoop<MaybeDirs>([&tmp,
						      &self,
						      &rc...,
						      &site](const auto&...d)
				{
				  tmp(((d+1)%NDIM)...)=self(site,rc...,d...);
				},self.getDynamicSizes());
				
				self(site,rc...)=tmp;
				
			      },self.getDynamicSizes());
			    });
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
      subExpr=std::move(oth.subExpr);
      haloIsValid=oth.haloIsValid;
      haloPresence=oth.haloPresence;
      
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
    
    /// Internal storage type
    using Data=
      DynamicTens<Comps,Fund,MT,IsRef>;
    
    /// Executes where Data exec
    static constexpr ExecSpace execSpace=
      Data::execSpace;
    
    /////////////////////////////////////////////////////////////////
    
    /// Computes the size to allocate
    static Site nSitesToAllocate(const HaloPresence& haloPresence)
    {
      return lat->getLocVol()+
	((haloPresence>=WITH_HALO)?
	 2*lat->getSurfSize():0);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Total allocated sites
    Site nTotalAllocatedSites;
    
    /// Number of internal degrees of freedom - this will be made dynamic
    static constexpr int nInternalDegs=
      indexMaxValue<C...>();
    
    /// Storage data
    mutable Data subExpr;
    
#define PROVIDE_GET_DATA(ATTRIB)			\
    CUDA_HOST_AND_DEVICE constexpr INLINE_FUNCTION	\
    ATTRIB Data& getData() ATTRIB			\
    {							\
      return subExpr;					\
    }
    
    PROVIDE_GET_DATA(const);
    
    PROVIDE_GET_DATA(/* not const */);
    
#undef PROVIDE_GET_DATA
    
    /// Presence of halo
    HaloPresence haloPresence;
    
    /// States whether the halo is updated
    mutable bool haloIsValid;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr
    auto getDynamicSizes() const
    {
      return std::make_tuple(lat->getLocVol());
    }
    
#define PROVIDE_EVAL(ATTRIB)					\
    template <typename...U>					\
    CUDA_HOST_AND_DEVICE constexpr INLINE_FUNCTION		\
    ATTRIB Fund& eval(const U&...cs) ATTRIB			\
    {								\
      return subExpr(cs...);					\
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
      Field<CompsList<>,_Fund,FL,MT> buf;
      
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
      Site n=lat->getLocVol();
      
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
      Field<CompsList<C...>,Fund,FieldLayout::CPU,MT> buf(*this);
      buf.selfReduce();
      
      StackTens<CompsList<C...>,Fund> res;
      memcpy<MemoryType::CPU,MT>(res.storage,buf.subExpr.storage,res.nElements*sizeof(Fund));
      
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
    auto getRef() ATTRIB	\
    {									\
      return Field<CompsList<C...>,ATTRIB _Fund,FL,MT,true>(*this,(_CopyConstructInternalDispatcher*)nullptr); \
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
    /// Type obtained reinterpreting the fund
    template <typename NFund>
    using ReinterpretFund=
      Field<CompsList<C...>,NFund,FL,MT,IsRef>;
    
    /////////////////////////////////////////////////////////////////
    
    using FlattenedInnerComp=MergedComp<InnerComps>;
    
#define PROVIDE_FLATTEN(ATTRIB)						\
    ATTRIB auto& flatten() ATTRIB					\
    {									\
      return *((Field<CompsList<FlattenedInnerComp>,ATTRIB _Fund,FL,MT,true>*)this); \
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
      Field<CompsList<C...>,std::decay_t<_Fund>,FL,MT>;
    
    /// Parameters to recreate an equivalent storage
    auto getEquivalentStoragePars() const
    {
      return std::make_tuple(haloPresence);
    }
    
    /// Creates a copy
    ClosingType createEquivalentStorage() const
    {
      return haloPresence;
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
    
    void allocate(const HaloPresence& _haloPresence=WITHOUT_HALO)
    {
      nTotalAllocatedSites=nSitesToAllocate(haloPresence);
      subExpr.allocate(std::make_tuple(nTotalAllocatedSites));
      haloPresence=_haloPresence;
      
      invalidateHalo();
    }
    
    /// Create a field
    Field(const HaloPresence& haloPresence=WITHOUT_HALO) :
      haloPresence(haloPresence)
    {
      static_assert(not IsRef,"Can allocate only if not a reference");
      
      allocate(haloPresence);
    }
    
    /// Assign another expression
    template <DerivedFromNode O>
    INLINE_FUNCTION
    Field(const O& oth) :
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
    
    /// Copy constructor, internal implementation
    template <typename O>
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    Field(O&& oth,
	   _CopyConstructInternalDispatcher*) :
      nTotalAllocatedSites(oth.nTotalAllocatedSites),
      subExpr(oth.subExpr),
      haloPresence(oth.haloPresence),
      haloIsValid(oth.haloIsValid)
    {
#ifndef COMPILING_FOR_DEVICE
      if constexpr(not IsRef)
	verbosity_lv3_master_printf("Using copy constructor of Field, isRef: %d\n",IsRef);
#endif
    }
    
    /// Return a copy on the given memory space
    template <MemoryType OMT>
    Field<CompsList<C...>,_Fund,FL,OMT> copyToMemorySpace() const
    {
      Field<CompsList<C...>,_Fund,FL,OMT> res(haloPresence);
      res.subExpr=subExpr;
      
      return res;
    }
    
    /* Return a copy on the given memory space, only if needed */
#define PROVIDE_COPY_TO_MEMORY_SPACE_IF_NEEDED(ATTRIB,REF,RETURNED_IN_SAME_MT_CASE)	\
    template <MemoryType OMT>						\
    std::conditional_t<OMT==MT,						\
		       RETURNED_IN_SAME_MT_CASE,			\
		       Field<CompsList<C...>,Fund,FL,OMT>>		\
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
    explicit Field(O&& oth) :
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
      subExpr(std::move(oth.subExpr)),
      haloPresence(oth.haloPresence),
      haloIsValid(oth.haloIsValid)
    {
#ifndef COMPILING_FOR_DEVICE
      verbosity_lv3_master_printf("Using move constructor of Field\n");
#endif
    }
    
    /// Construct from another exec space and/or field layout and/or memory layout
    template <FieldLayout OFL,
	      MemoryType OMT,
	      bool OIR>
    requires(OFL!=FL or OMT!=MT)
    INLINE_FUNCTION
    Field(const Field<CompsList<C...>,_Fund,OFL,OMT,OIR>& oth) :
      Field(oth.haloPresence)
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
    Field(const Field<CompsList<C...>,_Fund,OFL,MT,OIR>& oth) :
      Field(oth.haloPresence)
    {
      (*this)=oth;
    }
    
    /// Set halo as invalid
    void invalidateHalo()
    {
      haloIsValid=false;
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
      fillSendingBufWith([lat=lat->getRef()] CUDA_DEVICE(const Site& i) INLINE_ATTRIBUTE
      {
	return lat.getSurfSiteOfHaloSite(i);
      },2*lat->getSurfSize());
      
      // for(size_t i=0;i<bord_vol*StackTens<CompsList<C...>,Fund>::nElements;i++)
      // 	master_printf("s %zu %lg\n",i,((Fund*)send_buf)[i]);
    }
    
    /// Fill the surface using the data from the buffer
    template <typename B,
	      typename F>
    void fillSurfaceWithReceivingBuf(const F& f)
    {
      for(int bf=0;bf<2;bf++)
	for(Dir mu=0;mu<NDIM;mu++)
	  PAR(0,lat->getSurfSizePerDir(mu),
	      CAPTURE(f,bf,mu,
		      n=lat->getSurfSize(),
		      off=lat->getSurfOffsetOfDir(mu),
		      t=this->getWritable()),
	      iHaloOriDir,
	      {
		const LocLxSite iHalo=
		  bf*n+
		  iHaloOriDir+off;
		
		const LocLxSite iSurf=
		  lat->getSurfSiteOfHaloSite(iHalo);
		
		f(t[iSurf],
		  ((B*)recv_buf)[iHalo],
		  bf,
		  mu);
	      });
    }
    
    /// Fill the sending buf using the data with a given function
    INLINE_FUNCTION
    void fillHaloWithReceivingBuf() const
    {
      PAR(0,2*lat->getSurfSize(),
	  CAPTURE(subExpr=this->subExpr.getWritable(),
		  offset=lat->getLocVol(),
		  dynamicSizes=this->getDynamicSizes()),
	  i,
	  {
	    compsLoop<InnerComps>([offset,
				   i,
				   &subExpr,
				   dynamicSizes] CUDA_DEVICE(const auto&...c) INLINE_ATTRIBUTE
	    {
	      const auto internalDeg=index(dynamicSizes,c...);
	      
	      asMutable(subExpr(offset+i,c...))=
		((Fund*)recv_buf)[internalDeg+nInternalDegs*i()];
	    },dynamicSizes);
	  });
      
      // for(size_t i=0;i<bord_vol*StackTens<CompsList<C...>,Fund>::nElements;i++)
      // 	master_printf("r %zu %lg\n",i,((Fund*)recv_buf)[i]);
    }
    
    /// Fills the sending buffer with the halo, compressing into elements of B using f
    template <typename B,
	      typename F>
    void fillSendingBufWithHalo(const F& f) const
    {
      for(int bf=0;bf<2;bf++)
	for(Dir mu=0;mu<NDIM;mu++)
	  PAR(0,lat->getSurfSizePerDir()[mu],
	      CAPTURE(bf,mu,f,
		      n=lat->getSurfSize(),
		      l=lat->getLocVol(),
		      off=lat->getSurfOffsetOfDir(mu),
		      t=this->getReadable()),
	      iHaloOriDir,
	      {
		const int iHalo=
		  bf*n+
		  iHaloOriDir+off;
		
		f(((B*)send_buf)[iHalo],
		  t[l+iHalo],
		  bf,
		  mu);
	      });
    }
    
    /////////////////////////////////////////////////////////////////
    
    static std::vector<MPI_Request> startBufHaloNeighExchange(const size_t& bps)
    {
    /// Pending requests
      std::vector<MPI_Request> requests;
      
      for(Dir mu=0;mu<NDIM;mu++)
	if(isDirParallel[mu])
	  for(Ori sendOri=0;sendOri<2;sendOri++)
	    {
	      const auto sendOrRecv=
		[&mu,
		 &bps,
		 ns=lat->getSurfSize()(),
		 &requests,
		 off=lat->getSurfOffsetOfDir(mu)(),
		 messageTag=sendOri()+2*mu()]
		(const char* oper,
		 const auto sendOrRecv,
		 auto* ptr,
		 const Ori& ori)
		{
		  const size_t offset=(off+ns*ori)*bps;
		  const MpiRank neighRank=neighRanks(ori,mu);
		  const size_t messageLength=2*lat->getSurfSizePerDir()[mu]()*bps;
		  // printf("rank %d %s ori %d dir %d, corresponding rank: %d, tag: %d length: %zu\n"
		  //        ,rank,oper,ori,mu,neighRank,messageTag,messageLength);
		  
		  MPI_Request request;
		  
		  sendOrRecv(ptr+offset,messageLength,MPI_CHAR,neighRank(),
			     messageTag,MPI_COMM_WORLD,&request);
		  
		  requests.push_back(request);
		};
	      
	      const Ori recvOri=1-sendOri;
	      
	      sendOrRecv("send",MPI_Isend,send_buf,sendOri);
	      sendOrRecv("recv",MPI_Irecv,recv_buf,recvOri);
	    }
      
      return requests;
    }
    
    /// Start the communications of halo
    static std::vector<MPI_Request> startHaloAsyincComm()
    {
      return startBufHaloNeighExchange(nInternalDegs*sizeof(Fund));
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
      
      assertCanCommunicate(2*lat->getSurfSize());
      
      //fill the communicator buffer, start the communication and take time
      fillSendingBufWithSurface();
      requests=startHaloAsyincComm();
      
      return requests;
    }
    
    /// Finalize communications
    void finishCommunicatingHalo(std::vector<MPI_Request> requests) const
    {
      //take note of passed time and write some debug info
      verbosity_lv3_master_printf("Finish communication of halo\n");
      
      //wait communication to finish, fill back the vector and take time
      waitAsyncCommsFinish(requests);
      fillHaloWithReceivingBuf();
      
      haloIsValid=true;
    }
    
    /// Crash if the halo is not allocated
    void assertHasHalo() const
    {
      if(not (haloPresence>=WITH_HALO))
	crash("needs halo allocated!");
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
  };
  
  template <typename T,
	    DerivedFromComp...Ci>
  INLINE_FUNCTION constexpr
  auto Node<T,CompsList<Ci...>>::closeToField() const
    requires(_canCloseToField())
  {
    return (Field<TupleFilterAllTypes<typename T::Comps,CompsList<LocLxSite>>,
	    typename T::Fund>)*this;
  }
}

#endif
