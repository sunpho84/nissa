#ifndef _FIELD2_HPP
#define _FIELD2_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file grill/field.hpp

#include <expr/comps.hpp>
#include <expr/dynamicTens.hpp>
#include <expr/fieldDeclaration.hpp>

// #include <lattice/parityProvider.hpp>
// #include <lattice/lattice.hpp>
// #include <resources/Mpi.hpp>

namespace nissa
{
  /// Start the communications of buffer interpreted as halo
  std::vector<MPI_Request> startBufHaloNeighExchange(const int& divCoeff,
						     const size_t& bps);
  
  DECLARE_UNTRANSPOSABLE_COMP(Parity,int,2,createParity);
  DECLARE_UNTRANSPOSABLE_COMP(Dir,int,NDIM,createDir);
  
  DECLARE_UNTRANSPOSABLE_COMP(Ori,int,2,createOri);
  
  /// Backward, see real imag comment
#define bw Ori(0)
  
  /// Forward
#define fw Ori(1)
  
  DECLARE_PARALLELIZABLE_COMP(LocLxSite,int,locLxSite);
  DECLARE_PARALLELIZABLE_COMP(LocEoSite,int,locEoSite);
  DECLARE_PARALLELIZABLE_COMP(LocEvnSite,int,locEvnSite);
  DECLARE_PARALLELIZABLE_COMP(LocOddSite,int,locOddSite);
  
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
  Field2<CompsList<C...>,_Fund,FC,FL,MT,IsRef>
  
#define BASE					\
  Node<THIS>
  
  /// Defines a field
  template <typename...C,
	    typename _Fund,
	    FieldCoverage FC,
	    FieldLayout FL,
	    MemoryType MT,
	    bool IsRef>
  struct THIS :
    DynamicCompsProvider<FIELD_COMPS>,
    DetectableAsField2,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    // /// Importing assignment operator from BaseTens
    // using Base::operator=;
    
    /// Aassign from another expression
    template <typename O>
    INLINE_FUNCTION
    void assign(O&& oth)
    {
#define LOOP(LOOP_TYPE)				\
      LOOP_TYPE(0,this->nSites(),		\
	   CAPTURE(self=this->getWritable(),	\
		   TO_READ(oth)),		\
	   site,				\
	   {							\
	     using RhsComps=typename std::decay_t<O>::Comps;	\
									\
	     /*! we need to take care that the rhs might not have the site (such in the case of a scalar) */ \
									\
	     if constexpr(tupleHasType<RhsComps,Site>)			\
	       self(site)=oth(site);					\
	     else							\
	       self(site)=oth;						\
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
    
    /// Copy assign
    INLINE_FUNCTION
    Field2& operator=(const Field2& oth)
    {
      assign(oth);
      
      return *this;
    }
    
    /// Assigns another node
    template <typename O>
    INLINE_FUNCTION
    Field2& operator=(const Node<O>& oth)
    {
      assign(*oth);
      
      return *this;
    }
    
    /// Move assign
    INLINE_FUNCTION
    Field2& operator=(Field2&& oth)
    {
      std::swap(data,oth.data);
      
      return *this;
    }
    
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
    
    /// Executes where allocated
    static constexpr MemoryType execSpace=MT;
    
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
    const Site nTotalAllocatedSites;
    
    /// Number of internal degrees of freedom - this will be made dynamic
    static constexpr int nInternalDegs=indexMaxValue<C...>();
    
    /// Storage data
    mutable Data data;
    
    /// Presence of halo and edges
    const HaloEdgesPresence haloEdgesPresence;
    
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
    
    /// Describe a field2
    void describe(const std::string& pref="") const
    {
      master_printf("%sField2 %s address %p\n",pref.c_str(),demangle(typeid(*this).name()).c_str(),this);
      master_printf("%s Components:\n",pref.c_str());
      (master_printf("%s  %s\n",pref.c_str(),demangle(typeid(C).name()).c_str()),...);
      master_printf("%s Site type: %s\n",pref.c_str(),demangle(typeid(decltype(this->nSites())).name()).c_str());
      master_printf("%s Fund: %s\n",pref.c_str(),demangle(typeid(_Fund).name()).c_str());
      master_printf("%s FieldCoverage: %d\n",pref.c_str(),FC);
      master_printf("%s FieldLayout: %d\n",pref.c_str(),FL);
      master_printf("%s MemoryType: %d\n",pref.c_str(),MT);
      master_printf("%s FieldLayout: %d\n",pref.c_str(),FL);
      master_printf("%s IsRef: %d\n",pref.c_str(),IsRef);
      master_printf("%sEnd of Field2\n",pref.c_str());
    }
    /////////////////////////////////////////////////////////////////
    
    /// Computes the squared norm, overloading default expression
    Fund norm2() const
    {
      Field2<CompsList<>,_Fund,FC,FL,MT> buf;
      
      PAR(0,this->nSites(),
	  CAPTURE(t=this->getReadable(),
		  TO_WRITE(buf)),
	  site,
	  {
	    buf(site)=t(site).norm2();
	  });
      
      Fund res;
      glb_reduce(&res,buf,this->nSites()());
      
      return res;
    }
    
    /////////////////////////////////////////////////////////////////
    
    void locReduce()
    {
      const Site nOri=nSites();
      Site n=nSites();
      //verbosity_lv2_master_printf("n: %d, nori: %d\n",n(),nOri());
      
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
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)						\
    Field2<CompsList<C...>,ATTRIB _Fund,FC,FL,MT,true> getRef() ATTRIB	\
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
      return *((Field2<CompsList<FlattenedInnerComp>,ATTRIB _Fund,FC,FL,MT,true>*)this); \
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
    
    /////////////////////////////////////////////////////////////////
    
// #define PROVIDE_SIMDIFY(ATTRIB)						\
//     INLINE_FUNCTION							\
//     auto simdify() ATTRIB						\
//     {									\
//       if constexpr(FL==FieldLayout::SIMDIFIABLE)			\
// 	return Field<CompsList<C...>,ATTRIB _Fund,L,LC,FieldLayout::SIMDIFIED,ES,true> \
// 	  (*lattice,haloPresence,(void*)data.storage,data.nElements,data.getDynamicSizes()); \
//       else								\
// 	{								\
// 	  using Traits=CompsListSimdifiableTraits<CompsList<C...>,_Fund>; \
// 	  								\
// 	  using SimdFund=typename Traits::SimdFund;			\
// 	  								\
// 	  return Field<typename Traits::Comps,ATTRIB SimdFund,L,LC,FieldLayout::SERIAL,ES,true> \
// 	    (*lattice,haloPresence,(ATTRIB void*)data.storage,data.nElements,data.getDynamicSizes()); \
// 	}								\
//     }
    
//     PROVIDE_SIMDIFY(const);
    
//     PROVIDE_SIMDIFY(/* non const */);
    
// #undef PROVIDE_SIMDIFY
    
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
    
    /// Create a field
    Field2(const HaloEdgesPresence& haloEdgesPresence=WITHOUT_HALO) :
      nTotalAllocatedSites(FieldSizes<fieldCoverage>::nSitesToAllocate(haloEdgesPresence)),
      data(std::make_tuple(nTotalAllocatedSites)),
      haloEdgesPresence(haloEdgesPresence),
      haloIsValid(false),
      edgesAreValid(false)
    {
      static_assert(not IsRef,"Can allocate only if not a reference");
      
      invalidateHalo();
    }
    
    /// Used to dispatch the copy constructor
    struct _CopyConstructInternalDispatcher;
    
    /// Copy constructor, internal implementation
    template <typename O>
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    Field2(O&& oth,
	   _CopyConstructInternalDispatcher*) :
      nTotalAllocatedSites(oth.nTotalAllocatedSites),
      data(oth.data),
      haloEdgesPresence(oth.haloEdgesPresence),
      haloIsValid(oth.haloIsValid),
      edgesAreValid(oth.edgesAreValid)
    {
    }
    
    /// Copy construct, taking as input a non-reference when this is a reference
    template <typename O,
	      bool B=IsRef,
	      ENABLE_THIS_TEMPLATE_IF(B and isField2<O>)>
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    Field2(O&& oth) :
      Field2(std::forward<O>(oth),(_CopyConstructInternalDispatcher*)nullptr)
    {
    }
    
    /// Copy constructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    Field2(const Field2& oth) :
      Field2(oth,(_CopyConstructInternalDispatcher*)nullptr)
    {
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
	      
	      ((Fund*)send_buf)[internalDeg+nInternalDegs*i()]=
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
	      
	      data(offset+i,c...)=
		((Fund*)recv_buf)[internalDeg+nInternalDegs*i()];
	    },dynamicSizes);
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
      
      assertCanCommunicate(Field2::nHaloSites());
      
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
      
      assertCanCommunicate(Field2::nEdgesSites());
      
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
