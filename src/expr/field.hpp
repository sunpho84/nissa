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
  DECLARE_UNTRANSPOSABLE_COMP(Parity,int,2,parity);
  
  DECLARE_UNTRANSPOSABLE_COMP(LocLxSite,int,0,locLxSite);
  DECLARE_UNTRANSPOSABLE_COMP(LocEoSite,int,0,locEoSite);
  DECLARE_UNTRANSPOSABLE_COMP(LocEvnSite,int,0,locEvnSite);
  DECLARE_UNTRANSPOSABLE_COMP(LocOddSite,int,0,locOddSite);
  
  /// Has or not the halo and the edges
  /// Actual components to be exposed for a field
  ///
  /// Forward declaration of internal implementation
  template <typename TP,
	    typename F,
	    SitesCoverage SC,
	    FieldLayout FL>
  struct FieldCompsProvider;
  
#define PROVIDE_FIELD_COMPS_PROVIDER(COVERAGE,LAYOUT,SITE,TYPES...)	\
  									\
  template <typename...C,						\
	    typename F>							\
  struct FieldCompsProvider<CompsList<C...>,F,				\
			    SitesCoverage::COVERAGE,			\
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
  
  PROVIDE_DETECTABLE_AS(Field2);
  
#define FIELD_COMPS_PROVIDER FieldCompsProvider<CompsList<C...>,_Fund,SC,FL>
  
#define FIELD_COMPS typename FIELD_COMPS_PROVIDER::Comps
  
#define THIS						\
  Field2<CompsList<C...>,_Fund,SC,FL,MT,IsRef>
  
#define BASE					\
  Node<THIS>
  
  /// Defines a field
  template <typename...C,
	    typename _Fund,
	    SitesCoverage SC,
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
    
    /// Importing assignment operator from BaseTens
    using Base::operator=;
    
    /// Copy assign
    INLINE_FUNCTION
    Field2& operator=(const Field2& oth)
    {
#warning crash("critical");
      // Base::operator=(oth);
      
      return *this;
    }
    
    /// Move assign
    INLINE_FUNCTION
    Field2& operator=(Field2&& oth)
    {
      std::swap(data,oth.data);
      
      return *this;
    }
    
    static constexpr SitesCoverage sitesCoverage=SC;
    
    static constexpr FieldLayout fieldLayout=FL;
    
    /// Components
    using Comps=
      FIELD_COMPS;
    
    /// Import dynamic comps
    using DynamicComps=
      typename DynamicCompsProvider<FIELD_COMPS>::DynamicComps;
    
    /// Type used for the site
    using Site=typename FIELD_COMPS_PROVIDER::Site;
    
#undef FIELD_COMPS

#undef FIELD_COMPS_PROVIDER
    
    /// Fundamental tye
    using Fund=
      typename FieldCompsProvider<CompsList<C...>,_Fund,SC,FL>::Fund;
    
    /// Internal storage type
    using Data=
      DynamicTens<Comps,Fund,MT,IsRef>;
    
    /// Executes where allocated
    static constexpr MemoryType execSpace=MT;
    
    /////////////////////////////////////////////////////////////////
    
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static void assertHasDefinedCoverage()
    {
      static_assert(sitesCoverage==FULL_SPACE or sitesCoverage==EVEN_SITES or sitesCoverage==ODD_SITES,"Trying to probe some feature of a field with unknwon coverage. If you are accessing an EvnOrOddField subfield, do it without subscribing them, or cast to the specific subspace");
    }
    
    /// Number of sites covered by the field
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr const Site& nSites()
    {
      if constexpr(sitesCoverage==FULL_SPACE)
	return locVol;
      else
	return locVolh;
    }
    
    /// Number of sites in the halo of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr const Site& nHaloSites()
    {
      if constexpr(sitesCoverage==FULL_SPACE)
	return bord_vol;
      else
	return bord_volh;
    }
    
    /// Number of sites in the edges of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr const Site& nEdgesSites()
    {
      if constexpr(sitesCoverage==FULL_SPACE)
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
    static auto surfSiteOfHaloSite(const Site& iHalo)
    {
      assertHasDefinedCoverage();
      
      if constexpr(sitesCoverage==FULL_SPACE)
	return surflxOfBordlx[iHalo];
      else
	if constexpr(sitesCoverage==EVEN_SITES or sitesCoverage==ODD_SITES)
	  return surfeo_of_bordeo[sitesCoverage][iHalo];
    }
    
    /// Surface site of a site in the e
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static auto surfSiteOfEdgeSite(const Site& iEdge)
    {
      assertHasDefinedCoverage();
      
      if constexpr(sitesCoverage==FULL_SPACE)
	return surflxOfEdgelx[iEdge];
      else
	if constexpr(sitesCoverage==EVEN_SITES or sitesCoverage==ODD_SITES)
	  return surfeo_of_edgeo[sitesCoverage][iEdge];
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Storage data
    Data data;
    
    /// Presence of halo and edges
    const HaloEdgesPresence haloEdgesPresence;
    
    /// Total allocated sites
    const Site externalSize;
    
    /// States whether the halo is updated
    mutable bool haloIsValid;
    
    /// States whether the edges are updated
    mutable bool edgesAreValid;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    auto getDynamicSizes() const
    {
      return std::make_tuple(externalSize);
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
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)						\
    Field2<CompsList<C...>,ATTRIB _Fund,SC,FL,MT,true> getRef() ATTRIB	\
    {									\
      return *this;							\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /// Gets a writeable reference
    constexpr INLINE_FUNCTION
    auto getWritable()
    {
      return this->getRef();
    }
    
    /// Gets a read-only reference
    constexpr INLINE_FUNCTION
    auto getReadable() const
    {
      return this->getRef();
    }
    
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
      data(std::make_tuple(externalSize)),
      haloEdgesPresence(haloEdgesPresence),
      externalSize(FieldSizes<sitesCoverage>::nSitesToAllocate(haloEdgesPresence)),
      haloIsValid(false),
      edgesAreValid(false)
    {
      static_assert(not IsRef,"Can allocate only if not a reference");
      
      invalidateHalo();
    }
    
    /// Create a refence to a field
    template <typename O,
	      bool B=IsRef,
	      ENABLE_THIS_TEMPLATE_IF(B and isField2<O>)>
    Field2(O&& oth) :
      data(oth.data),
      haloEdgesPresence(oth.haloEdgesPresence),
      externalSize(oth.externalSize),
      haloIsValid(oth.haloIsValid),
      edgesAreValid(oth.edgesAreValid)
    {
    }
    
    // /// Copy constructor
    // INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    // Field(const Field& oth) :
    //   lattice(oth.lattice),
    //   data(oth.data),
    //   haloPresence(oth.haloPresence),
    //   haloStatus(oth.haloStatus)
    // {
    // }
    
    
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
  };
}

#endif
