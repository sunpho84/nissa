#ifndef _EOFIELD_HPP
#define _EOFIELD_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <expr/eoFieldDeclaration.hpp>

namespace nissa
{
  /// Specifies the order of components
  template <typename TP,
	    typename F,
	    FieldLayout FL>
  struct EoField2CompsProvider;
  
#define PROVIDE_EO_FIELD2_COMPS_PROVIDER(LAYOUT,SITE,TYPES...)		\
  									\
  template <typename...C,						\
	    typename F>							\
  struct EoField2CompsProvider<CompsList<C...>,F,			\
			       FieldLayout::LAYOUT>			\
  {									\
    using Comps=CompsList<TYPES>;					\
    									\
    using Site=SITE;							\
    									\
    using Fund=F;							\
  }
  
  PROVIDE_EO_FIELD2_COMPS_PROVIDER(CPU,LocEoSite,Parity,LocEoSite,C...);
  PROVIDE_EO_FIELD2_COMPS_PROVIDER(GPU,LocEoSite,Parity,C...,LocEoSite);
  
#undef PROVIDE_EO_FIELD2_COMPS_PROVIDER
  
  /////////////////////////////////////////////////////////////////
  
  PROVIDE_DETECTABLE_AS(EoField2);
  
#define EO_FIELD2_COMPS_PROVIDER EoField2CompsProvider<CompsList<C...>,_Fund,FL>
  
#define FIELD_COMPS typename EO_FIELD2_COMPS_PROVIDER::Comps
  
#define THIS						\
  EoField2<CompsList<C...>,_Fund,FL,MT,IsRef>
  
#define BASE					\
  Node<THIS>
  
  /// Structure to hold an even/old field
  template <typename...C,
	    typename _Fund,
	    FieldLayout FL,
	    MemoryType MT,
	    bool IsRef>
  struct THIS :
    DynamicCompsProvider<FIELD_COMPS>,
    DetectableAsEoField2,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    /// Type representing a pointer to type T
    template <FieldCoverage EO>
    using _F=Field2<CompsList<C...>,_Fund,EO,FL,MT,IsRef>;
    
    using Fevn=_F<EVEN_SITES>;
    
    using Fodd=_F<ODD_SITES>;
    
    Fevn evenPart;
    
    Fodd oddPart;
    
    static constexpr FieldLayout fieldLayout=FL;
    
    /// Components
    using Comps=
      FIELD_COMPS;
    
    /// Import dynamic comps
    using DynamicComps=
      typename DynamicCompsProvider<FIELD_COMPS>::DynamicComps;
    
    /// Type used for the site
    using Site=typename EO_FIELD2_COMPS_PROVIDER::Site;
    
    /// Fundamental tye
    using Fund=typename EO_FIELD2_COMPS_PROVIDER::Fund;
    
#undef FIELD_COMPS
    
#undef EO_FIELD2_COMPS_PROVIDER
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    auto getDynamicSizes() const
    {
      return std::apply([](const auto&...c)
      {
	return std::make_tuple(compCast<LocEoSite,LocEvnSite>(c)...);
      },evenPart.data.getDynamicSizes());
    }
    
    /// Provide evaluator
#define PROVIDE_EVAL(ATTRIB)						\
    template <typename...U>						\
    CUDA_HOST_AND_DEVICE constexpr INLINE_FUNCTION			\
    ATTRIB Fund& eval(const U&...cs) ATTRIB				\
    {									\
      const auto c=std::make_tuple(cs...);				\
      									\
      const Parity parity=std::get<Parity>(c);				\
      									\
      const Site site=std::get<Site>(c);				\
      									\
      if(parity==0)							\
	return evenPart(locEvnSite(site()),std::get<C...>(c));		\
      else								\
	return oddPart(locOddSite(site()),std::get<C...>(c));		\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/* non const */);
    
#undef PROVIDE_EVAL
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=
      not std::is_const_v<Fund>;
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(CONST_ATTRIB)					\
    INLINE_FUNCTION constexpr						\
    EoField2<CompsList<C...>,CONST_ATTRIB _Fund,FL,MT,true> getRef() CONST_ATTRIB\
    {							\
      return {evenPart.getRef(),oddPart.getRef()};	\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* not const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
    /// We keep referring to the original object
    static constexpr bool storeByRef=not IsRef;
    
    /// Returns that can assign
    constexpr INLINE_FUNCTION
    bool canAssign()
    {
      return canAssignAtCompileTime;
    }
    
    /// Constructor
    EoField2(Fevn&& ev,
	     Fodd&& od) :
      evenPart(ev),
      oddPart(od)
    {
    }
    
    /// Constructor
    EoField2(const HaloEdgesPresence& haloEdgesPresence=WITHOUT_HALO) :
      evenPart(haloEdgesPresence),
      oddPart(haloEdgesPresence)
    {
    }
    
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    ~EoField2()
    {
    }
    
    // /// Reset
    // INLINE_FUNCTION
    // void reset()
    // {
    //   forBothParities([this](const auto& par)
    //   {
    // 	(*this)[par].reset();
    //   });
    // }
    
    // /// Update the halo of both parities
    // INLINE_FUNCTION
    // void updateHalo() const
    // {
    //   forBothParities([this](const auto& par)
    //   {
    // 	(*this)[par].updateHalo();
    //   });
    // }
    
    // /// Update the edges of both parities
    // INLINE_FUNCTION
    // void updateEdges() const
    // {
    //   forBothParities([this](const auto& par)
    //   {
    // 	(*this)[par].updateEdges();
    //   });
    // }
    
    // /// Invalidate the halo of both parities
    // INLINE_FUNCTION
    // void invalidateHalo()
    // {
    //   forBothParities([this](const auto& par)
    //   {
    // 	(*this)[par].invalidateHalo();
    //   });
    // }
    
    // /// Invalidate the edges of both parities
    // INLINE_FUNCTION
    // void invalidateEdges()
    // {
    //   forBothParities([this](const auto& par)
    //   {
    // 	(*this)[par].invalidateEdges();
    //   });
    // }
    
    // /// Compare
    // INLINE_FUNCTION
    // bool operator==(const EoField& oth) const
    // {
    //   return evenPart==oth.evenPart and oddPart==oth.oddPart;
    // }
    
    // /// Negate comparison
    // INLINE_FUNCTION
    // bool operator!=(const EoField& oth) const
    // {
    //   return not (*this)==oth;
    // }
    
    // /// Assign
    // INLINE_FUNCTION
    // EoField& operator=(const EoField& oth)
    // {
    //   evenPart=oth.evenPart;
    //   oddPart=oth.oddPart;
      
    //   return *this;
    // }
  };
  
  // /// Loop on both parities
  // template <typename F>
  // INLINE_FUNCTION void forBothParities(F&& f)
  // {
  //   f(Par<0>{});
  //   f(Par<1>{});
  // }
}

#endif
