#ifndef _EOFIELD_HPP
#define _EOFIELD_HPP

#include <base/field.hpp>
#include <base/memory_manager.hpp>
#include <expr/field.hpp>

namespace nissa
{
  /// Structure to hold an even/old field
  template <typename C,
	    typename Fund,
	    FieldLayout FL=defaultFieldLayout,
	    MemoryType MT=defaultMemoryType,
	    typename Fevn=Field2<C,Fund,EVEN_SITES,FL,MT>,
	    typename Fodd=Field2<C,Fund,ODD_SITES,FL,MT>>
  struct EoField2
  {
    /// Type representing a pointer to type T
    template <SitesCoverage EO>
    using F=Field<T,EO,FL>;
    
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
	CONST Field<T,EVEN_OR_ODD_SITES,FL>;				\
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
	    FieldRef<Field<T,ODD_SITES,FL>>>
    getWritable()
    {
      return {evenPart,oddPart};
    }
    
    constexpr INLINE_FUNCTION
    EoField<T,FL,
	    FieldRef<const Field<T,EVEN_SITES,FL>>,
	    FieldRef<const Field<T,ODD_SITES,FL>>>
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
      return not (*this)==oth;
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
}

#endif
