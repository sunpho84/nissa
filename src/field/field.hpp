#ifndef _FIELD_HPP
#define _FIELD_HPP

#include "component.hpp"
#include "tensor.hpp"

namespace nissa
{
  template <typename T>
  struct SpaceTimeProps
  {
    inline static T loc();
  };
  
#define PROVIDE_LOC_PROP(TYPE,VAL)		\
  template <>					\
  inline TYPE SpaceTimeProps<TYPE>::loc()	\
  {						\
    return TYPE{VAL};				\
  }
  
  PROVIDE_LOC_PROP(LocVolIdx,loc_vol)
  PROVIDE_LOC_PROP(LocVolEvnIdx,loc_volh)
  PROVIDE_LOC_PROP(LocVolOddIdx,loc_volh)
  
  /// Field type
  template <typename Tc,
	    typename F=double>
  struct Field : public Subscribable<Field<Tc,F>>
  {
    using CT=CompsTraits<Field<Tc,F>>;
    
    /// Storing data
    Tens<Tc,F> data;
    
    /// Create taking the dynamical sizes as argument
    template <typename...S>
    Field(S&&...s) :
      data(SpaceTimeProps<typename CT::SpaceTimeComp>::loc(),std::forward<S>(s)...)
    {
    }
    
    /// Evaluate
    template <typename...A>
    const F& eval(A&&...a) const
    {
      return data(std::forward<A>(a)...);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD(eval);
  };
  
  /// Properties of components of field
  template <typename Tc,
	    typename F>
  struct CompsTraits<Field<Tc,F>>
  {
    /// Actual type
    using Type=Field<Tc,F>;
    
    /// List of valid types
    using SpaceTimeTypes=std::tuple<LocVolIdx,LocVolEvnIdx,LocVolOddIdx>;
    
    /// Components
    using Comps=Tc;
    
    /// Search the spacetime
    using _SpaceTimeComps=TupleCommonTypes<SpaceTimeTypes,Tc>;
    
    static_assert(std::tuple_size<_SpaceTimeComps>::value==1,"More or no spacetime component present");
    
    /// Detect the spacetime type
    using SpaceTimeComp=std::remove_reference_t<decltype(std::get<0>(*(_SpaceTimeComps*)nullptr))>;
  };
}

#endif
