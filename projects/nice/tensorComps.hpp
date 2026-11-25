#ifndef _NICE_TENSOR_COMPS_HPP
#define _NICE_TENSOR_COMPS_HPP

#include "base/debug.hpp"
#include "new_types/su3.hpp"
#include "routines/math_routines.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace nissa::nice
{
  /// Numerical scalar carried at each combination of tensor indices
  enum class ScalarKind{Real,Complex};
  
  /// Numeric identity of a component kind within one compiled NICE program
  using CompKind=
    size_t;
  
  /// Fixed identities for the components used by the initial NISSA kernels
  namespace StandardComp
  {
    /// Extent supplied by geometry, a run setting, or a script declaration
    inline constexpr size_t ContextDependentSize=0;
    
    /// Compile-time extent of a standard kind, or ContextDependentSize
    template<CompKind Kind>
    inline constexpr size_t compSize=ContextDependentSize;
    
    /// Whether a standard kind requires row or column orientation
    template<CompKind Kind>
    inline constexpr bool compRequiresOrientation=false;
    
    /// Define one fixed kind together with its index, extent, and orientation rule
#define PROVIDE_COMP_KIND(NAME,SIZE,INDEX,ORIENTED) \
    inline constexpr CompKind NAME=INDEX; \
    template<> \
    inline constexpr size_t compSize<NAME> =SIZE; \
    template<> \
    inline constexpr bool compRequiresOrientation<NAME> =ORIENTED
    
    /// Local lattice-site index
    PROVIDE_COMP_KIND(Site,ContextDependentSize,0,false);
    
    /// Dirac spin index
    PROVIDE_COMP_KIND(Spin,NDIRAC,1,true);
    
    /// SU(3) color index
    PROVIDE_COMP_KIND(Color,NCOL,2,true);
    
    /// Four-dimensional Lorentz index
    PROVIDE_COMP_KIND(Lorentz,NDIM,3,false);
    
    /// Independent copy index
    PROVIDE_COMP_KIND(Copy,ContextDependentSize,4,false);
    
    /// Time separation of a two-point result
    PROVIDE_COMP_KIND(Time,ContextDependentSize,5,false);
    
    /// Index into the ordered gamma-pair list of a contraction
    PROVIDE_COMP_KIND(GammaPair,ContextDependentSize,6,false);
    
#undef PROVIDE_COMP_KIND
    
    /// Number of registered standard comps
    inline constexpr CompKind NStandardComps=GammaPair+1;
    
    /// Sentinel that cannot identify a registered component
    inline constexpr CompKind Invalid=std::numeric_limits<CompKind>::max();
    
    /// Fixed extents indexed by the numeric identity of each standard kind
    inline constexpr std::array<size_t,NStandardComps> compSizes=
      []<size_t... I>(std::index_sequence<I...>)
      {
        return std::array<size_t,sizeof...(I)>{compSize<I>...};
      }(std::make_index_sequence<NStandardComps>{});
    
    /// Orientation requirements indexed by standard component kind
    inline constexpr std::array<bool,NStandardComps> compOrientationRequired=
      []<size_t... I>(std::index_sequence<I...>)
      {
        return std::array<bool,sizeof...(I)>{compRequiresOrientation<I>...};
      }(std::make_index_sequence<NStandardComps>{});
  }
  
  /// Assign numeric component kinds while retaining their names for diagnostics
  struct CompRegistry
  {
    /// Names in component-kind order, including the fixed standard kinds
    std::vector<std::string> names{"Site","Spin","Color","Lorentz","Copy","Time","GammaPair"};
    
    /// Reuse an existing kind or assign a new one at script compilation time
    CompKind intern(const std::string& name)
    {
      if(name.empty())
        CRASH("NICE component kind requires a name");
      
      for(CompKind kind=0;kind<names.size();kind++)
        if(names[kind]==name)
          return kind;
      
      const CompKind kind=names.size();
      names.push_back(name);
      
      return kind;
    }
    
    /// Resolve a kind back to its name for script diagnostics
    const std::string& name(const CompKind kind) const
    {
      if(kind>=names.size())
        CRASH("NICE component kind %zu is not registered",kind);
      
      return names[kind];
    }
  };
  
  /// Orientation used when determining products and contractions
  enum class CompSide{None,Row,Column};
  
  /// Identity shared by a component and a coordinate selecting it
  struct CompIdentity
  {
    /// Numeric identity of its index space in the program's component registry
    CompKind kind;
    
    /// Row, column, or an unoriented index
    CompSide side;
    
    /// Distinguish repetitions of the same kind and orientation
    size_t label{0};
    
    /// Compare the identities without comparing extent or coordinate
    bool sameAs(const CompIdentity& other) const
    {
      return kind==other.kind and side==other.side and label==other.label;
    }
    
    /// Resolve the component name through its registry when available
    std::string kindName(const CompRegistry* registry=nullptr) const
    {
      if(registry!=nullptr)
        return registry->name(kind);
      
      static const CompRegistry standardRegistry;
      
      if(kind<standardRegistry.names.size())
        return standardRegistry.name(kind);
      
      return "kind "+std::to_string(kind);
    }
  };
  
  /// One interned tensor component with its local extent
  struct Comp :
    CompIdentity
  {
    /// Number of locally stored values along this component
    size_t extent;
    
    /// Keep the existing kind, side, extent, label construction order
    Comp(const CompKind kind,
         const CompSide side,
         const size_t extent,
         const size_t label=0) :
      CompIdentity{kind,side,label},
      extent(extent)
    {
    }
    
    /// Validate this component independently of other tensor components
    void validate(const CompRegistry* registry=nullptr) const
    {
      if(kind==StandardComp::Invalid)
        CRASH("NICE component has an invalid kind");
      
      if(registry!=nullptr and kind>=registry->names.size())
        CRASH("NICE component kind %zu is not registered",kind);
      
      if(side!=CompSide::None and side!=CompSide::Row and side!=CompSide::Column)
        CRASH("NICE component %s has an unknown orientation",kindName(registry).c_str());
      
      if(kind<StandardComp::compSizes.size() and
         StandardComp::compSizes[kind]!=StandardComp::ContextDependentSize and
         extent!=StandardComp::compSizes[kind])
        CRASH("NICE component %s has the wrong fixed extent",kindName(registry).c_str());
      
      if(kind<StandardComp::compOrientationRequired.size() and
         (side!=CompSide::None)!=StandardComp::compOrientationRequired[kind])
        CRASH("NICE component %s has an invalid orientation",kindName(registry).c_str());
    }
  };
  
  /// Coordinate of one component identified by kind, orientation, and label
  struct CompIndex :
    CompIdentity
  {
    /// Coordinate within the selected component
    size_t coordinate;
    
    /// Keep the kind, side, coordinate, label selection order
    CompIndex(const CompKind kind,
              const CompSide side,
              const size_t coordinate,
              const size_t label=0) :
      CompIdentity{kind,side,label},
      coordinate(coordinate)
    {
    }
  };
  
  /// Shape and storage order of one tensor
  struct TensorDescriptor
  {
    /// Type of each tensor element
    ScalarKind scalar;
    
    /// Components in tensor storage order
    std::vector<Comp> comps;
    
    /// Adjacent logical elements per real or imaginary plane; one means interleaved
    size_t scalarComponentStride{1};
    
    /// Validate each component and reject repeated kind, orientation, and label
    void validate(const CompRegistry* registry=nullptr) const
    {
      if(scalarComponentStride==0 or
         (scalar==ScalarKind::Real and scalarComponentStride!=1))
        CRASH("NICE tensor has an invalid scalar component stride");
      
      for(size_t i=0;i<comps.size();i++)
        {
          const Comp& comp=comps[i];
          
          comp.validate(registry);
          
          for(size_t earlier=0;earlier<i;earlier++)
            if(comps[earlier].sameAs(comp))
              CRASH("NICE tensor repeats component %s with the same orientation and label",
                    comp.kindName(registry).c_str());
        }
    }
    
    /// Number of locally stored scalar values, checking multiplication overflow
    size_t nElements() const
    {
      validate();
      
      for(const Comp& comp : comps)
        if(comp.extent==0)
          return 0;
      
      size_t count=1;
      
      for(const Comp& comp : comps)
        count=checkedMultiply(count,comp.extent,"NICE tensor element count overflow");
      
      return count;
    }
  };
  
}

#endif
