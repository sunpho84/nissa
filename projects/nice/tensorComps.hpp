#ifndef _NICE_TENSOR_COMPS_HPP
#define _NICE_TENSOR_COMPS_HPP

#include "base/debug.hpp"
#include "new_types/su3.hpp"
#include "routines/math_routines.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace nissa::nice
{
  /// Index of a component definition in the runtime registry
  using CompKind=
    size_t;
  
  /// Definition shared by every occurrence of a registered component
  struct CompDefinition
  {
    /// Name used for registration and diagnostics
    std::string name;
    
    /// Fixed extent, or zero to use each occurrence's local extent
    size_t size{0};
  };
  
  /// Append-only component definitions with stable numeric identities
  struct CompRegistry
  {
    /// Definitions in registration order; existing entries must remain unchanged
    std::vector<CompDefinition> definitions;
    
    /// Reuse or register a definition, rejecting conflicting explicit sizes
    CompKind intern(const std::string& name,
                    const std::optional<size_t> size=std::nullopt)
    {
      if(name.empty())
        CRASH("NICE component kind requires a name");
      
      for(CompKind kind=0;kind<definitions.size();kind++)
        if(definitions[kind].name==name)
          {
            if(size and definitions[kind].size!=*size)
              CRASH("NICE component %s is already registered with a different size",name.c_str());
            
            return kind;
          }
      
      const CompKind kind=definitions.size();
      definitions.push_back({name,size.value_or(0)});
      
      return kind;
    }
    
    /// Resolve a registered component definition
    const CompDefinition& definition(const CompKind kind) const
    {
      if(kind>=definitions.size())
        CRASH("NICE component kind %zu is not registered",kind);
      
      return definitions[kind];
    }
    
    /// Resolve a component name for diagnostics
    const std::string& name(const CompKind kind) const
    {
      return definition(kind).name;
    }
  };
  
  /// Registry shared by initialization and subsequent script execution
  inline CompRegistry& compRegistry()
  {
    static CompRegistry registry;
    
    return registry;
  }
  
  /// Standard components registered through the same interface as script components
  namespace StandardComp
  {
    /// Local lattice-site index
    inline const CompKind Site=compRegistry().intern("Site",0);
    
    /// Dirac spin index
    inline const CompKind Spin=compRegistry().intern("Spin",NDIRAC);
    
    /// SU(3) color index
    inline const CompKind Color=compRegistry().intern("Color",NCOL);
    
    /// Four-dimensional Lorentz index
    inline const CompKind Lorentz=compRegistry().intern("Lorentz",NDIM);
    
    /// Independent copy index
    inline const CompKind Copy=compRegistry().intern("Copy",0);
    
    /// Time separation of a two-point result
    inline const CompKind Time=compRegistry().intern("Time",0);
    
    /// Index into the ordered gamma-pair list
    inline const CompKind GammaPair=compRegistry().intern("GammaPair",0);
    
    /// Real and imaginary parts of a complex value
    inline const CompKind Complex=compRegistry().intern("Complex",2);
  }
  
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
    
    /// Resolve the registered component name
    const std::string& kindName() const
    {
      return compRegistry().name(kind);
    }
  };
  
  /// One occurrence of a registered component in a tensor
  struct Comp :
    CompIdentity
  {
    /// Local extent, ignored when the registered definition has a fixed size
    size_t extent;
    
    /// Construct an occurrence with an optional local extent and tag
    Comp(const CompKind kind,
         const CompSide side=CompSide::None,
         const size_t extent=0,
         const size_t label=0) :
      CompIdentity{kind,side,label},
      extent(extent)
    {
    }
    
    /// Resolve the fixed size or the occurrence's local extent
    size_t size() const
    {
      const size_t fixed=compRegistry().definition(kind).size;
      
      return fixed!=0?fixed:extent;
    }
    
    /// Check registration and the occurrence's orientation
    void validate() const
    {
      (void)compRegistry().definition(kind);
      
      if(side!=CompSide::None and side!=CompSide::Row and side!=CompSide::Column)
        CRASH("NICE component %s has an unknown orientation",kindName().c_str());
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
    /// Components in storage order, including Complex when present
    std::vector<Comp> comps;
    
    /// Validate each component and reject repeated kind, orientation, and label
    void validate() const
    {
      for(size_t i=0;i<comps.size();i++)
        {
          const Comp& comp=comps[i];
          
          comp.validate();
          
          for(size_t earlier=0;earlier<i;earlier++)
            if(comps[earlier].sameAs(comp))
              CRASH("NICE tensor repeats component %s with the same orientation and label",
                    comp.kindName().c_str());
        }
    }
    
    /// Number of stored real values, checking multiplication overflow
    size_t nElements() const
    {
      validate();
      
      for(const Comp& comp : comps)
        if(comp.size()==0)
          return 0;
      
      size_t count=1;
      
      for(const Comp& comp : comps)
        count=checkedMultiply(count,comp.size(),"NICE tensor element count overflow");
      
      return count;
    }
  };
}

#endif
