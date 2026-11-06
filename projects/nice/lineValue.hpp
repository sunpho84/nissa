#ifndef _NICE_LINE_VALUE_HPP
#define _NICE_LINE_VALUE_HPP

#include "runOptions.hpp"
#include "sField.hpp"

#include <cstddef>
#include <optional>
#include <variant>
#include <vector>

namespace nissa::nice
{
  /// Realized source or propagator with the policy that determines its outer components
  struct LineValue
  {
    /// Independent copies retained in the line
    size_t nCopies;
    
    /// Whether the source spin column is present
    bool diluteSpin;
    
    /// Whether the source color column is present
    bool diluteColor;
    
    /// Field owning every copy and diluted source component
    SField field;
    
    /// Capture the effective run settings for one line realization
    LineValue(const RunOptions& options,
              const MemorySpace space=defaultMemorySpace,
              const SpaceTimeLayout layout=defaultSpaceTimeLayout) :
      nCopies(checkedCopies(options)),
      diluteSpin(options.diluteSpin),
      diluteColor(options.diluteColor),
      field(makeLineField(nCopies,diluteSpin,diluteColor,space,layout))
    {
    }
    
    /// Convert the run count only after checking that it is positive
    static size_t checkedCopies(const RunOptions& options)
    {
      if(options.nCopies<=0)
        CRASH("NICE line requires a positive number of copies");
      
      return static_cast<size_t>(options.nCopies);
    }
    
    /// Check the policy needed to combine or extend two line values
    bool sameShapeAs(const LineValue& other) const
    {
      return nCopies==other.nCopies and diluteSpin==other.diluteSpin and
        diluteColor==other.diluteColor and field.nSites==other.field.nSites;
    }
  };
  
  /// One general tensor, non-line lattice field, or realized line
  using TensorValue=
    std::variant<Tensor,SField,LineValue>;
  
  /// Per-hit owned results indexed by stable task identity
  using TaskValues=
    std::vector<std::optional<TensorValue>>;
}

#endif
