#ifndef _DYNAMICCOMPSPROVIDER_HPP
#define _DYNAMICCOMPSPROVIDER_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/dynamicCompsProvider.hpp

#include <expr/comps.hpp>
#include <metaprogramming/crtp.hpp>
#include <tuples/tupleDiscriminate.hpp>
#include <tuples/tupleSubset.hpp>

namespace nissa
{
  /// Provides the dynamic or static components
  template <typename C>
  struct DynamicCompsProvider
  {
    using DynamicStaticComps=
      TupleDiscriminate<SizeIsKnownAtCompileTime,C>;
    
    /// List of all statically allocated components
    using StaticComps=
      typename DynamicStaticComps::Valid;
    
    /// List of all dynamically allocated components
    using DynamicComps=
      typename DynamicStaticComps::Invalid;
    
    /// Number of dynamic components
    static constexpr int nDynamicComps=
      std::tuple_size_v<DynamicComps>;
    
    /// Check if the expression has any dynamic component
    static constexpr bool hasDynamicComps=
      nDynamicComps!=0;
    
    /// Returns dynamic comps from a tuple
    template <typename...T>
    INLINE_FUNCTION
    static DynamicComps filterDynamicComps(const CompsList<T...>& tc)
    {
      return tupleGetSubset<DynamicComps>(tc);
    }
    
    /// Returns dynamic comps from a list
    template <DerivedFromComp...T>
    INLINE_FUNCTION
    static DynamicComps filterDynamicComps(const T&...td)
    {
      return filterDynamicComps(std::make_tuple(td...));
    }
  };
  
  /// Provide no dynamics comp
  struct ProvideNoDynamicsComp
  {
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    CompsList<> getDynamicSizes() const
    {
      return {};
    }
  };
}

#endif
