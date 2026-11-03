#ifndef _NICE_GAUGE_RECIPES_HPP
#define _NICE_GAUGE_RECIPES_HPP

#include "configurationSources.hpp"
#include "base/debug.hpp"
#include "base/field.hpp"
#include "io/reader.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/APE.hpp"

#include <cmath>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

struct Evaluator;

namespace nissa::nice
{
  /// Gauge storage shared by loading and preparation operations
  using GaugeField=
    LxField<quad_su3>;
  
  /// Index of a slot containing both field storage and its producer
  /// All handles refer to the single ordered recipe collection of the current program
  using GaugeHandle=
    size_t;
  
  struct GaugeRecipe;
  
  /// Slots in dependency order, retained across configuration entries
  using GaugeRecipes=
    std::vector<GaugeRecipe>;
  
  /// Field storage paired with the recipe that fills it for the current configuration
  struct GaugeRecipe
  {
    /// Owned field, empty before preparation and after release
    std::unique_ptr<GaugeField> field;
    
    /// Produce a field using the selected entry and already filled earlier slots
    std::function<std::unique_ptr<GaugeField>(const ConfigurationEntry&,const GaugeRecipes&)> produce;
  };
  
  /// Resolve an already produced field, checking bounds and availability
  inline const GaugeField& gaugeField(const GaugeRecipes& recipes,
                                     const GaugeHandle handle)
  {
    if(handle>=recipes.size() or not recipes[handle].field)
      CRASH("NICE gauge field %zu is not available for this configuration",handle);
    
    return *recipes[handle].field;
  }
  
  /// Append a deferred producer and return its target index
  inline GaugeHandle addGaugeRecipe(GaugeRecipes& recipes,
                                    std::function<std::unique_ptr<GaugeField>(const ConfigurationEntry&,const GaugeRecipes&)> produce)
  {
    if(not produce)
      CRASH("NICE gauge recipe requires a producer");
    
    const GaugeHandle target=recipes.size();
    recipes.push_back({nullptr,std::move(produce)});
    
    return target;
  }
  
  /// Check callable presence after direct recipe edits
  inline void validateGaugeRecipes(const GaugeRecipes& recipes)
  {
    for(size_t slot=0;slot<recipes.size();slot++)
      if(not recipes[slot].produce)
        CRASH("NICE gauge recipe %zu has no producer",slot);
  }
  
  /// Record loading the selected entry's gauge file without performing I/O
  inline GaugeHandle loadedConf(GaugeRecipes& recipes)
  {
    return addGaugeRecipe(recipes,[](const ConfigurationEntry& entry,
                                    const GaugeRecipes&)
    {
      if(entry.path.empty())
        CRASH("NICE loadedConf requires an input configuration path");
      
      auto field=std::make_unique<GaugeField>("nice loaded conf",WITH_HALO_EDGES);
      read_ildg_gauge_conf(*field,entry.path);
      
      return field;
    });
  }
  
  /// Record spatial APE smearing of a previously registered target
  inline GaugeHandle apeSmearConf(GaugeRecipes& recipes,
                                  const GaugeHandle input,
                                  const double alpha,
                                  const int n)
  {
    if(input>=recipes.size())
      CRASH("NICE gauge handle is outside the registered recipes");
    
    if(not std::isfinite(alpha) or n<0)
      CRASH("NICE apeSmearConf requires finite alpha and nonnegative n");
    
    return addGaugeRecipe(recipes,[input,alpha,n](const ConfigurationEntry&,
                                                const GaugeRecipes& current)
    {
      const auto& original=gaugeField(current,input);
      auto field=std::make_unique<GaugeField>("nice APE conf",WITH_HALO_EDGES);
      ape_spatial_smear_conf(*field,original,alpha,n);
      
      return field;
    });
  }
  
  /// Release the current configuration's fields while retaining all producers
  inline void clearGaugeFields(GaugeRecipes& recipes)
  {
    for(auto& recipe : recipes)
      recipe.field.reset();
  }
  
  /// Clear previous fields and fill each slot for one entry after geometry setup
  inline void prepareGaugeFields(GaugeRecipes& recipes,
                                 const ConfigurationEntry& entry)
  {
    validateGaugeRecipes(recipes);
    
    if(not lxGeomInited)
      CRASH("NICE gauge preparation requires initialized geometry");
    
    // Clear every slot before the first producer, preventing stale inputs from a previous entry
    clearGaugeFields(recipes);
    
    for(size_t slot=0;slot<recipes.size();slot++)
      {
        auto& recipe=recipes[slot];
        recipe.field=recipe.produce(entry,recipes);
        
        if(not recipe.field)
          CRASH("NICE gauge recipe %zu returned no field",slot);
      }
  }
  
  /// Expose gauge handles and deferred recipe builders to the script evaluator
  /// The recipe collection must outlive the registered functions
  void registerGaugeRecipes(Evaluator& evaluator,
                            GaugeRecipes& recipes);
}

#endif
