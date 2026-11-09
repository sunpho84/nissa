#ifndef _NICE_GAUGE_RECIPES_HPP
#define _NICE_GAUGE_RECIPES_HPP

#include "configurationSources.hpp"
#include "base/debug.hpp"
#include "base/field.hpp"
#include "io/reader.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/APE.hpp"
#include "rngLayout.hpp"
#include "runOptions.hpp"

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
  
  /// Index of a deferred producer and its per-configuration field slot
  /// The index refers to a recipe, never directly to a prepared gauge field
  using GaugeRecipeId=
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
    
    /// Report this producer's double count for one configuration with the effective run settings
    std::function<uint64_t(const RunOptions&)> rngDoubles;
    
    /// Produce a field using the selected entry and already filled earlier slots
    std::function<std::unique_ptr<GaugeField>(const ConfigurationEntry&,
                                              const GaugeRecipes&,
                                              const ActionRngRange&)> produce;
  };
  
  /// Resolve an already produced field, checking bounds and availability
  inline const GaugeField& gaugeField(const GaugeRecipes& recipes,
                                     const GaugeRecipeId recipeId)
  {
    if(recipeId>=recipes.size() or not recipes[recipeId].field)
      CRASH("NICE gauge recipe %zu has no field for this configuration",recipeId);
    
    return *recipes[recipeId].field;
  }
  
  /// Append a deferred producer and return its target index
  inline GaugeRecipeId addGaugeRecipe(GaugeRecipes& recipes,
                                    std::function<std::unique_ptr<GaugeField>(const ConfigurationEntry&,
                                                                              const GaugeRecipes&,
                                                                              const ActionRngRange&)> produce,
                                    std::function<uint64_t(const RunOptions&)> rngDoubles={})
  {
    if(not produce)
      CRASH("NICE gauge recipe requires a producer");
    
    const GaugeRecipeId target=recipes.size();
    recipes.push_back({nullptr,std::move(rngDoubles),std::move(produce)});
    
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
  inline GaugeRecipeId loadedConf(GaugeRecipes& recipes)
  {
    return addGaugeRecipe(recipes,[](const ConfigurationEntry& entry,
                                    const GaugeRecipes&,
                                    const ActionRngRange&)
    {
      if(entry.path.empty())
        CRASH("NICE loadedConf requires an input configuration path");
      
      auto field=std::make_unique<GaugeField>("nice loaded conf",WITH_HALO_EDGES);
      read_ildg_gauge_conf(*field,entry.path);
      
      return field;
    });
  }
  
  /// Record spatial APE smearing of a previously registered target
  inline GaugeRecipeId apeSmearConf(GaugeRecipes& recipes,
                                  const GaugeRecipeId inputRecipe,
                                  const double alpha,
                                  const int n)
  {
    if(inputRecipe>=recipes.size())
      CRASH("NICE input gauge recipe is outside the registered recipes");
    
    if(not std::isfinite(alpha) or n<0)
      CRASH("NICE apeSmearConf requires finite alpha and nonnegative n");
    
    return addGaugeRecipe(recipes,[inputRecipe,alpha,n](const ConfigurationEntry&,
                                                       const GaugeRecipes& current,
                                                       const ActionRngRange&)
    {
      const auto& original=gaugeField(current,inputRecipe);
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
                                 const ConfigurationEntry& entry,
                                 const Sitmo::Rng& rng,
                                 const RngLayout& rngLayout,
                                 const uint64_t iConfiguration)
  {
    validateGaugeRecipes(recipes);
    rngLayout.validate();
    
    if(not lxGeomInited)
      CRASH("NICE gauge preparation requires initialized geometry");
    
    if(rngLayout.gaugeActions.size()!=recipes.size())
      CRASH("NICE RNG layout does not match the registered gauge recipes");
    
    // Clear every slot before the first producer, preventing stale inputs from a previous entry
    clearGaugeFields(recipes);
    
    for(size_t slot=0;slot<recipes.size();slot++)
      {
        auto& recipe=recipes[slot];
        const ActionRngRange range=
          rngLayout.actionRange(rng,rngLayout.gaugeActions[slot],iConfiguration);
        recipe.field=
          recipe.produce(entry,recipes,range);
        
        if(not recipe.field)
          CRASH("NICE gauge recipe %zu returned no field",slot);
      }
  }
  
  /// Expose gauge recipe IDs and deferred builders to the script evaluator
  /// The recipe collection must outlive the registered functions
  void registerGaugeRecipes(Evaluator& evaluator,
                            GaugeRecipes& recipes);
}

#endif
