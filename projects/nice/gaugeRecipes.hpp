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

#include <bit>
#include <cmath>
#include <format>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

struct Evaluator;

namespace nissa::nice
{
  /// Gauge storage shared by loading and preparation operations
  using GaugeField=
    LxField<quad_su3>;
  
  /// Index of a deferred producer and its per-gauge configuration field slot
  /// The index refers to a recipe, never directly to a prepared gauge field
  using GaugeRecipeId=
    size_t;
  
  struct GaugeRecipe;
  
  /// Slots in dependency order, retained across gauge configuration entries
  using GaugeRecipes=
    std::vector<GaugeRecipe>;
  
  /// Field storage paired with the recipe that fills it for the current gauge configuration
  struct GaugeRecipe
  {
    /// Owned field, empty before preparation and after release
    std::unique_ptr<GaugeField> field;
    
    /// Stable description of the producer and its parameters for checkpoint compatibility
    std::string checkpointSignature;
    
    /// Report this producer's uniform-double consumption cost for one gauge configuration
    std::function<uint64_t(const RunOptions&)> doubleDrawCost;
    
    /// Produce a field using the selected entry and already filled earlier slots
    std::function<std::unique_ptr<GaugeField>(const ConfigurationEntry&,
                                              const GaugeRecipes&,
                                              const OperationRngRange&)> produce;
  };
  
  /// Resolve an already produced field, checking bounds and availability
  inline const GaugeField& gaugeField(const GaugeRecipes& recipes,
                                     const GaugeRecipeId recipeId)
  {
    if(recipeId>=recipes.size() or not recipes[recipeId].field)
      CRASH("NICE gauge recipe %zu has no field for this gauge configuration",recipeId);
    
    return *recipes[recipeId].field;
  }
  
  /// Append a deferred producer and return its target index
  inline GaugeRecipeId addGaugeRecipe(GaugeRecipes& recipes,
                                    std::string checkpointSignature,
                                    std::function<std::unique_ptr<GaugeField>(const ConfigurationEntry&,
                                                                              const GaugeRecipes&,
                                                                              const OperationRngRange&)> produce,
                                    std::function<uint64_t(const RunOptions&)> doubleDrawCost={})
  {
    if(not produce or checkpointSignature.empty())
      CRASH("NICE gauge recipe requires a producer and checkpoint signature");
    
    const GaugeRecipeId target=recipes.size();
    recipes.push_back({nullptr,std::move(checkpointSignature),
                       std::move(doubleDrawCost),std::move(produce)});
    
    return target;
  }
  
  /// Check callable presence after direct recipe edits
  inline void validateGaugeRecipes(const GaugeRecipes& recipes)
  {
    for(size_t slot=0;slot<recipes.size();slot++)
      if(not recipes[slot].produce or recipes[slot].checkpointSignature.empty())
        CRASH("NICE gauge recipe %zu has no producer or checkpoint signature",slot);
  }
  
  /// Record loading the selected entry's gauge file without performing I/O
  inline GaugeRecipeId loadedConf(GaugeRecipes& recipes)
  {
    return addGaugeRecipe(recipes,"loadedConf",[](const ConfigurationEntry& entry,
                                    const GaugeRecipes&,
                                    const OperationRngRange&)
    {
      if(entry.path.empty())
        CRASH("NICE loadedConf requires an input gauge configuration path");
      
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
    
    return addGaugeRecipe(recipes,
      std::format("apeSmearConf:{}:{}:{}",inputRecipe,std::bit_cast<uint64_t>(alpha),n),
      [inputRecipe,alpha,n](const ConfigurationEntry&,
                                                       const GaugeRecipes& current,
                                                       const OperationRngRange&)
    {
      const auto& original=gaugeField(current,inputRecipe);
      auto field=std::make_unique<GaugeField>("nice APE conf",WITH_HALO_EDGES);
      ape_spatial_smear_conf(*field,original,alpha,n);
      
      return field;
    });
  }
  
  /// Release the current gauge configuration's fields while retaining all producers
  inline void clearGaugeFields(GaugeRecipes& recipes)
  {
    for(auto& recipe : recipes)
      recipe.field.reset();
  }
  
  /// Fill gauge fields in dependency order, checking a claim after the initial producer
  /// Return false after releasing the initial field when the claim was lost
  inline bool prepareGaugeFields(GaugeRecipes& recipes,
                                 const ConfigurationEntry& entry,
                                 const Sitmo::Rng& rng,
                                 const RngLayout& rngLayout,
                                 const uint64_t iGaugeConf,
                                 const std::function<bool()>& checkClaim={})
  {
    validateGaugeRecipes(recipes);
    rngLayout.validate();
    
    if(not lxGeomInited)
      CRASH("NICE gauge preparation requires initialized geometry");
    
    if(rngLayout.gaugeOperations.size()!=recipes.size())
      CRASH("NICE RNG layout does not match the registered gauge recipes");
    
    // Clear every slot before the first producer, preventing stale inputs from a previous entry
    clearGaugeFields(recipes);
    
    bool claimed=true;
    for(size_t slot=0;slot<recipes.size() and claimed;slot++)
      {
        auto& recipe=recipes[slot];
        const OperationRngRange range=
          rngLayout.operationRange(rng,rngLayout.gaugeOperations[slot],iGaugeConf);
        recipe.field=
          recipe.produce(entry,recipes,range);
        
        if(not recipe.field)
          CRASH("NICE gauge recipe %zu returned no field",slot);
	
        if(slot==0 and checkClaim)
          claimed=checkClaim();
      }
    
    if(recipes.empty() and checkClaim)
      claimed=checkClaim();
    
    if(not claimed)
      {
        MASTER_PRINTF("NICE: ownership claim changed for gauge configuration '%s' in output '%s' %s; skipping remaining gauge preparation\n",
                      entry.path.c_str(),entry.outputDirectory.c_str(),
                      recipes.empty()?"before any recipe":"after the initial recipe");
        clearGaugeFields(recipes);
      }
    
    return claimed;
  }
  
  /// Expose gauge recipe IDs and deferred builders to the script evaluator
  /// The recipe collection must outlive the registered functions
  void registerGaugeRecipes(Evaluator& evaluator,
                            GaugeRecipes& recipes);
}

#endif
