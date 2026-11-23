#ifndef _NICE_SOURCE_GENERATION_HPP
#define _NICE_SOURCE_GENERATION_HPP

#include "corrParts.hpp"
#include "runOptions.hpp"
#include "new_types/su3.hpp"

namespace nissa::nice
{
  /// Prepared source generator and its exact uniform-double consumption
  /// Site is the fastest RNG index, followed by real component, color, spin, and copy
  /// Dilution projects one underlying noise field and does not increase its cost
  struct SourceGeneration
  {
    /// Supported random distributions and deterministic diagnostic fills
    enum class Noise{Z4,Z2,Gauss,AllPlusOne,Point};
    
    /// Distribution selected before entering numerical kernels
    Noise noise;
    
    /// Global spatial volume, also the site count of a selected time slice
    uint64_t spatialVolume;
    
    /// Temporal extent used to validate a selected slice
    int temporalExtent;
    
    /// Selected absolute time, or -1 for the whole lattice
    int tSelect;
    
    /// Independent source copies
    uint64_t nCopies;
    
    /// Project spin onto the requested outer source column
    bool diluteSpin;
    
    /// Project color onto the requested outer source column
    bool diluteColor;
    
    /// Imported lines consume no source-generation draws
    bool imported;
    
    /// Decode a script-visible distribution name without passing strings to kernels
    static Noise parseNoise(const std::string& name)
    {
      if(name=="Z4")
        return Noise::Z4;
      else if(name=="Z2")
        return Noise::Z2;
      else if(name=="Gauss")
        return Noise::Gauss;
      else if(name=="AllPlusOne")
        return Noise::AllPlusOne;
      else if(name=="Point")
        return Noise::Point;
      else
        CRASH("Unknown NICE source noise '%s': expected Z4, Z2, Gauss, AllPlusOne, or Point",name.c_str());
    }
    
    /// Capture geometry and dilution, validating support and total draw arithmetic once
    SourceGeneration(const Source& source,
                     const RunOptions& options) :
      noise(parseNoise(source.noise)),
      spatialVolume(options.physics.volume()/options.physics.T),
      temporalExtent(options.physics.T),
      tSelect(source.tSelect),
      nCopies(options.physics.nCopies),
      diluteSpin(options.physics.diluteSpin),
      diluteColor(options.physics.diluteColor),
      imported(source.ildgPath.has_value())
    {
      options.validate();
      
      if(tSelect<-1 or tSelect>=temporalExtent)
        CRASH("NICE source time %d is outside [-1,%d)",tSelect,temporalExtent);
      
      doubleDrawCost();
    }
    
    /// Number of uniform doubles transformed into each populated complex component
    CUDA_HOST_AND_DEVICE
    uint64_t doublesPerComplex() const
    {
      if(noise==Noise::AllPlusOne or noise==Noise::Point)
        return 0;
      else if(noise==Noise::Z2)
        return 1;
      else
        return 2;
    }
    
    /// Sites on which stochastic or AllPlusOne values are filled
    CUDA_HOST_AND_DEVICE
    uint64_t activeSites() const
    {
      return spatialVolume*(tSelect<0?temporalExtent:1);
    }
    
    /// Count only consumed uniforms, including every copy and underlying spin/color entry
    uint64_t doubleDrawCost() const
    {
      if(imported or doublesPerComplex()==0)
        return 0;
      else
        {
          const uint64_t perCopy=
            checkedMultiply(activeSites(),uint64_t{NDIRAC*NCOL}*doublesPerComplex(),
                            "NICE source draw cost overflow");
          
          return checkedMultiply(perCopy,nCopies,"NICE source copy draw cost overflow");
        }
    }
    
    /// Fill one global site, preserving identical values under MPI decomposition and dilution
    /// The caller supplies the operation's already-positioned random range
    CUDA_HOST_AND_DEVICE
    void fillSite(spincolor& out,
                  const OperationRngRange& range,
                  const uint64_t globalSite,
                  const int time,
                  const uint64_t iCopy,
                  const size_t sourceSpin,
                  const size_t sourceColor) const
    {
      const bool active=tSelect<0 or time==tSelect;
      const uint64_t site=tSelect<0?globalSite:globalSite%spatialVolume;
      
      for(size_t spin=0;spin<NDIRAC;spin++)
        for(size_t color=0;color<NCOL;color++)
          {
            complex value{};
            
            if(active and (not diluteSpin or spin==sourceSpin) and
               (not diluteColor or color==sourceColor))
              {
                if(noise==Noise::AllPlusOne)
                  value[RE]=1;
                else
                  {
                    const uint64_t component=(iCopy*NDIRAC+spin)*NCOL+color;
                    const uint64_t offset=component*doublesPerComplex()*activeSites()+site;
                    value[RE]=range.uniformAt(offset);
                    
                    if(noise==Noise::Z2)
                      z2Transform(value);
                    else
                      {
                        value[IM]=range.uniformAt(offset+activeSites());
                        
                        if(noise==Noise::Z4)
                          z4Transform(value);
                        else
                          {
                            const complex mean{0,0};
                            const complex sigma{M_SQRT1_2,M_SQRT1_2};
                            BoxMullerTransform(value,mean,sigma);
                          }
                      }
                  }
              }
            
            out[spin][color][RE]=value[RE];
            out[spin][color][IM]=value[IM];
          }
    }
  };
}

#endif
