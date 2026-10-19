#ifndef _CORRELATIONS_HPP
#define _CORRELATIONS_HPP

#include "corrParts.hpp"

#include <algorithm>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace nissa::nice
{
  /// Named contraction request and the line pairs contributing to it
  struct Correlation
  {
    /// Name of the requested result block
    std::string name;
    
    /// Sink and source gamma indices requested for the block
    std::set<std::pair<int,int>> gammaList;
    
    /// Backward and forward line pairs in insertion order
    std::vector<std::pair<Line,Line>> traces;
    
    /// Add a gamma pair, retaining each distinct pair once
    void addGammas(const int sink,
                   const int source)
    {
      gammaList.emplace(sink,source);
    }
    
    /// Add a trace pair, rejecting a duplicate structural pair
    void addTrace(const Line& backward,
                  const Line& forward)
    {
      const auto trace=
        std::pair{backward,forward};
      
      if(std::find(traces.begin(),traces.end(),trace)!=traces.end())
        CRASH("duplicate trace in correlation block %s",name.c_str());
      
      traces.push_back(trace);
    }
  };
  
  /// Correlation requests independent of graph lowering and scheduling
  struct Correlations
  {
    /// Requested blocks in creation order
    std::vector<Correlation> blocks;
    
    /// Create a named block, rejecting an empty or reused name
    Correlation& addBlock(const std::string& name)
    {
      if(name.empty())
        CRASH("correlation block name must not be empty");
      
      const auto alreadyPresent=
        std::find_if(blocks.begin(),blocks.end(),[&name](const Correlation& block)
        {
          return block.name==name;
        });
      
      if(alreadyPresent!=blocks.end())
        CRASH("correlation block %s already exists",name.c_str());
      
      blocks.push_back(Correlation{.name=name,.gammaList{},.traces{}});
      
      return blocks.back();
    }
    
    /// Find a previously created block by name
    Correlation& block(const std::string& name)
    {
      const auto found=
        std::find_if(blocks.begin(),blocks.end(),[&name](const Correlation& block)
        {
          return block.name==name;
        });
      
      if(found==blocks.end())
        CRASH("correlation block %s does not exist",name.c_str());
      
      return *found;
    }
  };
}

#endif
