#include "niceSetup.hpp"

#include "correlations.hpp"
#include "nASTy.hpp"

#include <string>

namespace nissa::nice
{
  void setupNice(Evaluator& evaluator,
                 Correlations& correlations)
  {
    registerHostFunction(evaluator,"corrBlock",{"name"},
                         [&correlations](const std::string& name)
    {
      correlations.addBlock(name);
    });
    
    registerHostFunction(evaluator,"corrGammas",{"name","sink","source"},
                         [&correlations](const std::string& name,
                                         const int sink,
                                         const int source)
    {
      correlations.block(name).addGammas(sink,source);
    });
    
    registerHostFunction(evaluator,"corrTrace",{"name","backward","forward"},
                         [&correlations](const std::string& name,
                                         const Line& backward,
                                         const Line& forward)
    {
      correlations.block(name).addTrace(backward,forward);
    });
  }
}
