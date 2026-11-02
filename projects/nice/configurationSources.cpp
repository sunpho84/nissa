#include "configurationSources.hpp"

#include "base/debug.hpp"
#include "nASTy.hpp"
#include "routines/mpi_routines.hpp"

#include <filesystem>
#include <fstream>
#include <iterator>
#include <system_error>
#include <unordered_set>

namespace nissa::nice
{
  void registerConfigurationSources(Evaluator& evaluator,
                                    ConfigurationSources& sources)
  {
    registerHostFunction(evaluator,"configurationsFromFile",{"path"},
                         [&sources](const std::string& path)
    {
      if(path.empty())
        CRASH("NICE configurationsFromFile requires a nonempty path");
      
      sources.lists.push_back(path);
    });
  }
  
  std::vector<ConfigurationEntry> parseConfigurationList(const std::string_view text,
                                                        const std::string& listPath)
  {
    namespace fs=std::filesystem;
    std::error_code error;
    const fs::path absolutePath=fs::absolute(fs::path(listPath),error);
    
    if(error)
      CRASH("NICE list '%s': cannot resolve path: %s",listPath.c_str(),error.message().c_str());
    
    const fs::path base=absolutePath.lexically_normal().parent_path();
    std::vector<ConfigurationEntry> entries;
    std::size_t start=0;
    std::size_t lineNumber=1;
    
    while(start<text.size())
      {
        const std::size_t end=text.find('\n',start);
        const std::size_t length=(end==std::string_view::npos?text.size():end)-start;
        std::string_view line=text.substr(start,length);
        
        if(not line.empty() and line.back()=='\r')
          line.remove_suffix(1);
        
        const auto first=line.find_first_not_of(" \t\v\f\r");
        
        if(first!=std::string_view::npos and line[first]!='#')
          {
            const std::size_t separator=line.find('\t');
            
            if(separator==std::string_view::npos or separator==0 or
               separator+1==line.size() or line.find('\t',separator+1)!=std::string_view::npos or
               line.find('\r')!=std::string_view::npos or line.find('\0')!=std::string_view::npos)
              CRASH("NICE list '%s' line %zu: expected one tab-separated configuration path and output directory",
                    listPath.c_str(),lineNumber);
            
            entries.push_back({(base/fs::path(line.substr(0,separator))).lexically_normal().string(),
                               (base/fs::path(line.substr(separator+1))).lexically_normal().string()});
          }
        
        start=end==std::string_view::npos?text.size():end+1;
        lineNumber++;
      }
    
    return entries;
  }
  
  std::vector<ConfigurationEntry> materializeAndBroadcastConfigurations(const ConfigurationSources& sources)
  {
    std::vector<ConfigurationEntry> entries;
    std::unordered_set<std::string> outputs;
    
    for(const std::string& path : sources.lists)
      {
        std::string text;
        std::string absolutePath;
        
        if(rank==master_rank)
          {
            std::error_code error;
            const auto resolved=std::filesystem::absolute(std::filesystem::path(path),error);
            
            if(error)
              CRASH("NICE list '%s': cannot resolve path: %s",path.c_str(),error.message().c_str());
            
            absolutePath=resolved.lexically_normal().string();
            std::ifstream file(absolutePath,std::ios::binary);
            
            if(not file)
              CRASH("NICE list '%s': cannot open file",path.c_str());
            
            text.assign(std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>());
            
            if(file.bad())
              CRASH("NICE list '%s': cannot read file",path.c_str());
          }
        
        broadcast(absolutePath);
        broadcast(text);
        const auto parsed=parseConfigurationList(text,absolutePath);
        
        for(const auto& entry : parsed)
          {
            if(not outputs.insert(entry.outputDirectory).second)
              CRASH("NICE list '%s': duplicate output directory '%s'",
                    path.c_str(),entry.outputDirectory.c_str());
            
            entries.push_back(entry);
          }
      }
    
    return entries;
  }
}
