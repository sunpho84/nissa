#include "configurationSources.hpp"

#include "base/debug.hpp"
#include "nASTy.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"

#include <filesystem>
#include <sstream>
#include <system_error>
#include <unordered_map>

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
    
    std::vector<ConfigurationEntry> entries;
    size_t start=0;
    size_t rowNumber=1;
    
    while(start<text.size())
      {
        const size_t end=text.find('\n',start);
        const size_t length=(end==std::string_view::npos?text.size():end)-start;
        std::string_view row=text.substr(start,length);
        
        if(not row.empty() and row.back()=='\r')
          row.remove_suffix(1);
        
        const auto first=row.find_first_not_of(" \t\v\f\r");
        
        if(first!=std::string_view::npos and row[first]!='#')
          {
            const size_t separator=row.find('\t');
            std::string input;
            std::string output;
            bool malformed=row.find('\r')!=std::string_view::npos or
              row.find('\0')!=std::string_view::npos;
            
            if(separator!=std::string_view::npos)
              {
                // A single tab separates fields that may themselves contain spaces
                malformed=malformed or row.find('\t',separator+1)!=std::string_view::npos;
                input=row.substr(0,separator);
                output=row.substr(separator+1);
                malformed=malformed or input.find_first_not_of(" \v\f")==std::string::npos or
                  output.find_first_not_of(" \v\f")==std::string::npos;
              }
            else
              {
                // Without a tab, exactly two whitespace-delimited paths are allowed
                std::istringstream columns{std::string(row)};
                std::string extra;
                malformed=malformed or not (columns>>input>>output) or (columns>>extra);
              }
            
            if(malformed)
              CRASH("NICE list '%s' row %zu: expected two columns; use one tab when a path contains spaces",
                    listPath.c_str(),rowNumber);
            
            const fs::path inputPath=fs::absolute(fs::path(input),error);

            if(error)
              CRASH("NICE list '%s' row %zu: cannot resolve configuration path: %s",
                    listPath.c_str(),rowNumber,error.message().c_str());

            const fs::path outputPath=fs::absolute(fs::path(output),error);

            if(error)
              CRASH("NICE list '%s' row %zu: cannot resolve output directory: %s",
                    listPath.c_str(),rowNumber,error.message().c_str());

            entries.push_back({inputPath.lexically_normal().string(),
                               outputPath.lexically_normal().string(),
                               absolutePath.lexically_normal().string(),rowNumber});
          }
        
        start=end==std::string_view::npos?text.size():end+1;
        rowNumber++;
      }
    
    return entries;
  }
  
  std::vector<ConfigurationEntry> materializeAndBroadcastConfigurations(const ConfigurationSources& sources)
  {
    std::vector<ConfigurationEntry> entries;
    std::unordered_map<std::string,size_t> outputs;
    
    for(const std::string& path : sources.lists)
      {
        std::string absolutePath;
        
        if(rank==master_rank)
          {
            std::error_code error;
            const auto resolved=std::filesystem::absolute(std::filesystem::path(path),error);
            
            if(error)
              CRASH("NICE list '%s': cannot resolve path: %s",path.c_str(),error.message().c_str());
            
            absolutePath=resolved.lexically_normal().string();
          }
        
        broadcast(absolutePath);
        const std::string text=readFileBytes(absolutePath,FileReadMode::Collective);
        const auto parsed=parseConfigurationList(text,absolutePath);
        
        for(const auto& entry : parsed)
          {
            const auto [previous,inserted]=
	      outputs.emplace(entry.outputDirectory,entries.size());
            
            if(not inserted)
              {
                const auto& earlier=entries[previous->second];
                CRASH("NICE duplicate output directory '%s': '%s' row %zu and '%s' row %zu",
                      entry.outputDirectory.c_str(),earlier.sourceList.c_str(),earlier.sourceRow,
                      entry.sourceList.c_str(),entry.sourceRow);
              }
            
            entries.push_back(entry);
          }
      }
    
    return entries;
  }
}
