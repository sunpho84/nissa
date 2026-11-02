#ifndef _NICE_CONFIGURATION_SOURCES_HPP
#define _NICE_CONFIGURATION_SOURCES_HPP

#include <string>
#include <string_view>
#include <vector>

struct Evaluator;

namespace nissa::nice
{
  /// List paths registered by script calls in evaluation order
  struct ConfigurationSources
  {
    /// Paths interpreted relative to the process working directory
    std::vector<std::string> lists;
  };
  
  /// One configuration and its independent output directory
  struct ConfigurationEntry
  {
    /// Input gauge path, resolved against the list directory
    std::string path;
    
    /// Output directory, resolved against the list directory
    std::string outputDirectory;
    
    bool operator==(const ConfigurationEntry&) const=default;
  };
  
  /// Register a list source without reading it during script evaluation
  void registerConfigurationSources(Evaluator& evaluator,
                                    ConfigurationSources& sources);
  
  /// Parse two-column rows, using a tab to distinguish paths containing spaces
  std::vector<ConfigurationEntry> parseConfigurationList(std::string_view text,
                                                        const std::string& listPath);
  
  /// Read each list on the master, broadcast its bytes, and validate the combined list
  /// MPI and rank information must already be initialized
  std::vector<ConfigurationEntry> materializeAndBroadcastConfigurations(const ConfigurationSources& sources);
}

#endif
