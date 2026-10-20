#ifndef _LEGACY_INPUT_WRITER_HPP
#define _LEGACY_INPUT_WRITER_HPP

#include "compiledCorrelations.hpp"

#include <string>

namespace nissa::nice
{
  /// Render scheduled correlations in the prepareInput diagnostic format
  std::string writeLegacyInput(const CompiledCorrelations& compiled);
}

#endif
