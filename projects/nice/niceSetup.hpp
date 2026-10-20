#ifndef _NICE_SETUP_HPP
#define _NICE_SETUP_HPP

struct Evaluator;

namespace nissa::nice
{
  struct Correlations;
  
  /// Set up correlation builders and the legacy-input diagnostic function
  /// The specification must outlive the evaluator's registered functions
  void setupNice(Evaluator& evaluator,
                 Correlations& correlations);
}

#endif
