#ifndef _BACKTRACE_HPP
#define _BACKTRACE_HPP

namespace nissa
{
  /// Print stack addresses and available source locations on selected ranks, or all for -1
  void print_backtrace_list(int which_rank=0);
  
  /// Cache ELF mappings and initialize unwinding and symbolization before trapping signals
  void initializeBacktrace();
  
  /// Install terminating signal diagnostics and nonterminating SIGUSR1 stack snapshots
  bool installSignalHandlers();
}

#endif
