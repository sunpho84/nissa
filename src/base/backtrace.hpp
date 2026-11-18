#ifndef _BACKTRACE_HPP
#define _BACKTRACE_HPP

namespace nissa
{
  /// Print stack addresses and available source locations on selected ranks, or all for -1
  void printBacktraceList(int which_rank=0);
  
  /// Cache ELF mappings and initialize unwinding and symbolization before trapping signals
  void initializeBacktrace();
  
  /// Install signal diagnostics: external requests report on rank zero unless EVERY_RANK_PRINT is set
  /// Fault signals report on the failing rank; SIGUSR1 snapshots resume execution
  bool installSignalHandlers();
}

#endif
