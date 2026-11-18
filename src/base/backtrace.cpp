#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/backtrace.hpp"
#include "geometry/geometry_lx.hpp"
#include "routines/ios.hpp"

#include <cerrno>
#include <atomic>
#include <cstdint>
#include <csignal>
#include <cstdlib>
#include <link.h>
#include <optional>
#include <string>
#include <ucontext.h>
#include <unistd.h>
#include <unwind.h>
#include <vector>

#ifdef __linux__
# include <sys/syscall.h>
#endif

#ifdef USE_BACKTRACE
# include <backtrace.h>
#endif

#ifdef USE_LIBUNWIND
# define UNW_LOCAL_ONLY
# include <libunwind.h>
#endif

namespace nissa
{
  namespace
  {
    /// Executable segment and its relocation bias, cached before installing handlers
    struct TraceObject
    {
      /// First mapped executable address
      uintptr_t start;
      
      /// First address beyond the segment
      uintptr_t end;
      
      /// Relocation bias; zero for a conventional non-PIE executable
      uintptr_t bias;
      
      /// ELF object containing debug information for these addresses
      std::string path;
    };
    
    /// Immutable snapshot read by the signal handler without calling the dynamic loader
    std::vector<TraceObject> traceObjects;
    
    /// Suppress recursive or simultaneous stack unwinding in one process
    std::atomic_flag handlingSignal=ATOMIC_FLAG_INIT;

    /// Enable individual unwinder step diagnostics when requested before initialization
    bool verboseBacktrace{false};
    
#ifdef USE_BACKTRACE
    /// Process-lifetime symbol cache, initialized before signals and shared across threads
    backtrace_state* traceState{nullptr};
#endif
    
    /// One bounded stderr record, formatted without allocation or stdio
    struct TraceRecord
    {
      /// Buffer kept smaller than Linux PIPE_BUF to reduce interleaved MPI output
      char data[4096];
      
      /// Number of initialized bytes
      size_t size{0};
      
      /// Append a character while leaving room for the terminating newline
      void character(const char value)
      {
        if(size<sizeof(data)-1)
          data[size++]=value;
      }
      
      /// Append literal text or an immutable cached object path
      void text(const char* value)
      {
        for(size_t i=0;value[i]!='\0';i++)
          character(value[i]);
      }
      
      /// Write an integer in decimal or hexadecimal without library formatting
      void number(uintptr_t value,
                  const unsigned base=10)
      {
        char digits[3*sizeof(uintptr_t)];
        size_t count=0;
        do
          {
            digits[count++]="0123456789abcdef"[value%base];
            value/=base;
          }
        while(value!=0);
        
        while(count!=0)
          {
            count--;
            character(digits[count]);
          }
      }
      
      /// Label every line by rank, process, and Linux thread without allocating
      void prefix(const char* kind)
      {
        text("NISSA ");
        text(kind);
        text(" rank=");
        number(rank);
        text(" pid=");
        number(getpid());
#ifdef __linux__
        text(" tid=");
        number(syscall(SYS_gettid));
#endif
      }
      
      /// Emit directly to stderr, retrying only an interrupted or partial write
      void emit()
      {
        data[size++]='\n';
        size_t offset=0;
        bool writable=true;
        
        while(offset<size and writable)
          {
            const ssize_t written=write(STDERR_FILENO,data+offset,size-offset);
            if(written>0)
              offset+=static_cast<size_t>(written);
            else if(written==0 or errno!=EINTR)
              writable=false;
          }
      }
    };
    
#ifdef USE_BACKTRACE
    /// Keep the raw frame when debug information is absent or cannot be read
    void ignoreSymbolError(void*,
                           const char*,
                           int)
    {
    }
    
    /// Resolve an initial address outside signal context to populate the library cache
    int warmupSymbols(void*,
                      uintptr_t,
                      const char*,
                      int,
                      const char*)
    {
      return 0;
    }
    
    /// Emit one source location for each physical or inlined function without stdio
    int emitSourceLocation(void* data,
                           uintptr_t,
                           const char* filename,
                           const int line,
                           const char* function)
    {
      if(filename!=nullptr or function!=nullptr)
        {
          TraceRecord record=*static_cast<const TraceRecord*>(data);
          if(function!=nullptr)
            {
              record.text(" function=\"");
              record.text(function);
              record.character('"');
            }
          
          if(filename!=nullptr)
            {
              record.text(" source=\"");
              record.text(filename);
              if(line>0)
                {
                  record.character(':');
                  record.number(line);
                }
              
              record.character('"');
            }
          
          record.emit();
        }
      
      return 0;
    }
#endif
    
    /// Cache executable PT_LOAD mappings using ELF load bias, including for PIE and DSOs
    int collectTraceObject(dl_phdr_info* info,
                           size_t,
                           void*)
    {
      char executable[4096];
      const char* path=info->dlpi_name;
      if(path==nullptr or *path=='\0')
        {
          const ssize_t length=readlink("/proc/self/exe",executable,sizeof(executable)-1);
          if(length>=0)
            {
              executable[length]='\0';
              path=executable;
            }
          else
            path="/proc/self/exe";
        }
      
      for(size_t i=0;i<info->dlpi_phnum;i++)
        {
          const auto& segment=info->dlpi_phdr[i];
          if(segment.p_type==PT_LOAD and (segment.p_flags&PF_X)!=0)
            traceObjects.push_back({info->dlpi_addr+segment.p_vaddr,
                                    info->dlpi_addr+segment.p_vaddr+segment.p_memsz,
                                    info->dlpi_addr,path});
        }
      
      return 0;
    }
    
    /// Keep interrupted instructions exact and move return addresses back into their caller
    uintptr_t lookupProgramCounter(const uintptr_t pc,
                                   const bool interrupted)
    {
      return not interrupted and pc!=0?pc-1:pc;
    }
    
    /// Emit an address without entering either the unwinder or the symbolizer
    void emitTraceFrame(const uintptr_t pc,
                        const bool interrupted,
                        const size_t frame)
    {
      const uintptr_t lookup=lookupProgramCounter(pc,interrupted);
      const TraceObject* object=nullptr;
      for(size_t i=0;i<traceObjects.size() and object==nullptr;i++)
        if(lookup>=traceObjects[i].start and lookup<traceObjects[i].end)
          object=&traceObjects[i];
      
      TraceRecord record;
      record.prefix("frame");
      record.text(interrupted?" kind=interrupted":" kind=return");
      record.text(" index=");
      record.number(frame);
      record.text(" pc=0x");
      record.number(pc,16);
      if(object!=nullptr)
        {
          record.text(" elf=0x");
          record.number(lookup-object->bias,16);
          record.text(" object=\"");
          record.text(object->path.c_str());
          record.character('"');
        }
      else
        record.text(" object=unknown");
      
      record.emit();
    }
    
#ifdef USE_BACKTRACE
    /// Resolve an already printed address after all collected raw frames have been emitted
    void emitTraceLocation(const uintptr_t pc,
                           const bool interrupted,
                           const size_t frame)
    {
      if(traceState!=nullptr)
        {
          TraceRecord location;
          location.prefix("source");
          location.text(interrupted?" kind=interrupted":" kind=return");
          location.text(" index=");
          location.number(frame);
          backtrace_pcinfo(traceState,lookupProgramCounter(pc,interrupted),
                           emitSourceLocation,ignoreSymbolError,&location);
        }
    }
#endif
    
    /// Identify the diagnostic phase so a blocked or failed library call is visible
    void emitTraceStage(const char* phase,
                        const char* state,
                        const int count=-1,
                        const std::optional<int> result={},
                        const char* backend=nullptr)
    {
      TraceRecord record;
      record.prefix(phase);
      record.text(" state=");
      record.text(state);
      if(count>=0)
        {
          record.text(" frames=");
          record.number(count);
        }
      
      if(result)
        {
          record.text(" result=");
          const int64_t value=*result;
          if(value<0)
            record.character('-');
          
          record.number(value<0?-value:value);
        }
      
      if(backend!=nullptr)
        {
          record.text(" backend=");
          record.text(backend);
        }
      
      record.emit();
    }
    
    /// Bounded address storage shared by unwinder callbacks without dynamic allocation
    struct UnwindTrace
    {
      /// Maximum number of collected caller frames
      static constexpr int capacity{128};
      
      /// Collected program counters, filled in order by the callback
      uintptr_t frames[capacity];
      
      /// Number of valid addresses collected so far
      int count{0};
      
      /// Disable diagnostic output during initialization
      bool report{true};
      
      /// Append one nonzero address while storage remains, reporting it immediately
      bool append(const uintptr_t pc)
      {
        const int frame=count;
        frames[count++]=pc;
        if(report)
          {
            if(verboseBacktrace)
              emitTraceStage("unwind","progress",count);
            emitTraceFrame(pc,false,frame);
          }
        
        const bool hasRoom=count<capacity;
        if(not hasRoom and report)
          emitTraceStage("unwind","limit",count);
        
        return hasRoom;
      }
    };
    
    /// Record and print each frame before the unwinder attempts the next step
    _Unwind_Reason_Code collectUnwindFrame(_Unwind_Context* context,
                                          void* data)
    {
      auto& trace=*static_cast<UnwindTrace*>(data);
      const uintptr_t pc=_Unwind_GetIP(context);
      _Unwind_Reason_Code result=_URC_NO_REASON;
      if(pc!=0 and not trace.append(pc))
        result=_URC_END_OF_STACK;
      
      return result;
    }
    
#ifdef USE_LIBUNWIND
    /// Start from the saved interrupted state, avoiding a walk through the signal trampoline
    int collectSignalFrames(UnwindTrace& trace,
                            void* signalContext)
    {
      unw_cursor_t cursor;

      if(verboseBacktrace)
        emitTraceStage("unwind","context-init",trace.count);

      int result=unw_init_local2(&cursor,static_cast<unw_context_t*>(signalContext),
                                 UNW_INIT_SIGNAL_FRAME);

      if(verboseBacktrace or result!=0)
        emitTraceStage("unwind","context-ready",trace.count,result);

      bool walking=result==0;
      while(walking)
        {
          // The initial interrupted PC was printed separately, so advance to its caller
          if(verboseBacktrace)
            emitTraceStage("unwind","step-begin",trace.count);

          result=unw_step(&cursor);

          if(verboseBacktrace or result<0)
            emitTraceStage("unwind","step-complete",trace.count,result);

          if(result>0)
            {
              unw_word_t pc=0;
              result=unw_get_reg(&cursor,UNW_REG_IP,&pc);
              if(result==0 and pc!=0)
                walking=trace.append(pc);
              else
                {
                  emitTraceStage("unwind","register-stop",trace.count,result);
                  walking=false;
                }
            }
          else
            walking=false;
        }
      
      return result;
    }
#endif
    
    /// Emit frames during unwinding so a later blocked step cannot hide earlier addresses
    void emitStackTrace(const uintptr_t interruptedPc=0,
                        void* signalContext=nullptr)
    {
      if(interruptedPc!=0)
        emitTraceFrame(interruptedPc,true,0);
      
      UnwindTrace trace;
      int result;
#ifdef USE_LIBUNWIND
      if(signalContext!=nullptr)
        {
          emitTraceStage("unwind","begin",0,{},"libunwind");
          result=collectSignalFrames(trace,signalContext);
        }
      else
#else
      (void)signalContext;
#endif
        {
          emitTraceStage("unwind","begin",0,{},"abi");
          result=_Unwind_Backtrace(collectUnwindFrame,&trace);
        }
      
      emitTraceStage("unwind","complete",trace.count,result);
      
#ifdef USE_BACKTRACE
      if(traceState!=nullptr)
        {
          if(verboseBacktrace)
            emitTraceStage("symbols","begin");
          
          if(interruptedPc!=0)
            emitTraceLocation(interruptedPc,true,0);
          
          for(int i=0;i<trace.count;i++)
            emitTraceLocation(trace.frames[i],false,i);
          
          if(verboseBacktrace)
            emitTraceStage("symbols","complete");
        }
#endif
    }
    
    /// Get the interrupted instruction pointer before entering the stack unwinder
    uintptr_t interruptedProgramCounter(void* context)
    {
      const auto* state=static_cast<const ucontext_t*>(context);
#if defined(__linux__) and defined(__x86_64__)
      return static_cast<uintptr_t>(state->uc_mcontext.gregs[REG_RIP]);
#elif defined(__linux__) and defined(__aarch64__)
      return static_cast<uintptr_t>(state->uc_mcontext.pc);
#else
      (void)state;
      return 0;
#endif
    }
    
    /// Report external requests on rank zero, but retain faults on whichever rank failed
    /// The interrupted PC is emitted first; subsequent unwinding remains best effort
    void diagnosticSignalHandler(const int signal,
                                  siginfo_t*,
                                  void* context)
    {
      const int savedErrno=errno;
      const bool alreadyUnwinding=handlingSignal.test_and_set(std::memory_order_relaxed);
      
      const bool externalRequest=signal==SIGTERM or signal==SIGINT or signal==SIGUSR1;
      const bool report=not externalRequest or rank==0 or everyRankPrint;
      
      if(report)
        {
          TraceRecord record;
          record.prefix("signal");
          record.text(" number=");
          record.number(signal);
          if(alreadyUnwinding)
            {
              record.text(" state=already-unwinding terminating exit=");
              record.number(128+signal);
            }
          else
            record.text(signal==SIGUSR1?" snapshot":" terminating");
          
          record.emit();
        }
      
      // Expose a second handler before it terminates a stack walk active on this process
      if(alreadyUnwinding)
        _exit(128+signal);
      
      if(report)
        emitStackTrace(interruptedProgramCounter(context),context);
      
      if(signal==SIGUSR1)
        {
          handlingSignal.clear(std::memory_order_relaxed);
          errno=savedErrno;
        }
      else
        _exit(128+signal);
    }
  }
  
  void initializeBacktrace()
  {
    verboseBacktrace=getenv("NISSA_BACKTRACE_VERBOSE")!=nullptr;
    
    // Load the unwinder before a signal can interrupt a memory allocation or loader lock
    UnwindTrace warmup;
    warmup.report=false;
    (void)_Unwind_Backtrace(collectUnwindFrame,&warmup);
    
#ifdef USE_LIBUNWIND
    unw_context_t context;
    unw_cursor_t cursor;
    if(unw_getcontext(&context)==0 and unw_init_local(&cursor,&context)==0)
      {
        int result=1;
        for(int i=0;i<UnwindTrace::capacity and result>0;i++)
          result=unw_step(&cursor);
      }
#endif
    
    traceObjects.clear();
    dl_iterate_phdr(collectTraceObject,nullptr);
    
#ifdef USE_BACKTRACE
    if(traceState==nullptr)
      traceState=backtrace_create_state(nullptr,1,ignoreSymbolError,nullptr);
    
    if(traceState!=nullptr)
      backtrace_pcinfo(traceState,reinterpret_cast<uintptr_t>(&initializeBacktrace),
                       warmupSymbols,ignoreSymbolError,nullptr);
#endif
  }
  
  bool installSignalHandlers()
  {
    if(traceObjects.empty())
      initializeBacktrace();
    
    struct sigaction handler{};
    handler.sa_sigaction=diagnosticSignalHandler;
    handler.sa_flags=SA_SIGINFO|SA_RESTART;
    sigfillset(&handler.sa_mask);
    const int signals[]{SIGBUS,SIGSEGV,SIGFPE,SIGXCPU,SIGABRT,SIGINT,SIGTERM,SIGUSR1};
    bool installed=true;
    for(size_t i=0;i<sizeof(signals)/sizeof(signals[0]) and installed;i++)
      installed=sigaction(signals[i],&handler,nullptr)==0;
    
    return installed;
  }
  
  void printBacktraceList(const int which_rank)
  {
    if(which_rank==-1 or which_rank==rank)
      emitStackTrace();
  }
}
