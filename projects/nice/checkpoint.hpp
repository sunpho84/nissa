#ifndef _NICE_CHECKPOINT_HPP
#define _NICE_CHECKPOINT_HPP

#include "nice.hpp"
#include "checkpointIdentity.hpp"
#include "programExecution.hpp"
#include "rngLayout.hpp"
#include "wallTime.hpp"

#include "io/input.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "routines/fnv1a64.hpp"

#include <atomic>
#include <bit>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <iterator>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <variant>

namespace nissa::nice
{
  static_assert(sizeof(double)==sizeof(uint64_t));
  static_assert(std::numeric_limits<double>::is_iec559);
  
  /// Append one portable little-endian integer to a binary checkpoint
  inline void appendCheckpointUint64(std::string& bytes,
                                     const uint64_t value)
  {
    for(int iByte=0;iByte<8;iByte++)
      bytes.push_back(static_cast<char>(value>>(8*iByte)));
  }
  
  /// Decode one little-endian integer, rejecting truncated checkpoint data
  inline uint64_t readCheckpointUint64(const std::string_view bytes,
                                       size_t& offset)
  {
    if(offset>bytes.size() or bytes.size()-offset<8)
      CRASH("NICE binary checkpoint is truncated");
    
    uint64_t value=0;
    for(int iByte=0;iByte<8;iByte++)
      value|=static_cast<uint64_t>(static_cast<unsigned char>(bytes[offset++]))<<(8*iByte);
    
    return value;
  }
  
  /// Store the complex values of one hit in compiled output order
  inline std::string encodeHitValues(const RunOptions& options,
                                     const CompiledCorrelations& correlations,
                                     const TaskValues& values)
  {
    std::string bytes;
    
    for(const OutputBlock& output : correlations.outputs)
      for(const TaskId traceId : output.traces)
        {
          if(traceId>=values.size() or not values[traceId].has_value())
            CRASH("NICE output trace task %zu has no value",traceId);
          
          if(const Tensor* tensor=std::get_if<Tensor>(&*values[traceId]);tensor!=nullptr)
            for(size_t iGammaPair=0;iGammaPair<output.gammaList.size();iGammaPair++)
              {
                for(size_t iCopy=0;iCopy<static_cast<size_t>(options.physics.nCopies);iCopy++)
                  for(int deltaTime=0;deltaTime<glbSize[0];deltaTime++)
                    {
                      const size_t index=
                        (iGammaPair*static_cast<size_t>(glbSize[0])+deltaTime)*
                        static_cast<size_t>(options.physics.nCopies)+iCopy;
                      appendCheckpointUint64(bytes,std::bit_cast<uint64_t>(tensor->data[index*2]));
                      appendCheckpointUint64(bytes,std::bit_cast<uint64_t>(tensor->data[index*2+1]));
                    }
              }
          else
            CRASH("NICE output trace task %zu is not a tensor",traceId);
        }
    
    return bytes;
  }
  
  /// Format binary hit values as scientific rows for the final result
  inline std::string formatHitResults(const RunOptions& options,
                                      const CompiledCorrelations& correlations,
                                      const std::string_view bytes,
                                      const size_t iHit)
  {
    std::string rows;
    size_t offset=0;
    const int precision=std::numeric_limits<double>::max_digits10;
    
    for(const OutputBlock& output : correlations.outputs)
      for(const TaskId traceId : output.traces)
        for(const auto& [sinkGamma,sourceGamma] : output.gammaList)
          for(size_t iCopy=0;iCopy<static_cast<size_t>(options.physics.nCopies);iCopy++)
            for(int deltaTime=0;deltaTime<glbSize[0];deltaTime++)
              {
                const double real=std::bit_cast<double>(readCheckpointUint64(bytes,offset));
                const double imaginary=std::bit_cast<double>(readCheckpointUint64(bytes,offset));
                std::format_to(std::back_inserter(rows),
                               "{} {} {} {} {} {} {} {:.{}g} {:.{}g}\n",
                               iHit,output.name,traceId,sinkGamma,sourceGamma,
                               iCopy,deltaTime,real,precision,imaginary,precision);
              }
    
    if(offset!=bytes.size())
      CRASH("NICE binary hit checkpoint has extra values");
    
    return rows;
  }
  
  /// Read and share a complete checkpoint without changing its bytes
  inline std::string readCheckpointFile(const std::filesystem::path& path)
  {
    return readFileBytes(path.string(),FileReadMode::Collective);
  }
  
  /// Write all bytes to an already opened checkpoint or result file
  inline void writeCheckpointBytes(FILE* file,
                                   const std::filesystem::path& path,
                                   const std::string_view bytes)
  {
    if(file!=nullptr and not bytes.empty() and
       std::fwrite(bytes.data(),1,bytes.size(),file)!=bytes.size())
      CRASH("NICE cannot write '%s'",path.string().c_str());
  }
  
  /// Check the ownership tag collectively before modifying checkpoint output
  inline void requireCheckpointOwnership(lock_file_t<>& lock)
  {
    if(not lock.check_lock())
      CRASH("NICE ownership tag '%s' changed or could not be read; stopping before checkpoint write",
            lock.lockPath().c_str());
  }
  
  /// Test whether another run has recently updated its running file
  inline bool runningFileIsRecent(const std::filesystem::path& path,
                                  const std::chrono::seconds interval)
  {
    bool recent=false;
    std::error_code error;
    
    if(is_master_rank())
      {
        const std::filesystem::file_status status=
          std::filesystem::symlink_status(path,error);
        
        if(error==std::errc::no_such_file_or_directory)
          error.clear();
        else if(not error and std::filesystem::exists(status))
          {
            if(not std::filesystem::is_regular_file(status))
              error=std::make_error_code(std::errc::invalid_argument);
            else
              {
                const auto updated=std::filesystem::last_write_time(path,error);
                
                if(error==std::errc::no_such_file_or_directory)
                  error.clear();
                else if(not error)
                  recent=std::filesystem::file_time_type::clock::now()-updated<
                    2*interval;
              }
          }
      }
    
    checkCollectiveFileError(error,
                             "NICE cannot inspect running file '"+path.string()+"'");
    
    return broadcast(recent);
  }
  
  /// Write an ownership tag when the previous running file is absent or stale
  inline bool tryClaimRunningFile(lock_file_t<>& lock,
                                  const std::filesystem::path& path,
                                  const std::chrono::seconds interval)
  {
    bool claimed=false;
    
    if(not runningFileIsRecent(path,interval))
      {
        if(not lock.try_lock(path.string()))
          CRASH("NICE cannot write running file '%s'",path.string().c_str());
        
        claimed=true;
      }
    
    return claimed;
  }
  
  /// Refresh the running file while numerical work is in progress
  struct RunningFileUpdater
  {
    /// Ownership tag whose file is refreshed
    lock_file_t<>& lock;
    
    /// Configured refresh period
    std::chrono::milliseconds interval;
    
    /// Background updater on the master rank
    std::thread worker{};
    
    /// Protect the stop request and wake the updater promptly
    std::mutex mutex{};
    
    /// Wake the updater when the configuration finishes
    std::condition_variable wake{};
    
    /// True once the updater is asked to finish
    bool stopping{false};
    
    /// True when the running file could not be refreshed
    std::atomic<bool> updateFailed{false};
    
    /// Start the updater immediately after writing the ownership tag
    explicit RunningFileUpdater(lock_file_t<>& lock,
                                const std::chrono::milliseconds interval) :
      lock(lock),
      interval(interval)
    {
      if(is_master_rank())
        worker=std::thread([this]()
          {
            std::unique_lock guard(mutex);
            
            while(not stopping)
              {
                const bool requestedStop=
                  wake.wait_for(guard,this->interval,
                                [this]()
                                {
                                  return stopping;
                                });
                
                if(not requestedStop)
                  {
                    guard.unlock();
                    
                    const bool owned=this->lock.check_lock_on_master();
                    
                    if(owned)
                      {
                        std::error_code error;
                        std::filesystem::last_write_time(
                          this->lock.lockPath(),
                          std::filesystem::file_time_type::clock::now(),error);
                        
                        if(error)
                          {
                            updateFailed=true;
                            std::fprintf(stderr,
                                         "NICE: cannot refresh running file '%s': %s\n",
                                         this->lock.lockPath().c_str(),error.message().c_str());
                          }
                      }
                    else
                      std::fprintf(stderr,
                                   "NICE: running file '%s' changed ownership; stopping updates\n",
                                   this->lock.lockPath().c_str());
                    
                    guard.lock();
                    
                    if(not owned or updateFailed.load())
                      stopping=true;
                  }
              }
          });
    }
    
    /// Stop and join before the running file is removed
    void stop()
    {
      {
        std::lock_guard guard(mutex);
        stopping=true;
      }
      
      wake.notify_one();
      
      if(worker.joinable())
        worker.join();
    }
    
    /// Report an updater failure collectively before writing checkpoints
    void requireHealthy() const
    {
      if(broadcast(updateFailed.load()))
        CRASH("NICE running file '%s' could not be refreshed",lock.lockPath().c_str());
    }
    
    /// Stop the updater on every exit path
    ~RunningFileUpdater()
    {
      stop();
    }
  };
  
  /// Remove the running file only while its ownership tag still belongs to this run
  inline void releaseRunningFile(lock_file_t<>& lock)
  {
    requireCheckpointOwnership(lock);
    removeFile(lock.lockPath());
  }
  
  /// Write and rename on the same filesystem so interrupted writes stay invisible
  inline void commitCheckpointFile(const std::filesystem::path& path,
                                   const std::string& bytes,
                                   lock_file_t<>& lock)
  {
    const std::filesystem::path temporary=path.string()+".tmp";
    FILE* file=open_file(temporary.string(),"wb");
    writeCheckpointBytes(file,temporary,bytes);
    
    if(file!=nullptr)
      {
        // close_file does not report buffered write failures
        if(std::fclose(file)!=0)
          CRASH("NICE cannot write temporary checkpoint '%s'",temporary.string().c_str());
      }
    
    requireCheckpointOwnership(lock);
    mv(path.string(),temporary.string());
  }
  
  /// Construct a run-specific filename while preserving the file extension
  inline std::filesystem::path outputFilePath(const std::filesystem::path& directory,
                                               const RunOptions& options,
                                               const std::string_view stem,
                                               const std::string_view extension)
  {
    return directory/std::format("{}{}{}",stem,options.policy.outputSuffix,extension);
  }
  
  /// One immutable binary hit file, named by its stable zero-based hit index
  inline std::filesystem::path hitCheckpointPath(const std::filesystem::path& directory,
                                                 const RunOptions& options,
                                                 const size_t iHit)
  {
    return outputFilePath(directory,options,std::format("hit-{:08}",iHit),".bin");
  }
  
  /// Number of scientific rows required from one complete hit
  inline size_t expectedHitRows(const PreparedProgram& program)
  {
    size_t rows=0;
    
    for(const OutputBlock& output : program.correlations.outputs)
      rows+=output.traces.size()*output.gammaList.size()*
        static_cast<size_t>(program.options.physics.nCopies)*static_cast<size_t>(program.options.physics.T);
    
    return rows;
  }
  
  /// Construct the full compatibility description used by checkpoint and finish records
  inline std::string checkpointMetadata(const PreparedProgram& program,
                                        const RngLayout& rngLayout,
                                        const size_t iEntry);
  
  /// Hash the full compatibility description for compact binary records
  inline uint64_t checkpointMetadataHash(const PreparedProgram& program,
                                         const RngLayout& rngLayout,
                                         const size_t iEntry)
  {
    return fnv1a64(checkpointMetadata(program,rngLayout,iEntry));
  }
  
  /// Wrap one hit's binary values with a version, identity, count, and checksum
  inline std::string hitCheckpointBytes(const uint64_t identity,
                                        const size_t iHit,
                                        const std::string_view values)
  {
    if(values.size()%16!=0)
      CRASH("NICE hit values have incomplete complex data");
    
    std::string bytes("NICEHT3\0",8);
    appendCheckpointUint64(bytes,identity);
    appendCheckpointUint64(bytes,iHit);
    appendCheckpointUint64(bytes,values.size()/16);
    appendCheckpointUint64(bytes,fnv1a64(values));
    bytes+=values;
    
    return bytes;
  }
  
  /// Validate a binary hit checkpoint and return its complex-value payload
  inline std::string checkedHitValues(const std::filesystem::path& path,
                                      const uint64_t identity,
                                      const size_t iHit,
                                      const size_t requiredRows)
  {
    const std::string bytes=readCheckpointFile(path);
    if(bytes.size()<40 or bytes.compare(0,8,std::string("NICEHT3\0",8))!=0)
      CRASH("NICE checkpoint '%s' has an invalid binary header",path.string().c_str());
    
    size_t offset=8;
    const uint64_t storedIdentity=readCheckpointUint64(bytes,offset);
    const uint64_t storedHit=readCheckpointUint64(bytes,offset);
    const uint64_t storedRows=readCheckpointUint64(bytes,offset);
    const uint64_t storedHash=readCheckpointUint64(bytes,offset);
    if(storedIdentity!=identity or storedHit!=iHit or storedRows!=requiredRows)
      CRASH("NICE checkpoint '%s' does not match this run",path.string().c_str());
    
    const std::string values=bytes.substr(offset);
    if(requiredRows>std::numeric_limits<size_t>::max()/16 or
       values.size()!=requiredRows*16 or fnv1a64(values)!=storedHash)
      CRASH("NICE checkpoint '%s' is damaged",path.string().c_str());
    
    return values;
  }
  
  /// Hash an output file on master in bounded memory and share its checksum
  inline uint64_t checkpointFileHash(const std::filesystem::path& path)
  {
    uint64_t hash=fnv1a64Offset;
    if(is_master_rank())
      {
        std::ifstream file(path,std::ios::binary);
        if(not file)
          CRASH("NICE cannot open result '%s'",path.string().c_str());
        
        char buffer[65536];
        while(file.read(buffer,sizeof(buffer)) or file.gcount()!=0)
          hash=fnv1a64Update(hash,{buffer,static_cast<size_t>(file.gcount())});
        
        if(not file.eof() or file.bad())
          CRASH("NICE cannot read result '%s'",path.string().c_str());
      }
    
    return broadcast(hash);
  }
  
  /// Verify hits on all ranks and commit final output and marker on master
  inline void finalizeConfiguration(const PreparedProgram& program,
                                    const size_t iEntry,
                                    const RngLayout& rngLayout,
                                    lock_file_t<>& lock)
  {
    requireCheckpointOwnership(lock);
    
    const ConfigurationEntry& entry=program.configurationEntries[iEntry];
    const std::filesystem::path directory(entry.outputDirectory);
    const std::filesystem::path output=outputFilePath(directory,program.options,"correlations",".nice");
    const std::filesystem::path temporary=output.string()+".tmp";
    FILE* file=open_file(temporary.string(),"wb");
    
    writeCheckpointBytes(file,temporary,
                         std::format("# gauge configuration {}\n",entry.path));
    writeCheckpointBytes(file,temporary,
                         "# hit block trace sinkGamma sourceGamma copy deltaT real imaginary\n");
    
    for(size_t iHit=0;iHit<static_cast<size_t>(program.options.policy.nHits);iHit++)
      writeCheckpointBytes(file,temporary,
                           formatHitResults(program.options,program.correlations,
                                            checkedHitValues(hitCheckpointPath(directory,program.options,iHit),
                                                             checkpointEntryIdentity(program,iEntry),iHit,
                                                             expectedHitRows(program)),iHit));
    
    if(file!=nullptr)
      {
        if(std::fclose(file)!=0)
          CRASH("NICE cannot write result '%s'",temporary.string().c_str());
      }
    
    requireCheckpointOwnership(lock);
    mv(output.string(),temporary.string());
    
    const uint64_t outputHash=checkpointFileHash(output);
    // The marker commits both the result and full run metadata before cleanup
    requireCheckpointOwnership(lock);
    commitCheckpointFile(outputFilePath(directory,program.options,"finished",".nice"),
                         std::format("NICE-finished-v3 {} {} {} {}\n",
                                     checkpointEntryIdentity(program,iEntry),
                                     program.options.policy.nHits,
                                     checkpointMetadataHash(program,rngLayout,iEntry),
                                     outputHash),lock);
    
    if(not program.options.policy.keepCheckpoints)
      {
        for(size_t iHit=0;iHit<static_cast<size_t>(program.options.policy.nHits);iHit++)
          {
            requireCheckpointOwnership(lock);
            removeFile(hitCheckpointPath(directory,program.options,iHit).string());
          }
        
        requireCheckpointOwnership(lock);
        removeFile(outputFilePath(directory,program.options,"checkpoint",".nice").string());
        
        if(program.options.policy.outputSuffix.empty() and fileExists((directory/"checkpoint.bin").string()))
          {
            requireCheckpointOwnership(lock);
            removeFile((directory/"checkpoint.bin").string());
          }
      }
    
    MASTER_PRINTF("NICE: wrote results to '%s'\n",output.string().c_str());
  }
  
  /// Check a finished marker and its result before skipping numerical work
  inline bool configurationFinished(const PreparedProgram& program,
                                    const RngLayout& rngLayout,
                                    const size_t iEntry)
  {
    const std::filesystem::path directory=
      program.configurationEntries[iEntry].outputDirectory;
    const std::filesystem::path marker=outputFilePath(directory,program.options,"finished",".nice");
    bool finished=false;
    
    if(fileExists(marker.string()))
      {
        std::istringstream contents(readCheckpointFile(marker));
        std::string version;
        std::string extra;
        uint64_t identity{};
        int nHits{};
        uint64_t metadataHash{};
        uint64_t outputHash{};
        if(not (contents>>version>>identity>>nHits>>metadataHash>>outputHash) or
           contents>>extra or version!="NICE-finished-v3" or
           identity!=checkpointEntryIdentity(program,iEntry) or
           nHits<=0 or nHits>program.options.policy.nHits or
           metadataHash!=checkpointMetadataHash(program,rngLayout,iEntry))
          CRASH("NICE finished marker '%s' does not match this run",marker.string().c_str());
        
        const std::filesystem::path output=outputFilePath(directory,program.options,"correlations",".nice");
        if(checkpointFileHash(output)!=outputHash)
          CRASH("NICE finished result '%s' is damaged",output.string().c_str());
        
        if(nHits<program.options.policy.nHits and
           not fileExists(outputFilePath(directory,program.options,"checkpoint",".nice").string()) and
           not (program.options.policy.outputSuffix.empty() and fileExists((directory/"checkpoint.bin").string())))
          CRASH("NICE cannot increase hits in '%s': completed hit checkpoints were removed",
                directory.string().c_str());
        
        finished=nHits==program.options.policy.nHits;
      }
    
    return finished;
  }

  /// Check a script-selected stop request collectively at a work boundary
  inline bool stopRequested(const RunOptions& options)
  {
    return not options.policy.stopFile.empty() and fileExists(options.policy.stopFile);
  }
  
  /// Record the settings that determine checkpoint compatibility and RNG positions
  inline std::string checkpointMetadata(const PreparedProgram& program,
                                 const RngLayout& rngLayout,
                                 const size_t iEntry)
  {
    const ConfigurationEntry& entry=program.configurationEntries[iEntry];
    const RunOptions& options=program.options;
    std::string metadata=
      std::format("NICE-checkpoint-v3\nidentity {}\nentry {}\ninput {}\noutput {}\n"
                  "geometry {} {}\nseed {}\ncopies {}\n"
                  "dilution {} {}\nboundary {}\norigin {} {} {} {} {}\n"
                  "reservation {}\ntasks {}\nrows-per-hit {}\norigin-cost {}\n",
                  checkpointEntryIdentity(program,iEntry),iEntry,entry.path,entry.outputDirectory,
                  options.physics.L,options.physics.T,options.physics.seed,
                  options.physics.nCopies,options.physics.diluteSpin,options.physics.diluteColor,
                  static_cast<int>(options.physics.temporalBoundary),
                  options.physics.drawGlobalOriginTask.enabled,
                  options.physics.drawGlobalOriginTask.fixedOrigin[0],
                  options.physics.drawGlobalOriginTask.fixedOrigin[1],
                  options.physics.drawGlobalOriginTask.fixedOrigin[2],
                  options.physics.drawGlobalOriginTask.fixedOrigin[3],
                  options.physics.reservedRngDoublesPerSite,
                  program.correlations.instructions.size(),expectedHitRows(program),
                  rngLayout.origin.doubleDrawCost);
    
    for(const RngAllocation& allocation : rngLayout.gaugeOperations)
      metadata+=std::format("gauge-cost {}\n",allocation.doubleDrawCost);
    
    return metadata;
  }
  
  /// Encode legacy version 3 binary metadata for reading existing checkpoints
  inline std::string checkpointMetadataBytes(const PreparedProgram& program,
                                             const RngLayout& rngLayout,
                                             const size_t iEntry)
  {
    std::string bytes("NICECP3\0",8);
    appendCheckpointUint64(bytes,checkpointEntryIdentity(program,iEntry));
    appendCheckpointUint64(bytes,checkpointMetadataHash(program,rngLayout,iEntry));
    
    return bytes;
  }
  
  /// Consecutive attempts at one unfinished hit, independent of compatibility metadata
  struct CheckpointRetryState
  {
    /// Hit whose computation was attempted; older committed progress resets the count
    uint64_t hit{0};
    
    /// Attempts started without successfully committing this hit
    uint64_t nTrials{0};
  };
  
  /// Accept old metadata without retry state and strictly validate new mutable records
  inline CheckpointRetryState parseCheckpointRetryState(const std::string& contents,
                                                        const std::string& expected)
  {
    if(not contents.starts_with(expected))
      CRASH("NICE checkpoint metadata does not match this run");
    
    CheckpointRetryState state;
    if(contents.size()!=expected.size())
      {
        std::istringstream input(contents.substr(expected.size()));
        std::string version;
        std::string hit;
        std::string trials;
        std::string extra;
        
        if(not (input>>version>>hit>>trials) or input>>extra or version!="retry-state-v1")
          CRASH("NICE checkpoint retry state is damaged");
        
        const auto parse=[](const std::string& text,
                            uint64_t& result)
        {
          const auto parsed=std::from_chars(text.data(),text.data()+text.size(),result);
          
          if(parsed.ec!=std::errc{} or parsed.ptr!=text.data()+text.size())
            CRASH("NICE invalid checkpoint retry counter '%s'",text.c_str());
        };
        parse(hit,state.hit);
        parse(trials,state.nTrials);
      }
    
    return state;
  }
  
  /// Persist mutable retry state atomically while retaining unchanged compatibility bytes
  inline void writeCheckpointRetryState(const PreparedProgram& program,
                                        const RngLayout& rngLayout,
                                        const size_t iEntry,
                                        const CheckpointRetryState& state,
                                        lock_file_t<>& lock)
  {
    const auto path=outputFilePath(program.configurationEntries[iEntry].outputDirectory,
                                   program.options,"checkpoint",".nice");
    commitCheckpointFile(path,checkpointMetadata(program,rngLayout,iEntry)+
      std::format("retry-state-v1 {} {}\n",state.hit,state.nTrials),lock);
  }
  
  /// Restore attempts, resetting stale counts when a hit was committed before a crash
  inline CheckpointRetryState readCheckpointRetryState(const PreparedProgram& program,
                                                       const RngLayout& rngLayout,
                                                       const size_t iEntry,
                                                       const size_t firstHit)
  {
    CheckpointRetryState state;
    const auto path=outputFilePath(program.configurationEntries[iEntry].outputDirectory,
                                   program.options,"checkpoint",".nice");
    if(fileExists(path.string()))
      state=parseCheckpointRetryState(readCheckpointFile(path),
                                      checkpointMetadata(program,rngLayout,iEntry));
    
    if(state.hit>firstHit)
      CRASH("NICE checkpoint retry state refers to missing completed hits");
    else if(state.hit<firstHit)
      state={firstHit,0};
    
    return state;
  }
  
  /// Refuse incompatible resume state, creating metadata only with an ownership claim
  inline void validateCheckpointMetadata(const PreparedProgram& program,
					 const RngLayout& rngLayout,
					 const size_t iEntry,
					 lock_file_t<>* lock=nullptr)
  {
    const std::filesystem::path directory=
      program.configurationEntries[iEntry].outputDirectory;
    const std::filesystem::path path=outputFilePath(directory,program.options,"checkpoint",".nice");
    const std::string expected=checkpointMetadata(program,rngLayout,iEntry);
    if(fileExists(path.string()))
      {
        (void)parseCheckpointRetryState(readCheckpointFile(path),expected);
      }
    else if(program.options.policy.outputSuffix.empty() and fileExists((directory/"checkpoint.bin").string()))
      {
        if(readCheckpointFile(directory/"checkpoint.bin")!=checkpointMetadataBytes(program,rngLayout,iEntry))
          CRASH("NICE legacy checkpoint metadata does not match this run");
      }
    else
      {
        if(fileExists(hitCheckpointPath(directory,program.options,0).string()) or
           fileExists(outputFilePath(directory,program.options,"finished",".nice").string()) or
           fileExists(outputFilePath(directory,program.options,"correlations",".nice").string()))
          CRASH("NICE output directory '%s' has results without checkpoint metadata",
                directory.string().c_str());
        
        if(lock==nullptr)
          CRASH("NICE cannot create checkpoint metadata without ownership");
	
        commitCheckpointFile(path,expected,*lock);
      }
  }
  
  /// Find the first unfinished hit and reject gaps or damaged committed hits
  inline size_t firstUnfinishedHit(const PreparedProgram& program,
				   const size_t iEntry)
  {
    const std::filesystem::path directory=
      program.configurationEntries[iEntry].outputDirectory;
    size_t first=static_cast<size_t>(program.options.policy.nHits);
    
    for(size_t iHit=0;iHit<static_cast<size_t>(program.options.policy.nHits);iHit++)
      {
        const std::filesystem::path path=hitCheckpointPath(directory,program.options,iHit);
        if(fileExists(path.string()))
          {
            if(first!=static_cast<size_t>(program.options.policy.nHits))
              CRASH("NICE checkpoint '%s' follows an unfinished hit",path.string().c_str());
            
            checkedHitValues(path,checkpointEntryIdentity(program,iEntry),iHit,
                             expectedHitRows(program));
          }
        else if(first==static_cast<size_t>(program.options.policy.nHits))
          first=iHit;
      }
    
    return first;
  }
  
  /// Result of attempting one configuration's checkpointed work
  enum class ConfigurationRunOutcome
  {
    Completed,
    AlreadyFinished,
    ClaimedElsewhere,
    LostOwnership,
    TrialLimitReached,
    Stopped
  };
  
  /// Resume one configuration collectively, reporting skips and clean stops
  inline ConfigurationRunOutcome runCheckpointedConfiguration(PreparedProgram& program,
                                                                const RngLayout& rngLayout,
                                                                const Sitmo::Rng& rng,
                                                                const size_t iEntry,
                                                                WallTimeBudget& wallTime)
  {
    bool stopped=false;
    ConfigurationRunOutcome outcome=ConfigurationRunOutcome::Completed;
    const std::filesystem::path directory=
      program.configurationEntries[iEntry].outputDirectory;
    if(create_directories(directory.string())!=0)
      CRASH("NICE cannot create output directory '%s'",directory.string().c_str());
    
    if(configurationFinished(program,rngLayout,iEntry))
      {
        MASTER_PRINTF("NICE: skipping completed gauge configuration %zu\n",iEntry);
        outcome=ConfigurationRunOutcome::AlreadyFinished;
      }
    else
      {
        const std::filesystem::path runningPath=
          outputFilePath(directory,program.options,".nice",".lock");
        const std::chrono::seconds updateInterval(program.options.policy.runningFileUpdateSeconds);
        
        lock_file_t<> lock;
        lock.init();
        if(not tryClaimRunningFile(lock,runningPath,updateInterval))
          {
            MASTER_PRINTF("NICE: skipping gauge configuration %zu with recent running file '%s'\n",
                          iEntry,runningPath.string().c_str());
            outcome=ConfigurationRunOutcome::ClaimedElsewhere;
          }
        else
          {
            RunningFileUpdater updater(lock,updateInterval);
            
            // Initial gauge work makes competing running-file writes visible
            const bool claimed=
              prepareConfiguration(program,iEntry,rng,rngLayout,iEntry,
                                   [&lock]()
                                   {
                                     return lock.check_lock();
                                   });
            updater.requireHealthy();
            
            if(not claimed)
              {
                MASTER_PRINTF("NICE: skipping gauge configuration %zu after lost ownership of '%s'\n",
                              iEntry,lock.lockPath().c_str());
                outcome=ConfigurationRunOutcome::LostOwnership;
              }
            else if(configurationFinished(program,rngLayout,iEntry))
              {
                MASTER_PRINTF("NICE: skipping completed gauge configuration %zu\n",iEntry+1);
                outcome=ConfigurationRunOutcome::AlreadyFinished;
                updater.stop();
                releaseRunningFile(lock);
              }
            else
              {
                validateCheckpointMetadata(program,rngLayout,iEntry,&lock);
                const size_t firstHit=firstUnfinishedHit(program,iEntry);
                CheckpointRetryState retry=
                  readCheckpointRetryState(program,rngLayout,iEntry,firstHit);
                const bool trialLimitReached=
                  firstHit<static_cast<size_t>(program.options.policy.nHits) and
                  program.options.policy.maxTrials!=0 and
                  retry.nTrials>=static_cast<uint64_t>(program.options.policy.maxTrials);
                
                if(firstHit==0 and fileExists(outputFilePath(directory,program.options,"correlations",".nice").string()))
                  CRASH("NICE output '%s' exists without checkpoints; refusing to overwrite it",
                        outputFilePath(directory,program.options,"correlations",".nice").string().c_str());
                
                if(firstHit<static_cast<size_t>(program.options.policy.nHits) and not trialLimitReached)
                  {
                    for(size_t iHit=firstHit;
                        iHit<static_cast<size_t>(program.options.policy.nHits) and not stopped;iHit++)
                      {
                        stopped=stopRequested(program.options) or not wallTime.allowsWork();
                        
                        if(not stopped)
                          {
                            // Persist before work so an abort or killed solver consumes this trial
                            retry.hit=iHit;
                            retry.nTrials=checkedAdd(retry.nTrials,uint64_t{1},
                                                     "NICE retry counter overflow");
                            writeCheckpointRetryState(program,rngLayout,iEntry,retry,lock);
                            MASTER_PRINTF("NICE: trial %llu for hit %zu of gauge configuration %zu (limit %d, zero means unlimited)\n",
                                          static_cast<unsigned long long>(retry.nTrials),iHit,iEntry,
                                          program.options.policy.maxTrials);
                            
                            const TaskValues values=
                              executePreparedHit(program,rng,rngLayout,iEntry,iHit);
                            updater.requireHealthy();
                            const std::string hitValues=
                              encodeHitValues(program.options,program.correlations,values);
                            
                            requireCheckpointOwnership(lock);
                            commitCheckpointFile(hitCheckpointPath(directory,program.options,iHit),
                              hitCheckpointBytes(checkpointEntryIdentity(program,iEntry),iHit,hitValues),lock);
                            wallTime.recordCompletedHit();
                            
                            // The hit file is authoritative if a crash precedes this reset
                            retry={iHit+1,0};
                            writeCheckpointRetryState(program,rngLayout,iEntry,retry,lock);
                            
                            MASTER_PRINTF("NICE: committed hit %zu for gauge configuration %zu\n",
                                          iHit,iEntry);
                          }
                      }
                  }
                
                if(trialLimitReached)
                  {
                    updater.stop();
                    releaseRunningFile(lock);
                    MASTER_PRINTF("NICE: skipping gauge configuration %zu '%s': maximum number of trials reached (%llu attempts at hit %zu, limit %d), output '%s'\n",
                                  iEntry,program.configurationEntries[iEntry].path.c_str(),
                                  static_cast<unsigned long long>(retry.nTrials),firstHit,
                                  program.options.policy.maxTrials,directory.string().c_str());
                    outcome=ConfigurationRunOutcome::TrialLimitReached;
                  }
                else if(stopped)
                  {
                    updater.stop();
                    releaseRunningFile(lock);
                    MASTER_PRINTF("NICE: stopping before the next hit of gauge configuration %zu\n",iEntry);
                    outcome=ConfigurationRunOutcome::Stopped;
                  }
                else
                  {
                    updater.requireHealthy();
                    finalizeConfiguration(program,iEntry,rngLayout,lock);
                    updater.stop();
                    releaseRunningFile(lock);
                  }
              }
          }
	
        clearGaugeFields(program.gaugeRecipes);
      }
    
    return outcome;
  }
}

#endif
