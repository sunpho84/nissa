#ifndef _NICE_CHECKPOINT_HPP
#define _NICE_CHECKPOINT_HPP

#include "nice.hpp"
#include "checkpointIdentity.hpp"
#include "programExecution.hpp"
#include "rngLayout.hpp"

#include "io/input.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "routines/fnv1a64.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
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
                for(size_t iCopy=0;iCopy<static_cast<size_t>(options.nCopies);iCopy++)
                  for(int deltaTime=0;deltaTime<glbSize[0];deltaTime++)
                    {
                      const size_t index=
                        (iGammaPair*static_cast<size_t>(glbSize[0])+deltaTime)*
                        static_cast<size_t>(options.nCopies)+iCopy;
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
          for(size_t iCopy=0;iCopy<static_cast<size_t>(options.nCopies);iCopy++)
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
  
  /// One immutable binary hit file, named by its stable zero-based hit index
  inline std::filesystem::path hitCheckpointPath(const std::filesystem::path& directory,
                                          const size_t iHit)
  {
    return directory/std::format("hit-{:08}.bin",iHit);
  }
  
  /// Number of scientific rows required from one complete hit
  inline size_t expectedHitRows(const PreparedProgram& program)
  {
    size_t rows=0;
    
    for(const OutputBlock& output : program.correlations.outputs)
      rows+=output.traces.size()*output.gammaList.size()*
        static_cast<size_t>(program.options.nCopies)*static_cast<size_t>(program.options.T);
    
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
    const std::filesystem::path output=directory/"correlations.nice";
    const std::filesystem::path temporary=output.string()+".tmp";
    FILE* file=open_file(temporary.string(),"wb");
    
    writeCheckpointBytes(file,temporary,
                         std::format("# gauge configuration {}\n",entry.path));
    writeCheckpointBytes(file,temporary,
                         "# hit block trace sinkGamma sourceGamma copy deltaT real imaginary\n");
    
    for(size_t iHit=0;iHit<static_cast<size_t>(program.options.nHits);iHit++)
      writeCheckpointBytes(file,temporary,
                           formatHitResults(program.options,program.correlations,
                                            checkedHitValues(hitCheckpointPath(directory,iHit),
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
    commitCheckpointFile(directory/"finished.nice",
                         std::format("NICE-finished-v3 {} {} {} {}\n",
                                     checkpointEntryIdentity(program,iEntry),
                                     program.options.nHits,
                                     checkpointMetadataHash(program,rngLayout,iEntry),
                                     outputHash),lock);
    
    if(not program.options.keepCheckpoints)
      {
        for(size_t iHit=0;iHit<static_cast<size_t>(program.options.nHits);iHit++)
          {
            requireCheckpointOwnership(lock);
            removeFile(hitCheckpointPath(directory,iHit).string());
          }
        
        requireCheckpointOwnership(lock);
        removeFile((directory/"checkpoint.bin").string());
      }
    
    requireCheckpointOwnership(lock);
    removeFile(lock.lockPath());
    MASTER_PRINTF("NICE: wrote results to '%s'\n",output.string().c_str());
  }
  
  /// Check a finished marker and its result before skipping numerical work
  inline bool configurationFinished(const PreparedProgram& program,
                                    const RngLayout& rngLayout,
                                    const size_t iEntry)
  {
    const std::filesystem::path directory=
      program.configurationEntries[iEntry].outputDirectory;
    const std::filesystem::path marker=directory/"finished.nice";
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
           nHits<=0 or nHits>program.options.nHits or
           metadataHash!=checkpointMetadataHash(program,rngLayout,iEntry))
          CRASH("NICE finished marker '%s' does not match this run",marker.string().c_str());
        
        const std::filesystem::path output=directory/"correlations.nice";
        if(checkpointFileHash(output)!=outputHash)
          CRASH("NICE finished result '%s' is damaged",output.string().c_str());
        
        if(nHits<program.options.nHits and
           not fileExists((directory/"checkpoint.bin").string()))
          CRASH("NICE cannot increase hits in '%s': completed hit checkpoints were removed",
                directory.string().c_str());
        
        finished=nHits==program.options.nHits;
      }
    
    return finished;
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
                  options.L,options.T,options.seed,
                  options.nCopies,options.diluteSpin,options.diluteColor,
                  static_cast<int>(options.temporalBoundary),
                  options.drawGlobalOriginTask.enabled,
                  options.drawGlobalOriginTask.fixedOrigin[0],
                  options.drawGlobalOriginTask.fixedOrigin[1],
                  options.drawGlobalOriginTask.fixedOrigin[2],
                  options.drawGlobalOriginTask.fixedOrigin[3],
                  options.reservedRngDoublesPerSite,
                  program.correlations.instructions.size(),expectedHitRows(program),
                  rngLayout.origin.doubleDrawCost);
    
    for(const RngAllocation& allocation : rngLayout.gaugeOperations)
      metadata+=std::format("gauge-cost {}\n",allocation.doubleDrawCost);
    
    return metadata;
  }
  
  /// Encode a versioned binary checkpoint with fixed little-endian 64-bit fields
  inline std::string checkpointMetadataBytes(const PreparedProgram& program,
                                             const RngLayout& rngLayout,
                                             const size_t iEntry)
  {
    std::string bytes("NICECP3\0",8);
    appendCheckpointUint64(bytes,checkpointEntryIdentity(program,iEntry));
    appendCheckpointUint64(bytes,checkpointMetadataHash(program,rngLayout,iEntry));
    
    return bytes;
  }
  
  /// Refuse incompatible resume state, creating metadata only with an ownership claim
  inline void validateCheckpointMetadata(const PreparedProgram& program,
					 const RngLayout& rngLayout,
					 const size_t iEntry,
					 lock_file_t<>* lock=nullptr)
  {
    const std::filesystem::path directory=
      program.configurationEntries[iEntry].outputDirectory;
    const std::filesystem::path path=directory/"checkpoint.bin";
    const std::string expected=checkpointMetadataBytes(program,rngLayout,iEntry);
    if(fileExists(path.string()))
      {
        if(readCheckpointFile(path)!=expected)
          CRASH("NICE checkpoint metadata '%s' does not match this run",path.string().c_str());
      }
    else
      {
        if(fileExists(hitCheckpointPath(directory,0).string()) or
           fileExists((directory/"checkpoint.nice").string()) or
           fileExists((directory/"finished.nice").string()) or
           fileExists((directory/"correlations.nice").string()))
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
    size_t first=static_cast<size_t>(program.options.nHits);
    
    for(size_t iHit=0;iHit<static_cast<size_t>(program.options.nHits);iHit++)
      {
        const std::filesystem::path path=hitCheckpointPath(directory,iHit);
        if(fileExists(path.string()))
          {
            if(first!=static_cast<size_t>(program.options.nHits))
              CRASH("NICE checkpoint '%s' follows an unfinished hit",path.string().c_str());
            
            checkedHitValues(path,checkpointEntryIdentity(program,iEntry),iHit,
                             expectedHitRows(program));
          }
        else if(first==static_cast<size_t>(program.options.nHits))
          first=iHit;
      }
    
    return first;
  }
  
  /// Resume one configuration collectively from its first missing hit
  inline void runCheckpointedConfiguration(PreparedProgram& program,
                                    const RngLayout& rngLayout,
                                    const Sitmo::Rng& rng,
                                    const size_t iEntry)
  {
    const std::filesystem::path directory=
      program.configurationEntries[iEntry].outputDirectory;
    if(create_directories(directory.string())!=0)
      CRASH("NICE cannot create output directory '%s'",directory.string().c_str());
    
    if(configurationFinished(program,rngLayout,iEntry))
      MASTER_PRINTF("NICE: skipping completed gauge configuration %zu\n",iEntry);
    else
      {
        lock_file_t<> lock;
        lock.init();
        if(not lock.try_lock((directory/".nice.lock").string()))
          CRASH("NICE cannot create ownership tag in '%s'",directory.string().c_str());
	
        // The initial recipe provides the filesystem synchronization before checking ownership
        const bool claimed=
          prepareConfiguration(program,iEntry,rng,rngLayout,iEntry,
                               [&lock]()
                               {
                                 return lock.check_lock();
                               });
        if(not claimed)
          MASTER_PRINTF("NICE: skipping gauge configuration %zu after lost ownership of '%s'\n",
                        iEntry,lock.lockPath().c_str());
        else if(configurationFinished(program,rngLayout,iEntry))
          {
            MASTER_PRINTF("NICE: skipping completed gauge configuration %zu\n",iEntry);
          }
        else
          {
            validateCheckpointMetadata(program,rngLayout,iEntry,&lock);
            const size_t firstHit=firstUnfinishedHit(program,iEntry);
	    
            if(firstHit==0 and fileExists((directory/"correlations.nice").string()))
              CRASH("NICE output '%s' exists without checkpoints; refusing to overwrite it",
                    (directory/"correlations.nice").string().c_str());
	    
            if(firstHit<static_cast<size_t>(program.options.nHits))
              {
                for(size_t iHit=firstHit;
                    iHit<static_cast<size_t>(program.options.nHits);iHit++)
                  {
                    const TaskValues values=
                      executePreparedHit(program,rng,rngLayout,iEntry,iHit);
                    const std::string hitValues=
                      encodeHitValues(program.options,program.correlations,values);
		    
                    requireCheckpointOwnership(lock);
                    commitCheckpointFile(hitCheckpointPath(directory,iHit),
                      hitCheckpointBytes(checkpointEntryIdentity(program,iEntry),iHit,hitValues),lock);
		    
                    MASTER_PRINTF("NICE: committed hit %zu for gauge configuration %zu\n",
                                  iHit,iEntry);
                  }
              }
	    
            finalizeConfiguration(program,iEntry,rngLayout,lock);
          }
	
        clearGaugeFields(program.gaugeRecipes);
      }
  }
}

#endif
