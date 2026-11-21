#ifndef _RECURRING_WORK_HPP
#define _RECURRING_WORK_HPP

#include "base/debug.hpp"

#include <chrono>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <utility>

namespace nissa
{
  /// Run periodic work on a background thread, waiting one interval before each call
  /// The caller selects participating ranks; callbacks must not call MPI
  /// Return false from the callback to stop; callbacks must finish before stop can join
  /// Start, stop, and destruction belong to the owning thread, never the callback
  struct RecurringWork
  {
    /// Thread executing the callback
    std::thread worker{};
    
    /// Protect the stop request
    std::mutex mutex{};
    
    /// Interrupt the wait when stopping
    std::condition_variable wake{};
    
    /// Whether the worker has been asked to stop
    bool stopping{true};
    
    /// Construct without starting a worker
    RecurringWork()=default;
    
    /// Stop previous work and start a callback with a strictly positive interval
    void start(const std::chrono::milliseconds interval,
               std::function<bool()> callback)
    {
      if(interval.count()<=0 or not callback)
        CRASH("Recurring work requires a positive interval and a callback");
      
      stop();
      stopping=false;
      worker=std::thread([this,interval,callback=std::move(callback)]()
      {
        std::unique_lock guard(mutex);
        
        while(not stopping)
          {
            const bool requestedStop=wake.wait_for(guard,interval,[this]()
            {
              return stopping;
            });
            
            if(not requestedStop)
              {
                guard.unlock();
                const bool keepRunning=callback();
                guard.lock();
                stopping=stopping or not keepRunning;
              }
          }
      });
    }
    
    /// Interrupt the wait and join any in-flight callback; safe to call repeatedly
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
    
    /// Join before captured callback state can be destroyed
    ~RecurringWork()
    {
      stop();
    }
  };
}

#endif
