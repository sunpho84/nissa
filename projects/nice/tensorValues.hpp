#ifndef _NICE_TENSOR_VALUES_HPP
#define _NICE_TENSOR_VALUES_HPP

#include "tensorComps.hpp"
#include "base/memory_manager.hpp"

#include <limits>
#include <utility>

namespace nissa::nice
{
  /// Owned allocation of scalar elements in descriptor component order
  struct Tensor
  {
    /// Storage components and scalar type
    TensorDescriptor descriptor;
    
    /// Allocator space for the underlying doubles
    MemorySpace space;
    
    /// Number of doubles, twice the element count for complex tensors
    size_t nDoubles;
    
    /// Direct memory-manager allocation, moved rather than copied
    double* data;
    
    /// Allocate the tensor in the requested NISSA memory space
    Tensor(TensorDescriptor descriptor,
           const MemorySpace space) :
      descriptor(std::move(descriptor)),
      space(space),
      nDoubles(checkedMultiply(this->descriptor.nElements(),
                               this->descriptor.scalar==ScalarKind::Complex?size_t{2}:size_t{1},
                               "NICE tensor storage size overflow")),
      data(nullptr)
    {
      if(nDoubles>static_cast<size_t>(std::numeric_limits<Size>::max()))
        CRASH("NICE tensor allocation exceeds memory-manager size range");
      
      data=memoryManager(space)->provide<double>(static_cast<Size>(nDoubles));
    }
    
    /// Prevent accidental ownership duplication
    Tensor(const Tensor&)=delete;
    
    /// Move one allocation without copying its contents
    Tensor(Tensor&& other) :
      descriptor(std::move(other.descriptor)),
      space(other.space),
      nDoubles(std::exchange(other.nDoubles,0)),
      data(std::exchange(other.data,nullptr))
    {
    }
    
    /// Prevent accidental ownership duplication
    Tensor& operator=(const Tensor&)=delete;
    
    /// Release the old allocation and take ownership of another
    Tensor& operator=(Tensor&& other)
    {
      if(this!=&other)
        {
          memoryManager(space)->release(data);
          descriptor=std::move(other.descriptor);
          space=other.space;
          nDoubles=std::exchange(other.nDoubles,0);
          data=std::exchange(other.data,nullptr);
        }
      
      return *this;
    }
    
    /// Release the owned allocation through the same memory manager
    ~Tensor()
    {
      if(data!=nullptr)
        memoryManager(space)->release(data);
    }
  };
}

#endif
