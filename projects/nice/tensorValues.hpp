#ifndef _NICE_TENSOR_VALUES_HPP
#define _NICE_TENSOR_VALUES_HPP

#include "tensorComps.hpp"
#include "base/memory_manager.hpp"

#include <limits>
#include <utility>

namespace nissa::nice
{
  /// Borrow one real or complex scalar, whose numerical components may be separated
  struct ScalarView
  {
    /// Address of the real component
    double* data;
    
    /// Distance in doubles between real and imaginary components
    size_t componentStride;
    
    /// One component for real scalars, two for complex scalars
    size_t scalarWidth;
    
    /// Memory space containing the borrowed components
    MemorySpace space;
    
    /// Locate a numerical component without dereferencing device memory on the host
    CUDA_HOST_AND_DEVICE double* getPtrTo(const size_t component) const
    {
      if(component>=scalarWidth)
        CRASH("NICE scalar view component is outside its scalar type");
      
      return data+component*componentStride;
    }
    
    /// Access a numerical component on the processor owning this memory space
    CUDA_HOST_AND_DEVICE double& operator[](const size_t component) const
    {
      if(space!=currentMemorySpace)
        CRASH("NICE scalar view accessed from the wrong memory space");
      
      return *getPtrTo(component);
    }
    
    /// Access the real component through the same interface for both scalar types
    CUDA_HOST_AND_DEVICE double& real() const
    {
      return (*this)[0];
    }
    
    /// Access the imaginary component of a complex scalar
    CUDA_HOST_AND_DEVICE double& imag() const
    {
      return (*this)[1];
    }
  };
  
  /// Borrow a flattened tensor while preserving the physical scalar-component order
  struct TensorView
  {
    /// Start of the borrowed allocation
    double* data;
    
    /// Number of logical real or complex elements
    size_t nElements;
    
    /// Number of real-valued components per logical element
    size_t scalarWidth;
    
    /// Adjacent logical elements per real or imaginary plane
    size_t scalarComponentStride;
    
    /// Memory space containing the allocation
    MemorySpace space;
    
    /// Borrow one scalar without assuming that its real and imaginary parts are adjacent
    CUDA_HOST_AND_DEVICE ScalarView operator[](const size_t element) const
    {
      if(element>=nElements)
        CRASH("NICE tensor view index is outside the tensor");
      
      const size_t plane=element/scalarComponentStride;
      const size_t coordinate=element%scalarComponentStride;
      
      return {data+plane*scalarComponentStride*scalarWidth+coordinate,
              scalarComponentStride,scalarWidth,space};
    }
  };
  
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
      if((nDoubles/(this->descriptor.scalar==ScalarKind::Complex?2:1))%
         this->descriptor.scalarComponentStride!=0)
        CRASH("NICE tensor scalar planes do not divide its element count");
      
      if(nDoubles>static_cast<size_t>(std::numeric_limits<Size>::max()))
        CRASH("NICE tensor allocation exceeds memory-manager size range");
      
      data=memoryManager(space)->provide<double>(static_cast<Size>(nDoubles));
    }
    
    /// Borrow the tensor allocation without transferring ownership
    TensorView view() const
    {
      const size_t scalarWidth=descriptor.scalar==ScalarKind::Complex?2:1;
      
      return {data,nDoubles/scalarWidth,scalarWidth,descriptor.scalarComponentStride,space};
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
