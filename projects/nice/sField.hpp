#ifndef _NICE_S_FIELD_HPP
#define _NICE_S_FIELD_HPP

#include "tensorValues.hpp"
#include "base/field.hpp"
#include "geometry/geometry_lx.hpp"

#include <cstddef>
#include <initializer_list>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace nissa::nice
{
  /// Borrow a real-valued projection with explicit site and inner-component strides
  struct SFieldComponentView
  {
    /// Address of the first selected real-valued component
    double* data;
    
    /// Number of allocated sites, including requested halos
    size_t nSites;
    
    /// Number of inner elements per site
    size_t nInner;
    
    /// Distance in doubles between neighboring site coordinates
    size_t siteStride;
    
    /// Distance in doubles between neighboring inner coordinates
    size_t innerStride;
    
    /// Memory space containing the borrowed projection
    MemorySpace space;
    
    /// Locate a projected value without dereferencing it
    CUDA_HOST_AND_DEVICE double* getPtrTo(const size_t site,
                                          const size_t inner=0) const
    {
      if(site>=nSites or inner>=nInner)
        CRASH("NICE field component view index is outside the selected field");
      
      return data+site*siteStride+inner*innerStride;
    }
    
    /// Access one projected value on the processor owning this memory space
    CUDA_HOST_AND_DEVICE double& operator()(const size_t site,
                                            const size_t inner=0) const
    {
      assertRunningOnMemorySpace(space);
      
      return *getPtrTo(site,inner);
    }
  };
  
  /// Non-owning generic view of one selected set of outer components
  struct SFieldView
  {
    /// Start of the selected outer-component block in tensor storage
    double* data;
    
    /// Number of allocated sites, including any requested halo sites
    size_t nSites;
    
    /// Number of scalar elements in the inner components
    size_t nInner;
    
    /// Number of doubles per scalar element
    size_t scalarWidth;
    
    /// Physical ordering of site and inner components
    SpaceTimeLayout layout;
    
    /// Memory space containing the borrowed block
    MemorySpace space;
    
    /// Borrow one scalar using the same fundamental-component order as library fields
    CUDA_HOST_AND_DEVICE ScalarView operator()(const size_t site,
                                               const size_t inner=0) const
    {
      if(site>=nSites or inner>=nInner)
        CRASH("NICE field reference index is outside the selected field");
      
      const size_t offset=layout==SpaceTimeLayout::CPU?
        (site*nInner+inner)*scalarWidth:inner*scalarWidth*nSites+site;
      const size_t componentStride=layout==SpaceTimeLayout::CPU?1:nSites;
      
      return {data+offset,componentStride,scalarWidth,space};
    }
    
    /// Locate one numerical component; adjacent real and imaginary parts are not assumed
    CUDA_HOST_AND_DEVICE double* getPtrTo(const size_t site,
                                          const size_t inner=0,
                                          const size_t component=0) const
    {
      return (*this)(site,inner).getPtrTo(component);
    }
    
    /// Project one numerical component without copying or changing the logical indices
    CUDA_HOST_AND_DEVICE SFieldComponentView component(const size_t index) const
    {
      if(index>=scalarWidth)
        CRASH("NICE field projection is outside its scalar type");
      
      const bool cpu=layout==SpaceTimeLayout::CPU;
      
      return {data+index*(cpu?1:nSites),nSites,nInner,
              cpu?nInner*scalarWidth:1,cpu?scalarWidth:nSites*scalarWidth,space};
    }
    
    /// Borrow the real-valued field projection
    CUDA_HOST_AND_DEVICE SFieldComponentView real() const
    {
      return component(0);
    }
    
    /// Borrow the imaginary projection of a complex field
    CUDA_HOST_AND_DEVICE SFieldComponentView imag() const
    {
      return component(1);
    }
  };
  
  /// Lattice field over one tensor, with site between outer and inner components
  struct SField
  {
    /// Components selected to obtain a complete lattice field
    std::vector<Comp> outer;
    
    /// Components traversed by a numerical field kernel
    std::vector<Comp> inner;
    
    /// Sites allocated in each outer-component block
    size_t nSites;
    
    /// Scalar elements per site in each outer-component block
    size_t nInner;
    
    /// Number of doubles per scalar element
    size_t scalarWidth;
    
    /// Physical ordering of site and inner components
    SpaceTimeLayout layout;
    
    /// Halo allocation requested for the field
    HaloEdgesPresence halo;
    
    /// Sole owner of all outer blocks and their site/inner data
    Tensor tensor;
    
    /// Get allocated site count only after lattice geometry is ready
    static size_t allocatedSites(const HaloEdgesPresence halo)
    {
      if(not lxGeomInited)
        CRASH("NICE field requires initialized geometry");
      
      return static_cast<size_t>(FieldSizes<FULL_SPACE>::nSitesToAllocate(halo));
    }
    
    /// Construct the physical tensor shape from outer, site, and inner components
    static TensorDescriptor storageDescriptor(const ScalarKind scalar,
                                              const std::vector<Comp>& inner,
                                              const std::vector<Comp>& outer,
                                              const size_t nSites,
                                              const SpaceTimeLayout layout)
    {
      TensorDescriptor result{scalar,outer};
      const Comp site{StandardComp::Site,CompSide::None,nSites};
      
      for(const Comp& comp : outer)
        if(comp.extent==0)
          CRASH("NICE field outer components must have positive extent");
      
      for(const Comp& comp : inner)
        if(comp.extent==0)
          CRASH("NICE field inner components must have positive extent");
      
      if(layout==SpaceTimeLayout::CPU)
        {
          result.comps.push_back(site);
          result.comps.insert(result.comps.end(),inner.begin(),inner.end());
        }
      else if(layout==SpaceTimeLayout::GPU)
        {
          result.comps.insert(result.comps.end(),inner.begin(),inner.end());
          result.comps.push_back(site);
          result.scalarComponentStride=scalar==ScalarKind::Complex?nSites:1;
        }
      else
        CRASH("NICE field has an unsupported site layout");
      
      result.validate();
      
      return result;
    }
    
    /// Count inner elements without making site or outer components part of them
    static size_t nInnerElements(const std::vector<Comp>& inner)
    {
      size_t count=1;
      
      for(const Comp& comp : inner)
        count=checkedMultiply(count,comp.extent,"NICE field inner size overflow");
      
      return count;
    }
    
    /// Allocate one tensor with a site component added automatically
    SField(const ScalarKind scalar,
          std::vector<Comp> inner,
          std::vector<Comp> outer,
          const MemorySpace space=defaultMemorySpace,
          const SpaceTimeLayout layout=defaultSpaceTimeLayout,
          const HaloEdgesPresence halo=WITHOUT_HALO) :
      outer(std::move(outer)),
      inner(std::move(inner)),
      nSites(allocatedSites(halo)),
      nInner(nInnerElements(this->inner)),
      scalarWidth(scalar==ScalarKind::Complex?2:1),
      layout(layout),
      halo(halo),
      tensor(storageDescriptor(scalar,this->inner,this->outer,nSites,layout),space)
    {
    }
    
    /// Borrow one field after selecting each outer-component coordinate
    SFieldView view(const std::span<const size_t> coordinates) const
    {
      if(coordinates.size()!=outer.size())
        CRASH("NICE field view needs one coordinate per outer component");
      
      size_t block=0;
      
      for(size_t i=0;i<outer.size();i++)
        {
          if(coordinates[i]>=outer[i].extent)
            CRASH("NICE field view outer component %zu is out of range",i);
          
          block=block*outer[i].extent+coordinates[i];
        }
      
      const size_t blockDoubles=nSites*nInner*scalarWidth;
      
      return {tensor.data+block*blockDoubles,nSites,nInner,scalarWidth,
              layout,tensor.space};
    }
    
    /// Select outer components by kind, orientation, and label in any order
    SFieldView view(const std::span<const CompIndex> indices) const
    {
      if(indices.size()!=outer.size())
        CRASH("NICE field view needs one index per outer component");
      
      std::vector<size_t> coordinates(outer.size());
      std::vector<bool> selected(outer.size(),false);
      
      for(const CompIndex& index : indices)
        {
          size_t component=0;
	  
          while(component<outer.size() and not outer[component].sameAs(index))
            component++;
	  
          if(component==outer.size())
            CRASH("NICE field view index does not identify an outer component");
	  
          if(selected[component])
            CRASH("NICE field view repeats an outer component index");
	  
          selected[component]=true;
          coordinates[component]=index.coordinate;
        }
      
      return view(std::span<const size_t>{coordinates});
    }
    
    /// Convenience overload for a list of typed outer-component indices
    SFieldView view(const std::initializer_list<CompIndex> indices) const
    {
      return view(std::span<const CompIndex>{indices.begin(),indices.size()});
    }
    
    /// Materialize a library field borrowing one selected outer block
    template<typename T,
             SpaceTimeLayout STL,
             MemorySpace MS>
    nissa::LxField<T,STL,MS> shadow(const char* name,
                                    const std::span<const size_t> coordinates) const
    {
      using Fund=std::remove_all_extents_t<T>;
      static_assert(std::is_same_v<Fund,double>);
      
      if(layout!=STL or tensor.space!=MS)
        CRASH("NICE field shadow has incompatible layout or memory space");
      
      if(sizeof(T)!=nInner*scalarWidth*sizeof(double))
        CRASH("NICE field shadow has incompatible inner component size");
      
      const SFieldView selected=view(coordinates);
      
      return nissa::LxField<T,STL,MS>{name,selected.data,halo};
    }
  };
  
  /// Allocate a line field with source columns and copies as outer components
  inline SField makeLineField(const size_t nCopies,
                             const bool diluteSpin,
                             const bool diluteColor,
                             const MemorySpace space=defaultMemorySpace,
                             const SpaceTimeLayout layout=defaultSpaceTimeLayout)
  {
    if(nCopies==0)
      CRASH("NICE line field requires positive copies");
    
    std::vector<Comp> inner{{StandardComp::Spin,CompSide::Row,NDIRAC},
                            {StandardComp::Color,CompSide::Row,NCOL}};
    std::vector<Comp> outer;
    
    if(diluteSpin)
      outer.push_back({StandardComp::Spin,CompSide::Column,NDIRAC});
    
    outer.push_back({StandardComp::Copy,CompSide::None,nCopies});
    
    if(diluteColor)
      outer.push_back({StandardComp::Color,CompSide::Column,NCOL});
    
    return SField{ScalarKind::Complex,std::move(inner),std::move(outer),
                  space,layout,WITH_HALO};
  }
  
}

#endif
