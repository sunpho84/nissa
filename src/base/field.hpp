#ifndef _FIELD_HPP
#define _FIELD_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <cstddef>
#include <type_traits>

#include <base/bench.hpp>
#include <base/metaprogramming.hpp>
#include <base/vectors.hpp>
#include <communicate/communicate.hpp>
#include <geometry/geometry_lx.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  /// Subscription of the field
  ///
  /// Forward declaration
  template <typename F,
	    typename P>
  struct SubscribedField;
  
  /////////////////////////////////////////////////////////////////
  
  /// Memory layout
  enum FieldLayout{CPU_LAYOUT,GPU_LAYOUT};
  
  /// Coverage of the field
  enum SpaceTimeCoverage{FULL_SPACE,HALF_SPACE};
  
  /// Has or not the halo
  enum HaloPresence{WITHOUT_HALO,WITH_HALO};
  
  /// Predefinite memory layout
  constexpr FieldLayout DefaultFieldLayout=GPU_LAYOUT;
  
  /////////////////////////////////////////////////////////////////
  
  /// Number of sites contained in the field
  template <SpaceTimeCoverage spaceTimeCoverage,
	    HaloPresence haloPresence>
  struct FieldSizes
  {
    /// Number of sites covered by the field
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr int& nSites()
    {
      if constexpr(spaceTimeCoverage==FULL_SPACE)
	return locVol;
      else
	return locVolh;
    }
    
    /// Number of sites in the halo of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr int& nHaloSites()
    {
      if constexpr(spaceTimeCoverage==FULL_SPACE)
	return bord_vol;
      else
	return bord_volh;
    }
    
    /// Number of sites to be allocated
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static int nSitesToAllocate()
    {
      if constexpr(haloPresence==WITH_HALO)
	return nSites()+nHaloSites();
      else
	return nSites();
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Usable to recognize a field
  template <typename F>
  struct FieldFeat
  {
  };
  
  /// Field
  template <typename T,
	    SpaceTimeCoverage SC=FULL_SPACE,
	    HaloPresence HP=WITHOUT_HALO,
	    FieldLayout FL=DefaultFieldLayout>
  struct Field :
    FieldFeat<Field<T,SC,HP,FL>>,
    FieldSizes<SC,HP>
  {
    /// Name of the field
    const char* name;
    
    /// Fundamental type
    using Fund=
      std::remove_all_extents_t<T>;
    
    /// Components
    using Comps=T;
    
    /// Spacetime coverage
    static constexpr SpaceTimeCoverage spaceTimeCoverage=SC;
    
    /// Presence of the halo
    static constexpr HaloPresence haloPresence=HP;
    
    /// Memory layout of the field
    static constexpr FieldLayout fieldLayout=FL;
    
    /// Number of degrees of freedom
    static constexpr int nInternalDegs=
      sizeof(Comps)/sizeof(Fund);
    
    /// Size externally visible
    const int externalSize;
    
    /// Container for actual data
    Fund* data;
    
    /// States whether the border is allocated
    const bool borderIsAllocated;
    
    /// States whether the border is updated
    mutable bool borderIsValid;
    
    /// Computes the index of the data
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    int index(const int& site,
	      const int& internalDeg) const
    {
      if constexpr(FL==CPU_LAYOUT)
	return internalDeg+nInternalDegs*site;
      else
	return site+externalSize*internalDeg;
    }
    
    /// Constructor
    Field(const char *name,
	  const int& externalSize) :
      name(name),
      externalSize(externalSize),
      borderIsAllocated(externalSize==bord_vol or externalSize==bord_volh),
      borderIsValid(false)
    {
      master_printf("Allocating field\n");
      data=nissa_malloc(name,externalSize*nInternalDegs,Fund);
    }
    
    /// Destructor
    ~Field()
    {
      master_printf("Deallocating field\n");
      nissa_free(data);
    }
    
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    decltype(auto) operator[](const int& site) CONST			\
    {									\
      if constexpr(not std::is_array_v<T>)				\
	return								\
	  data[site];							\
      else								\
	if constexpr(FL==CPU_LAYOUT)					\
	  return ((CONST T*)data)[site];				\
	else								\
	  return							\
	    SubscribedField<CONST Field,				\
	    std::remove_extent_t<T>>(*this,site,nullptr);		\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* not const */);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    /// Set borders as invalid
    void set_borders_invalid()
    {
      borderIsValid=false;
    }
    
    /// Fill the sending buf using the data inside an lx vec
    void fill_sending_buf_with_lx_surface() const
    {
      verbosity_lv3_master_printf("filling filling filling\n");
      
      NISSA_PARALLEL_LOOP(ibord,0,bord_vol)
	for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	  ((Fund*)send_buf)[internalDeg+nInternalDegs*ibord]=
	    data[index(surflxOfBordlx[ibord],internalDeg)];
      NISSA_PARALLEL_LOOP_END;
    }
    
    void fill_lx_bord_with_receiving_buf() const
    {
      // for(int ibord=0;ibord<bord_vol;ibord++)
      NISSA_PARALLEL_LOOP(ibord,0,bord_vol)
	for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	  data[index(locVol+ibord,internalDeg)]=
	    ((Fund*)recv_buf)[internalDeg+nInternalDegs*ibord];
      NISSA_PARALLEL_LOOP_END;
    }
    
    /// Start communication using an lx border
    std::vector<MPI_Request> start_communicating_lx_borders() const
    {
      std::vector<MPI_Request> requests;
      
      if(not borderIsValid and nparal_dir>0)
	{
	  const size_t max=std::min(send_buf_size,recv_buf_size);
	  if(sizeof(T)*bord_vol>max)
	    crash("asking to create a communicator that needs %d large buffer (%d allocated)",
		  sizeof(T)*bord_vol,max);
	  
	  //take time and write some debug output
	  START_TIMING(tot_comm_time,ntot_comm);
	  verbosity_lv3_master_printf("Start communication of lx borders of %s\n",name);
	  
	  //fill the communicator buffer, start the communication and take time
	  fill_sending_buf_with_lx_surface();
          requests=comm_start(false);
          STOP_TIMING(tot_comm_time);
	}
      
      return requests;
    }
    
    void finish_communicating_lx_borders(std::vector<MPI_Request> requests) const
    {
      if(not borderIsValid and nparal_dir>0)
	{
	  //take note of passed time and write some debug info
	  START_TIMING(tot_comm_time,ntot_comm);
	  verbosity_lv3_master_printf("Finish communication of lx borders of %s\n",name);
	  
	  //wait communication to finish, fill back the vector and take time
	  comm_wait(requests);
	  fill_lx_bord_with_receiving_buf();
	  STOP_TIMING(tot_comm_time);
	  
	  borderIsValid=true;
      }
    }
    
    /// Communicate the borders
    void communicate_lx_borders() const
    {
      if(not borderIsValid)
	{
	  verbosity_lv3_master_printf("Sync communication of lx borders of %s\n",name);
	  
	  const std::vector<MPI_Request> requests=
	    start_communicating_lx_borders();
	  finish_communicating_lx_borders(requests);
      }
    }
    
    /// Start the communications
    static std::vector<MPI_Request> comm_start(const bool& isEoField)
    {
      std::vector<MPI_Request> requests(2*2*nparal_dir);
      
      const int div_coeff=(isEoField==0)?1:2; //dividing coeff
      
      int nRequests=0;
      
      for(int bf=0;bf<2;bf++)
	for(int mu=0;mu<NDIM;mu++)
	  if(paral_dir[mu])
	    {
	      const size_t sendOffset=(bord_offset[mu]+bord_volh*(!bf))*sizeof(T)/div_coeff;
	      const size_t recvOffset=(bord_offset[mu]+bord_volh*bf)*sizeof(T)/div_coeff;
	      const size_t messageLength=bord_dir_vol[mu]*sizeof(T)/div_coeff;
	      const int messageTag=bf+2*mu;
	      
	      MPI_Irecv(recv_buf+recvOffset,messageLength,MPI_CHAR,rank_neigh [bf][mu],
			messageTag,cart_comm,&requests[nRequests++]);
	      MPI_Isend(send_buf+sendOffset,messageLength,MPI_CHAR,rank_neigh[!bf][mu],
			messageTag,cart_comm,&requests[nRequests++]);
	    }
      
      return requests;
    }
    
    /// Wait for communications to finish
    static void comm_wait(std::vector<MPI_Request> requests)
    {
      verbosity_lv3_master_printf("Entering MPI comm wait\n");
      
      MPI_Waitall(requests.size(),&requests[0],MPI_STATUS_IGNORE);
    }
  };
  
  /// Hack
  template <typename T,
	    typename F>
  void set_borders_invalid(FieldFeat<F>& field)
  {
    static_cast<F*>(&field)->set_borders_invalid();
  }
  
  template <typename F,
	    typename P>
  struct SubscribedField
  {
    F& f;
    
    using Fund=
      typename F::Fund;
    
    using Comps=
      typename F::Comps;
    
    const int site;
    
    const P* ptr;
    
#define PROVIDE_EVAL(CONST)						\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr			\
    CONST ConstIf<std::is_const_v<std::remove_reference_t<F>>,Fund>& eval(const int& i) CONST \
    {									\
      const int internalDeg=						\
	(int)(size_t)(&ptr[i])/sizeof(Fund);				\
      									\
      return f.data[f.index(site,internalDeg)];				\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/* not const */);
    
#undef PROVIDE_EVAL
    
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    decltype(auto) operator[](const int& i) CONST			\
    {									\
      if constexpr(std::is_array_v<P>)					\
	return								\
	  SubscribedField<F,std::remove_extent_t<P>>(f,site,ptr[i]);	\
      else								\
	return eval(i);							\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* not const */);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    SubscribedField(F& f,
		    const int& site,
		    const P* ptr) :
      f(f),site(site),ptr(ptr)
    {
    }
  };
}

#endif
