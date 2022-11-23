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
  enum SitesCoverage{EVEN_SITES,ODD_SITES,FULL_SPACE,EVEN_OR_ODD_SITES};
  
  /// Has or not the halo and the edges
  enum HaloEdgesPresence{WITHOUT_HALO,WITH_HALO,WITH_HALO_EDGES};
  
  /// Predefinite memory layout
  constexpr FieldLayout DefaultFieldLayout=GPU_LAYOUT;
  
  /////////////////////////////////////////////////////////////////
  
  /// Number of sites contained in the field
  template <SitesCoverage sitesCoverage>
  struct FieldSizes
  {
    /// Number of sites covered by the field
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr int& nSites()
    {
      if constexpr(sitesCoverage==FULL_SPACE)
	return locVol;
      else
	return locVolh;
    }
    
    /// Number of sites in the halo of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr int& nHaloSites()
    {
      if constexpr(sitesCoverage==FULL_SPACE)
	return bord_vol;
      else
	return bord_volh;
    }
    
    /// Number of sites in the edges of the field (not necessarily allocated)
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static constexpr int& nEdgesSites()
    {
      if constexpr(sitesCoverage==FULL_SPACE)
	return edge_vol;
      else
	return edge_volh;
    }
    
    /// Computes the size to allocate
    static int nSitesToAllocate(const HaloEdgesPresence& haloEdgesPresence)
    {
      int res=locVol;
      
      if(haloEdgesPresence>=WITH_HALO)
	res+=nHaloSites();
      
      if(haloEdgesPresence>=WITH_HALO_EDGES)
	res+=nEdgesSites();
      
      return res;
    }
    
    /// Surface site of a site in the halo
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static auto surfSiteOfHaloSite(const int& iHalo)
    {
      if constexpr(sitesCoverage==FULL_SPACE)
	return surflxOfBordlx[iHalo];
      else
	if constexpr(sitesCoverage==EVEN_SITES or sitesCoverage==ODD_SITES)
	  return surfeo_of_bordeo[sitesCoverage][iHalo];
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
	    SitesCoverage SC,
	    FieldLayout FL=DefaultFieldLayout>
  struct Field :
    FieldFeat<Field<T,SC,FL>>,
    FieldSizes<SC>
  {
    /// Coefficient which divides the space time, if the field is covering only half the space
    static constexpr const int divCoeff=
      (SC==FULL_SPACE)?1:2;
    
    /// Name of the field
    const char* name;
    
    /// Fundamental type
    using Fund=
      std::remove_all_extents_t<T>;
    
    /// Components
    using Comps=T;
    
    /// Coverage of sites
    static constexpr SitesCoverage sitesCoverage=SC;
    
    /// Memory layout of the field
    static constexpr FieldLayout fieldLayout=FL;
    
    /// Number of degrees of freedom
    static constexpr int nInternalDegs=
      sizeof(Comps)/sizeof(Fund);
    
    /// Presence of halo and edges
    const HaloEdgesPresence haloEdgesPresence;
    
    /// Total allocated sites
    const int externalSize;
    
    /// Container for actual data
    Fund* data;
    
    /// States whether the halo is updated
    mutable bool haloIsValid;
    
    /// States whether the edges are updated
    mutable bool edgesAreValid;
    
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
	  const HaloEdgesPresence& haloEdgesPresence=WITHOUT_HALO) :
      name(name),
      haloEdgesPresence(haloEdgesPresence),
      externalSize(FieldSizes<sitesCoverage>::nSitesToAllocate(haloEdgesPresence)),
      haloIsValid(false),
      edgesAreValid(false)
    {
      master_printf("Allocating field %s\n",name);
      data=nissa_malloc(name,externalSize*nInternalDegs,Fund);
    }
    
    /// Destructor
    ~Field()
    {
      master_printf("Deallocating field %s\n",name);
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
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_FLATTENED_CALLER(CONST)					\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    CONST Fund& operator()(const int& site,const int& internalDeg) CONST \
    {									\
      return data[index(site,internalDeg)];				\
    }
    
    PROVIDE_FLATTENED_CALLER(const);
    
    PROVIDE_FLATTENED_CALLER(/* not const */);
    
#undef PROVIDE_FLATTENED_CALLER
    
    /////////////////////////////////////////////////////////////////
    
    /// Set halo as invalid
    void invalidateHalo()
    {
      haloIsValid=false;
    }
    
    /// Fill the sending buf using the data inside a vec
    void fillSendingBufWithSurface() const
    {
      verbosity_lv3_master_printf("filling filling filling\n");
      
      NISSA_PARALLEL_LOOP(iHalo,0,this->nHaloSites())
	for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	  ((Fund*)send_buf)[internalDeg+nInternalDegs*iHalo]=
	    data[index(this->surfSiteOfHaloSite(iHalo),internalDeg)];
      NISSA_PARALLEL_LOOP_END;
    }
    
    /// Fills the halo with the receivied buffer
    void fillHaloWithReceivingBuf() const
    {
      NISSA_PARALLEL_LOOP(iHalo,0,this->nHaloSites())
	for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	  data[index(this->nSites()+iHalo,internalDeg)]=
	    ((Fund*)recv_buf)[internalDeg+nInternalDegs*iHalo];
      NISSA_PARALLEL_LOOP_END;
    }
    
    /// Start the communications
    static std::vector<MPI_Request> startAsyincComm()
    {
      /// Pending requests
      std::vector<MPI_Request> requests(2*2*nparal_dir);
      
      int nRequests=0;
      
      for(int bf=0;bf<2;bf++)
	for(int mu=0;mu<NDIM;mu++)
	  if(paral_dir[mu])
	    {
	      const size_t sendOffset=(bord_offset[mu]+bord_volh*(!bf))*sizeof(T)/divCoeff;
	      const size_t recvOffset=(bord_offset[mu]+bord_volh*bf)*sizeof(T)/divCoeff;
	      const size_t messageLength=bord_dir_vol[mu]*sizeof(T)/divCoeff;
	      const int messageTag=bf+2*mu;
	      
	      MPI_Irecv(recv_buf+recvOffset,messageLength,MPI_CHAR,rank_neigh [bf][mu],
			messageTag,cart_comm,&requests[nRequests++]);
	      MPI_Isend(send_buf+sendOffset,messageLength,MPI_CHAR,rank_neigh[!bf][mu],
			messageTag,cart_comm,&requests[nRequests++]);
	    }
      
      return requests;
    }
    
    /// Wait for communications to finish
    static void waitAsyncCommsFinish(std::vector<MPI_Request> requests)
    {
      verbosity_lv3_master_printf("Entering MPI comm wait\n");
      
      MPI_Waitall(requests.size(),&requests[0],MPI_STATUS_IGNORE);
    }
    
    /// Start communication using halo
    std::vector<MPI_Request> startCommunicatingHalo() const
    {
      /// Pending requests
      std::vector<MPI_Request> requests;
      
      if(not haloIsValid and nparal_dir>0)
	{
	  /// Needed size in the buffer
	  const size_t neededBufSize=
	    sizeof(T)*this->nHaloSites();
	  
	  const size_t maxBufSize=
	    std::min(send_buf_size,recv_buf_size);
	  
	  if(neededBufSize>maxBufSize)
	    crash("asking to create a communicator that needs %d large buffer (%d allocated)",
		  neededBufSize,maxBufSize);
	  
	  //take time and write some debug output
	  START_TIMING(tot_comm_time,ntot_comm);
	  verbosity_lv3_master_printf("Start communication of borders of %s\n",name);
	  
	  //fill the communicator buffer, start the communication and take time
	  fillSendingBufWithSurface();
          requests=startAsyincComm();
          STOP_TIMING(tot_comm_time);
	}
      
      return requests;
    }
    
    /// Finalize communications
    void finishCommunicatingHalo(std::vector<MPI_Request> requests) const
    {
      if(not haloIsValid and nparal_dir>0)
	{
	  //take note of passed time and write some debug info
	  START_TIMING(tot_comm_time,ntot_comm);
	  verbosity_lv3_master_printf("Finish communication of halo of %s\n",name);
	  
	  //wait communication to finish, fill back the vector and take time
	  waitAsyncCommsFinish(requests);
	  fillHaloWithReceivingBuf();
	  STOP_TIMING(tot_comm_time);
	  
	  haloIsValid=true;
      }
    }
    
    /// Communicate the halo
    void updateHalo() const
    {
      if(not haloIsValid)
	{
	  verbosity_lv3_master_printf("Sync communication of halo of %s\n",name);
	  
	  const std::vector<MPI_Request> requests=
	    startCommunicatingHalo();
	  finishCommunicatingHalo(requests);
      }
    }
  };
  
  /// Hack
  template <typename F>
  void set_borders_invalid(FieldFeat<F>& field)
  {
    static_cast<F*>(&field)->invalidateHalo();
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Lexicographic field
  template <typename T,
	    FieldLayout FL=DefaultFieldLayout>
  using LxField=Field<T,FULL_SPACE,FL>;
  
  /// Field over even sites
  template <typename T,
	    FieldLayout FL=DefaultFieldLayout>
  using EvnField=Field<T,EVEN_SITES,FL>;
  
  /// Field over odd sites
  template <typename T,
	    FieldLayout FL=DefaultFieldLayout>
  using OddField=Field<T,ODD_SITES,FL>;
  
  /// Field over even or odd sites
  template <typename T,
	    FieldLayout FL=DefaultFieldLayout>
  using EvenOrOddField=Field<T,EVEN_OR_ODD_SITES,FL>;
  
  /////////////////////////////////////////////////////////////////
  
  template <int EO>
  struct Par
  {
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    operator int() const
    {
      return EO;
    }
    
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Par<1-EO> operator!() const
    {
      return {};
    }
  };
  
  /// Structure to hold an even/old field
  template <typename T,
	    FieldLayout FL=DefaultFieldLayout>
  struct EoField
  {
    /// Type representing a pointer to type T
    template <SitesCoverage EO>
    using F=Field<T,EO,FL>;
    
    F<EVEN_SITES> evenPart;
    
    F<ODD_SITES> oddPart;
    
    /////////////////////////////////////////////////////////////////
    
    /// Gets the even or odd part
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    template <int EO>							\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr			\
    CONST auto& operator[](Par<EO>) CONST				\
    {									\
      if constexpr(EO==EVN)						\
	return evenPart;						\
      else								\
	return oddPart;							\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* const*/ );
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    /////////////////////////////////////////////////////////////////
    
    /// Gets the even or odd part, non-constexpr version
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION				\
    CONST auto& operator[](const int eo) CONST				\
    {									\
      using EOOF=							\
	CONST Field<T,EVEN_OR_ODD_SITES,FL>;				\
      									\
      EOOF* t[2]={(EOOF*)&evenPart,(EOOF*)&oddPart};			\
      									\
      return *t[eo];							\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* const*/ );
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    /////////////////////////////////////////////////////////////////
    
    /// Constructor
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    EoField(const char* name) :
      evenPart(name),
      oddPart(name)
    {
    }
    
    /// Update the halo of both parities
    INLINE_FUNCTION
    void updateHalo() const
    {
      forBothParities([this](const auto& par)
      {
	(*this)[par].updateHalo();
      });
    }
    
    /// Invalidate the halo of both parities
    INLINE_FUNCTION
    void invalidateHalo()
    {
      forBothParities([this](const auto& par)
      {
	(*this)[par].invalidateHalo();
      });
    }
  };
  
  /// Loop on both parities
  template <typename F>
  void forBothParities(const F& f)
  {
    f(Par<0>{});
    f(Par<1>{});
  }
  
  /////////////////////////////////////////////////////////////////
  
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
