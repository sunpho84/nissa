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
  enum SitesCoverage{EVEN_SITES,ODD_SITES,FULL_SPACE};
  
  /// Has or not the halo
  enum HaloPresence{WITHOUT_HALO,WITH_HALO};
  
  /// Predefinite memory layout
  constexpr FieldLayout DefaultFieldLayout=GPU_LAYOUT;
  
  /////////////////////////////////////////////////////////////////
  
  /// Number of sites contained in the field
  template <SitesCoverage sitesCoverage,
	    HaloPresence HP>
  struct FieldSizes
  {
    /// States whether the halo is present
    static constexpr bool haloIsPresent=
      (HP==WITH_HALO);
    
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
    
    /// Number of sites to be allocated
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static int nSitesToAllocate()
    {
      if constexpr(haloIsPresent)
	return nSites()+nHaloSites();
      else
	return nSites();
    }
    
    /// Surface site of a site in the halo
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static int surfSiteOfHaloSite(const int& iHalo)
    {
      if constexpr(sitesCoverage==FULL_SPACE)
	return surflxOfBordlx[iHalo];
      else
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
	    SitesCoverage SC=FULL_SPACE,
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
    
    /// Coverage of sites
    static constexpr SitesCoverage sitesCoverage=SC;
    
    /// Presence of the halo
    static constexpr HaloPresence haloPresence=HP;
    
    /// Memory layout of the field
    static constexpr FieldLayout fieldLayout=FL;
    
    /// Number of degrees of freedom
    static constexpr int nInternalDegs=
      sizeof(Comps)/sizeof(Fund);
    
    /// Total allocated sites
    const int externalSize;
    
    /// Container for actual data
    Fund* data;
    
    /// States whether the halo is updated
    mutable bool haloIsValid;
    
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
    Field(const char *name) :
      name(name),
      externalSize(this->nSitesToAllocate()),
      haloIsValid(false)
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
    
    void fillHaloWithReceivingBuf() const
    {
      NISSA_PARALLEL_LOOP(iHalo,0,this->nHaloSites())
	for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	  data[index(locVol+iHalo,internalDeg)]=
	    ((Fund*)recv_buf)[internalDeg+nInternalDegs*iHalo];
      NISSA_PARALLEL_LOOP_END;
    }
    
    /// Start the communications
    static std::vector<MPI_Request> startAsyincComm()
    {
      std::vector<MPI_Request> requests(2*2*nparal_dir);
      
      const int div_coeff=
	(sitesCoverage==FULL_SPACE)?1:2; //dividing coeff
      
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
    static void waitAsyncCommsFinish(std::vector<MPI_Request> requests)
    {
      verbosity_lv3_master_printf("Entering MPI comm wait\n");
      
      MPI_Waitall(requests.size(),&requests[0],MPI_STATUS_IGNORE);
    }
    
    /// Start communication using halo
    std::vector<MPI_Request> startCommunicatingHalo() const
    {
      std::vector<MPI_Request> requests;
      
      if(not haloIsValid and nparal_dir>0)
	{
	  const size_t max=std::min(send_buf_size,recv_buf_size);
	  if(sizeof(T)*this->nHaloSites()>max)
	    crash("asking to create a communicator that needs %d large buffer (%d allocated)",
		  sizeof(T)*this->nHaloSites(),max);
	  
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
