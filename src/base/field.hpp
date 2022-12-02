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
#include <linalgs/reduce.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  /// Start the communications of buffer interpreted as halo
  std::vector<MPI_Request> startBufHaloNeighExchange(const int& divCoeff,
						     const size_t& bps);
  
  /// Start the communications of buffer interpreted as halo
  template <typename T>
  std::vector<MPI_Request> startBufHaloNeighExchange(const int& divCoeff)
  {
    return startBufHaloNeighExchange(divCoeff,sizeof(T));
  }
  
  /// Start the communications of buffer interpreted as edges
  std::vector<MPI_Request> startBufEdgesNeighExchange(const int& divCoeff,
						      const size_t& bps);
  
  /// Start the communications of buffer interpreted as edges
  template <typename T>
  std::vector<MPI_Request> startBufEdgesNeighExchange(const int& divCoeff)
  {
    return startBufEdgesNeighExchange(divCoeff,sizeof(T));
  }
  
  /// Wait for communications to finish
  inline void waitAsyncCommsFinish(std::vector<MPI_Request> requests)
  {
    verbosity_lv3_master_printf("Entering MPI comm wait\n");
    
    MPI_Waitall(requests.size(),&requests[0],MPI_STATUS_IGNORE);
  }
  
  /// Communicates the buffers as halo
  template <typename T>
  void exchangeHaloNeighBuf(const int& divCoeff)
  {
    waitAsyncCommsFinish(startBufHaloNeighExchange<T>(divCoeff));
  }
  
  /// Communicates the buffers as edges
  template <typename T>
  void exchangeEdgesNeighBuf(const int& divCoeff)
  {
    waitAsyncCommsFinish(startBufEdgesNeighExchange<T>(divCoeff));
  }
  
  /////////////////////////////////////////////////////////////////
  
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
    /// Assert that the coverage is definite
    static void assertHasDefinedCoverage()
    {
      static_assert(sitesCoverage==FULL_SPACE or sitesCoverage==EVEN_SITES or sitesCoverage==ODD_SITES,"Trying to probe some feature of a field with unknwon coverage. If you are accessing an EvnOrOddField subfield, do it without subscribing them, or cast to the specific subspace");
    }
    
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
      assertHasDefinedCoverage();
      
      if constexpr(sitesCoverage==FULL_SPACE)
	return surflxOfBordlx[iHalo];
      else
	if constexpr(sitesCoverage==EVEN_SITES or sitesCoverage==ODD_SITES)
	  return surfeo_of_bordeo[sitesCoverage][iHalo];
    }
    
    /// Surface site of a site in the e
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static auto surfSiteOfEdgeSite(const int& iEdge)
    {
      assertHasDefinedCoverage();
      
      if constexpr(sitesCoverage==FULL_SPACE)
	return surflxOfEdgelx[iEdge];
      else
	if constexpr(sitesCoverage==EVEN_SITES or sitesCoverage==ODD_SITES)
	  return surfeo_of_edgeo[sitesCoverage][iEdge];
    }
    
#define PROVIDE_NEIGH(UD)						\
									\
    /* Neighbor in the UD orientation */				\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION				\
    static auto locNeigh ## UD(const int& site,				\
			       const int& mu)				\
    {									\
      assertHasDefinedCoverage();					\
      									\
      if constexpr(sitesCoverage==FULL_SPACE)				\
	return loclxNeigh ## UD[site][mu];				\
      else								\
	if constexpr(sitesCoverage==EVEN_SITES or			\
		     sitesCoverage==ODD_SITES)				\
	  return loceo_neigh ## UD[sitesCoverage][site][mu];		\
    }
    
    PROVIDE_NEIGH(dw);
    
    PROVIDE_NEIGH(up);
    
#undef PROVIDE_NEIGH
    
    /////////////////////////////////////////////////////////////////
    
    /// Global coordinates
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static auto glbCoord(const int& site,
			 const int& mu)
    {
      assertHasDefinedCoverage();
      
      if constexpr(sitesCoverage==FULL_SPACE)
	return glbCoordOfLoclx[site][mu];
      else
	if constexpr(sitesCoverage==EVEN_SITES or sitesCoverage==ODD_SITES)
	  return glbCoordOfLoclx[loclx_of_loceo[sitesCoverage][site]][mu];
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  PROVIDE_FEATURE(Field);
  
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
    
    /// Exec the operation f on each site and degree of freedom
    template <typename F>
    Field& forEachSiteDeg(const F& f)
    {
      NISSA_PARALLEL_LOOP(site,0,this->nSites())
	for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	  f((*this)(site,internalDeg),site,internalDeg);
      NISSA_PARALLEL_LOOP_END;
      
      invalidateHalo();
      
      return *this;
    }
    
#define PROVIDE_SELFOP(OP)						\
    Field& operator OP ## =(const Field& oth)				\
    {									\
      NISSA_PARALLEL_LOOP(site,0,this->nSites())			\
	for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)	\
	  (*this)(site,internalDeg) OP ## =oth(site,internalDeg);	\
      NISSA_PARALLEL_LOOP_END;						\
									\
      invalidateHalo();							\
      									\
      return *this;							\
    }
    
    PROVIDE_SELFOP(+);
    PROVIDE_SELFOP(-);
    
#undef PROVIDE_SELFOP
    
    /// Reset to 0
    void reset()
    {
      NISSA_PARALLEL_LOOP(site,0,this->nSites())
	for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	  (*this)(site,internalDeg)=0.0;
      NISSA_PARALLEL_LOOP_END;
      
      invalidateHalo();
    }
    
    /// Squared norm
    double norm2() const
    {
      Field<Fund,SC> buf("buf");
      
      NISSA_PARALLEL_LOOP(site,0,this->nSites())
	{
	  double s2=0.0;
	  for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	    s2+=sqr((*this)(site,internalDeg));
	  buf[site]=s2;
	}
      NISSA_PARALLEL_LOOP_END;
      
      double res;
      glb_reduce(&res,buf,this->nSites());
      
      return res;
    }
    
    /// (*this,out)
    void scalarProdWith(complex& out,
			const Field& oth) const
    {
      Field<complex,SC> buf("buf");
      
      using NT=Fund[nInternalDegs][2];
      const auto& l=castComponents<NT>();
      const auto& r=oth.castComponents<NT>();
      
      NISSA_PARALLEL_LOOP(site,0,this->nSites())
	{
	  complex c;
	  complex_put_to_zero(c);
	  for(int internalDeg=0;internalDeg<nInternalDegs/2;internalDeg++)
	    complex_summ_the_conj1_prod(c,l[site][internalDeg],oth[site][internalDeg]);
	  complex_copy(buf[site],c);
	}
      NISSA_PARALLEL_LOOP_END;
      
      complex res;
      glb_reduce(&res,buf,this->nSites());
    }
    
    /// Re((*this,out))
    double realPartOfScalarProdWith(const Field& oth) const
    {
      Field<double,SC> buf("buf");
      
      NISSA_PARALLEL_LOOP(site,0,this->nSites())
	{
	  double r=0;
	  for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	    r+=(*this)(site,internalDeg)*oth(site,internalDeg);
	  buf[site]=r;
	}
      NISSA_PARALLEL_LOOP_END;
      
      double res;
      glb_reduce(&res,buf,this->nSites());
      
      return res;
    }

#define PROVIDE_CASTS(CONST)						\
    /* Cast to a different sitesCoverage */				\
    template <SitesCoverage NSC,					\
	      bool Force=false>						\
    CONST Field<T,NSC,FL>& castSitesCoverage() CONST			\
    {									\
      if constexpr(not (Force or SC==EVEN_OR_ODD_SITES))		\
	static_assert(NSC==EVEN_SITES or NSC==ODD_SITES,		\
		      "incompatible sitesCoverage! Force the change if needed"); \
									\
      return *(CONST Field<T,NSC,FL>*)this;				\
    }									\
									\
    /* Cast to a different comps */					\
    template <typename NT,						\
	      bool Force=false>						\
    CONST Field<NT,SC,FL>& castComponents() CONST			\
    {									\
      static_assert(Force or sizeof(T)==sizeof(NT),			\
		    "incompatible components! Force the change if needed"); \
      									\
      return *(CONST Field<NT,SC,FL>*)this;				\
    }
    
    PROVIDE_CASTS(const);
    
    PROVIDE_CASTS(/* not const */);
    
#undef PROVIDE_CASTS
    
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
    
    /// Move constructor
    Field(Field&& oth) :
      name(oth.name),
      haloEdgesPresence(oth.haloEdgesPresence),
      externalSize(oth.externalSize),
      data(oth.data),
      haloIsValid(oth.haloIsValid),
      edgesAreValid(oth.edgesAreValid)
    {
      oth.data=nullptr;
    }
    
    /// Copy constructor
    Field(const Field& oth) :
      Field(oth.name,oth.haloEdgesPresence)
    {
      *this=oth;
    }
    
    /// Destructor
    ~Field()
    {
      master_printf("Deallocating field %s\n",name);
      if(data)
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
    
    /// Set edges as invalid
    void invalidateEdges()
    {
      edgesAreValid=false;
    }
    
    /// Set halo as invalid
    void invalidateHalo()
    {
      haloIsValid=false;
      
      invalidateEdges();
    }
    
    /// Fill the sending buf using the data with a given function
    template <typename F>
    void fillSendingBufWith(F& f,
			    const int& n) const
    {
      NISSA_PARALLEL_LOOP(i,0,n)
	for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	  ((Fund*)send_buf)[internalDeg+nInternalDegs*i]=
	    data[index(f(i),internalDeg)];
      NISSA_PARALLEL_LOOP_END;
    }
    
    /// Fill the sending buf using the data on the surface of a field
    void fillSendingBufWithSurface() const
    {
      fillSendingBufWith(Field::surfSiteOfHaloSite,bord_vol/divCoeff);
    }
    
    /// Fill the sending buf using the data on the surface edge
    void fillSendingBufWithEdgesSurface() const
    {
      fillSendingBufWith(Field::surfSiteOfEdgeSite,edge_vol/divCoeff);
    }
    
    /// Fill the surface using the data from the buffer
    template <typename B,
	      typename F>
    void fillSurfaceWithReceivingBuf(const F& f)
    {
      for(int bf=0;bf<2;bf++)
	for(int mu=0;mu<NDIM;mu++)
	  NISSA_PARALLEL_LOOP(iHaloOriDir,0,bord_dir_vol[mu]/divCoeff)
	    {
	      const int iHalo=
		bf*bord_volh/divCoeff+
		iHaloOriDir+bord_offset[mu]/divCoeff;
	      
	      const int iSurf=
		Field::surfSiteOfHaloSite(iHalo);
	      
	      f((*this)[iSurf],
		((B*)recv_buf)[iHalo],
		bf,
		mu);
	    }
      NISSA_PARALLEL_LOOP_END;
    }
    
    /// Fill the sending buf using the data with a given function
    void fillHaloOrEdgesWithReceivingBuf(const int& offset,
					 const int& n) const
    {
      NISSA_PARALLEL_LOOP(i,0,n)
	for(int internalDeg=0;internalDeg<nInternalDegs;internalDeg++)
	  data[index(offset+i,internalDeg)]=
	    ((Fund*)recv_buf)[internalDeg+nInternalDegs*i];
      NISSA_PARALLEL_LOOP_END;
    }
    
    /// Fills the halo with the received buffer
    void fillHaloWithReceivingBuf() const
    {
      assertHasHalo();
      
      fillHaloOrEdgesWithReceivingBuf(locVol/divCoeff,bord_vol/divCoeff);
    }
    
    /// Fills the halo with the received buffer
    void fillEdgesWithReceivingBuf() const
    {
      assertHasEdges();
      
      fillHaloOrEdgesWithReceivingBuf((locVol+bord_vol)/divCoeff,edge_vol/divCoeff);
    }
    
    /// Fills the sending buffer with the halo, compressing into elements of B using f
    template <typename B,
	      typename F>
    void fillSendingBufWithHalo(const F& f) const
    {
      for(int bf=0;bf<2;bf++)
	for(int mu=0;mu<NDIM;mu++)
	  NISSA_PARALLEL_LOOP(iHaloOriDir,0,bord_dir_vol[mu]/divCoeff)
	    {
	      const int iHalo=
		bf*bord_volh/divCoeff+
		iHaloOriDir+bord_offset[mu]/divCoeff;
	      
	      f(((B*)send_buf)[iHalo],
		(*this)[locVol/divCoeff+iHalo],
		bf,
		mu);
	    }
          NISSA_PARALLEL_LOOP_END;
    }
    
    /// Start the communications of halo
    static std::vector<MPI_Request> startHaloAsyincComm()
    {
      return startBufHaloNeighExchange<T>(divCoeff);
    }
    
    /// Start the communications of edges
    static std::vector<MPI_Request> startEdgesAsyincComm()
    {
      return startBufEdgesNeighExchange<T>(divCoeff);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Assert that can communicate n sites
    static void assertCanCommunicate(const int& n)
    {
      /// Needed size in the buffer
      const size_t neededBufSize=
	sizeof(T)*n;
      
      const size_t maxBufSize=
	std::min(send_buf_size,recv_buf_size);
      
      if(neededBufSize>maxBufSize)
	crash("asking to create a communicator that needs %d large buffer (%d allocated)",
		  neededBufSize,maxBufSize);
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Start communication using halo
    std::vector<MPI_Request> startCommunicatingHalo() const
    {
      /// Pending requests
      std::vector<MPI_Request> requests;
      
      assertCanCommunicate(Field::nHaloSites());
      
      //take time and write some debug output
      START_TIMING(tot_comm_time,ntot_comm);
      verbosity_lv3_master_printf("Start communication of borders of %s\n",name);
      
      //fill the communicator buffer, start the communication and take time
      fillSendingBufWithSurface();
      requests=startHaloAsyincComm();
      STOP_TIMING(tot_comm_time);
      
      return requests;
    }
    
    /// Finalize communications
    void finishCommunicatingHalo(std::vector<MPI_Request> requests) const
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
    
    /// Crash if the halo is not allocated
    void assertHasHalo() const
    {
      if(not (haloEdgesPresence>=WITH_HALO))
	crash("needs halo allocated!");
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Start communication using edges
    std::vector<MPI_Request> startCommunicatingEdges() const
    {
      /// Pending requests
      std::vector<MPI_Request> requests;
      
      assertCanCommunicate(Field::nEdgesSites());
      
      //take time and write some debug output
      START_TIMING(tot_comm_time,ntot_comm);
      verbosity_lv3_master_printf("Start communication of edges of %s\n",name);
      
      //fill the communicator buffer, start the communication and take time
      fillSendingBufWithEdgesSurface();
      requests=startEdgesAsyincComm();
      STOP_TIMING(tot_comm_time);
      
      return requests;
    }
    
    /// Finalize communications
    void finishCommunicatingEdges(std::vector<MPI_Request> requests) const
    {
      //take note of passed time and write some debug info
      START_TIMING(tot_comm_time,ntot_comm);
      verbosity_lv3_master_printf("Finish communication of edges of %s\n",name);
      
      //wait communication to finish, fill back the vector and take time
      waitAsyncCommsFinish(requests);
      fillEdgesWithReceivingBuf();
      STOP_TIMING(tot_comm_time);
      
      haloIsValid=true;
    }
    
    /// Crash if the edges are not allocated
    void assertHasEdges() const
    {
      if(not (haloEdgesPresence>=WITH_HALO_EDGES))
	crash("needs edges allocated!");
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Communicate the halo
    void updateHalo(const bool& force=false) const
    {
      if(force or not haloIsValid)
	{
	  verbosity_lv3_master_printf("Sync communication of halo of %s\n",name);
	  
	  const std::vector<MPI_Request> requests=
	    startCommunicatingHalo();
	  finishCommunicatingHalo(requests);
      }
    }
    
    /// Communicate the edges
    void updateEdges(const bool& force=false) const
    {
      updateHalo(force);
      
      if(force or not edgesAreValid)
	{
	  verbosity_lv3_master_printf("Sync communication of edges of %s\n",name);
	  
	  const std::vector<MPI_Request> requests=
	    startCommunicatingEdges();
	  finishCommunicatingEdges(requests);
      }
    }
    
    /// Compare
    INLINE_FUNCTION
    bool operator==(const Field& oth) const
    {
      return data==oth.data;
    }
    
    /// Negate comparison
    INLINE_FUNCTION
    bool operator!=(const Field& oth) const
    {
      return not (*this)==oth;
    }
    
    /// Assign
    INLINE_FUNCTION
    Field& operator=(const Field& oth)
    {
      if(this!=&oth)
	{
	  NISSA_PARALLEL_LOOP(i,0,this->nSites()*nInternalDegs)
	    data[i]=oth.data[i];
	  NISSA_PARALLEL_LOOP_END;
	  
	  invalidateHalo();
	}
      
      return *this;
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
    int operator()() const
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
    EoField(const char* name,
	    const HaloEdgesPresence& haloEdgesPresence=WITHOUT_HALO) :
      evenPart(name,haloEdgesPresence),
      oddPart(name,haloEdgesPresence)
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
    
    /// Update the edges of both parities
    INLINE_FUNCTION
    void updateEdges() const
    {
      forBothParities([this](const auto& par)
      {
	(*this)[par].updateEdges();
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
    
    /// Invalidate the edges of both parities
    INLINE_FUNCTION
    void invalidateEdges()
    {
      forBothParities([this](const auto& par)
      {
	(*this)[par].invalidateEdges();
      });
    }
    
    /// Compare
    INLINE_FUNCTION
    bool operator==(const EoField& oth) const
    {
      return evenPart==oth.evenPart and oddPart==oth.oddPart;
    }
    
    /// Negate comparison
    INLINE_FUNCTION
    bool operator!=(const EoField& oth) const
    {
      return not (*this)==oth;
    }
    
    /// Assign
    INLINE_FUNCTION
    EoField& operator=(const EoField& oth)
    {
      evenPart=oth.evenPart;
      oddPart=oth.oddPart;
      
      return *this;
    }
  };
  
  /// Loop on both parities
  template <typename F>
  void forBothParities(F&& f)
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
