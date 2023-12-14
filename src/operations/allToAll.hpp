#ifndef _ALL_TO_ALL_HPP
#define _ALL_TO_ALL_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/allToAll.hpp

#include <algorithm>
#include <map>
#include <vector>

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/dynamicTens.hpp>
#include <expr/mirroredNode.hpp>
#include <threads/threads.hpp>
#include <tuples/invokeWithTypesOfTuple.hpp>
#include <tuples/tupleCat.hpp>
#include <tuples/tupleReplaceType.hpp>

namespace nissa
{
  /// All to all communicators
  template <DerivedFromComp CDst,
	    DerivedFromComp CSrc,
	    bool IsRef=false>
  struct AllToAllComm
  {
    /// Status: is the communicator inited?
    bool inited;
    
    DECLARE_DYNAMIC_COMP(BufComp);
    
    /// Default constructor
    AllToAllComm() :
      inited{false}
    {
    }
    
    /// Default move constructor
    AllToAllComm(AllToAllComm&& oth)=default;
    
    /// Move assign
    AllToAllComm& operator=(AllToAllComm&& oth)=default;
    
    /// Construct and initialize
    template <typename F>
    AllToAllComm(const CSrc &nSrc,
		 F&& f) :
      AllToAllComm()
    {
      init(nSrc,f);
    }
    
    /// Destination of the source, in the buffer
    MirroredTens<CompsList<CSrc>,BufComp,IsRef> outBufOfSrc;
    
    /// Destination of the buffer, in the final storage
    MirroredTens<CompsList<BufComp>,CDst,IsRef> dstOfInBuf;
    
    /// Keeps track of the number of elements to be sent to each given rank
    std::vector<std::pair<MpiRank,int64_t>> nSendToRank;
    
    /// Keeps track of the number of elements to be received from each given rank
    std::vector<std::pair<MpiRank,int64_t>> nRecvFrRank;
    
    /// Gets the destination size, which is equal to the in buffer size
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    CDst getNDst() const
    {
      return getInBufSize()();
    }
    
    /// Gets the incoming buffer size
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    constexpr BufComp getInBufSize() const
    {
      return dstOfInBuf.template getCompSize<BufComp>();
    }
    
    /// Gets the source size, which is equal to the out buffer size
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    CSrc getNSrc() const
    {
      return getOutBufSize()();
    }
    
    /// Gets the outgoing buffer size
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    constexpr BufComp getOutBufSize() const
    {
      return outBufOfSrc.template getCompSize<CSrc>()();
    }
    
    // Verify that the destination is assigned only once
    void verify() const
    {
      /// Buffer where to check destination usage
      std::vector<BufComp> inBufOfDest(getNDst()(),-1);
      
      for(BufComp i=0;i<getInBufSize();i++)
	{
	  auto& p=inBufOfDest[dstOfInBuf[i]()];
	  if(p!=-1)
	    crash("On rank %ld destination %ld is filled by %ld and at least %ld at the same time",
		  thisRank(),(int64_t)dstOfInBuf[i](),(int64_t)p(),(int64_t)i());
	  p=i;
	}
      
      /// Buffer where to check source usage
      std::vector<CSrc> srcOfOutBuf(getNSrc()(),-1);
      
      for(CSrc i=0;i<getNSrc();i++)
	{
	  auto& p=srcOfOutBuf[outBufOfSrc[i]()];
	  if(p!=-1)
	    crash("On rank %ld source %ld is filling %ld and at least %ld at the same time",
		  thisRank(),(int64_t)outBufOfSrc[i](),(int64_t)p(),(int64_t)i());
	  p=i;
	}
      
    }
    
    /// Initializes the AllToAll communicator
    template <typename F>
    void init(const CSrc& nSrc,
	      F&& f)
    {
      if(inited)
	crash("Cannot init twice");
      inited=true;
      
      outBufOfSrc.allocate(nSrc);
      
      /// For each destination rank, list all local sources
      std::vector<std::vector<CSrc>> locSrcsGroupedByDstRank(nRanks());
      
      /// For each destination rank, list all remote destinations
      std::vector<std::vector<CDst>> remDstsGroupedByDstRank(nRanks());
      for(CSrc locSrc=0;locSrc<nSrc;locSrc++)
	{
	  const auto [dstRank,remDst]=f(locSrc);
	  
	  if(dstRank>=nRanks or dstRank<0)
	    crash("destination rank %d does not exist!",dstRank);
	  
	  locSrcsGroupedByDstRank[dstRank()].emplace_back(locSrc);
	  remDstsGroupedByDstRank[dstRank()].emplace_back(remDst);
	}
      
      // for(MpiRank iRank=0;iRank<nranks;iRank++)
      // 	printf("AllToAll Rank %d transmitting to rank %d: %zu el\n",rank,iRank(),locSrcsGroupedByDstRank[iRank()].size());
      
      /// Progressive storage of the dst
      std::vector<CDst> buildDstOfInBuf;
      
      /// Fillable view on the mirrored node
      auto fillableOutBufOfSrc=
	outBufOfSrc.getFillable();
      
      BufComp nOutBuf=0;
      for(MpiRank dRank=0;dRank<nRanks;dRank++)
	{
	  /// Rank towards which to send
	  const MpiRank sendRank=
	    (thisRank+nRanks+dRank)%nRanks;
	  
	  for(const CSrc& locSrc : locSrcsGroupedByDstRank[sendRank()])
	    {
	      // printf("AllToAll Rank %d filling %zu with nOutBuf %zu to be sent to rank %d\n",rank,locSrc(),nOutBuf(),sendRank);
	      
	      fillableOutBufOfSrc[locSrc]=nOutBuf++;
	    }
	  
	  /// Rank from which to receive
	  const MpiRank recvRank=
	    (thisRank+nRanks-dRank)%nRanks;
	  
	  /// List of position in the ouput where to store the portion of the buffer relative to recvRank
	  const std::vector<CDst> dstOfBufFrRank=
	    mpiSendrecv(sendRank,remDstsGroupedByDstRank[sendRank()],recvRank);
	  
	  // printf("AllToAll Rank %d dRank %d sent to rank %d value
	  // nRemDstOfRank=%zu received from rank %d value
	  // %zu\n",rank,dRank,sendRank,remDstsGroupedByDstRank[sendRank].size(),recvRank,dstOfBufFrRank.size());
	  
	  // Increment the list of object to pull out from the buffer
	  for(const CDst& bufDst : dstOfBufFrRank)
	    buildDstOfInBuf.emplace_back(bufDst);
	  
	  // If the node has to receive, mark the node in the list with the size attached
	  if(const size_t s=dstOfBufFrRank.size();s)
	    nRecvFrRank.emplace_back(recvRank,s);
	  
	  // If the node has to send, mark the node to the list with the size attached
	  if(const size_t s=remDstsGroupedByDstRank[sendRank()].size();s)
	    nSendToRank.emplace_back(sendRank,s);
	}
      
      // Assigning the destination of all incoming buffer
      dstOfInBuf.allocate((BufComp)(int64_t)buildDstOfInBuf.size());
      
      /// Gets a fillable view on dstOfInBuf
      auto fillableDstOfInBuf=
	dstOfInBuf.getFillable();
      
      for(BufComp bc=0;const CDst& cDst : buildDstOfInBuf)
	fillableDstOfInBuf[bc++]=cDst;
      
      verify();
    }
    
    /// Communicate the expression
    template <DerivedFromNode _DstExpr,
	      DerivedFromNode SrcExpr>
    void communicate(_DstExpr&& out,
		     const SrcExpr& in) const
    {
      /// Actual type of the destination expression
      using DstExpr=std::decay_t<_DstExpr>;
      
      // Basic tests on presence of the asked components in the expressions
      static_assert(tupleHasType<typename DstExpr::Comps,CDst>,"Destination does not have the destination component");
      static_assert(tupleHasType<typename SrcExpr::Comps,CSrc>,"Source does not have the source component");
      
      /// Determine the execution space of the expression
      constexpr ExecSpace execSpace=
		  SrcExpr::execSpace*DstExpr::execSpace;
      
      static_assert(UniqueExecSpace<execSpace>,"Needs to have same exec space for src and dest");
      
      /// Components in the destination apart from CDSt
      using DstRedComps=
	TupleFilterAllTypes<typename DstExpr::Comps,std::tuple<CDst>>;
      
      /// Components in the source apart from CSrc
      using SrcRedComps=
	TupleFilterAllTypes<typename SrcExpr::Comps,std::tuple<CSrc>>;
      
      // if constexpr(not tupleHaveTypes<DstRedComps,SrcRedComps>)
      // 	{
      // 	  DstRedComps& aaa="let the error";
      // 	  SrcRedComps& bbb="get prompte!d";
      // 	}
      
      static_assert(tupleHaveTypes<DstRedComps,SrcRedComps>,"Components in the source must be the same apart from CDst and CSrc");
      
      /// Using destination Fund
      using Fund=
	DstExpr::Fund;
      
      /// Build the components of the buffer
      using BufComps=
	TupleCat<CompsList<BufComp>,DstRedComps>;
      
      /// Instantiate the buffer
      DynamicTens<BufComps,Fund,getMemoryType<execSpace>()> outBuf(std::tuple_cat(std::make_tuple(getOutBufSize()),in.getDynamicSizes()));
      
      // Set the out buffer to an awkward value to detect possible errors
      // if constexpr(std::is_same_v<Fund,double>)
      // outBuf=-10;
      
      /// Check that the out buffer size has the same length of the source
      const CSrc nSrc=
	in.template getCompSize<CSrc>();
      if(nSrc()!=getOutBufSize()())
	crash("Size of the in expression %ld different from expected %ld\n",(int64_t)nSrc(),(int64_t)getOutBufSize()());
      
      // Fills the output buffer
      PAR_ON_EXEC_SPACE(execSpace,
			0,
			nSrc,
			CAPTURE(TO_READ(in),
				outBufOfSrc=outBufOfSrc.template getRefForExecSpace<execSpace>(),
				TO_WRITE(outBuf)),
			iIn,
			{
			  outBuf(outBufOfSrc(iIn))=in(iIn);
			});
      
      /// Copy the output buffer to CPU, if needed
      decltype(auto) mergedHostOutBuf=
	outBuf.template copyToMemorySpaceIfNeeded<MemoryType::CPU>().template mergeComps<DstRedComps>();
      
      /// Extra components of the host copy of the output buffer
      using MergedHostOutBufExtraComps=
	TupleFilterAllTypes<typename std::decay_t<decltype(mergedHostOutBuf)>::Comps,CompsList<BufComp>>;
      
      /// Allocate input buffer
      DynamicTens<BufComps,Fund,MemoryType::CPU> hostInBuf(std::tuple_cat(std::make_tuple(getInBufSize()),out.getDynamicSizes()));
      
      // Set the in buffer to an awkward value to detect possible errors
      // if constexpr(std::is_same_v<Fund,double>)
      // hostInBuf=-9;
      
      /// Gets a merged view of the input buffer
      decltype(auto) mergedHostInBuf=
	hostInBuf.template mergeComps<DstRedComps>();
      
      /// Extra components of the merged input buffer
      using MergedHostInBufExtraComps=
	TupleFilterAllTypes<typename std::decay_t<decltype(mergedHostInBuf)>::Comps,CompsList<BufComp>>;
      
      /// Number of requests to be spawned by MPI
      const int nReqs=
	nSendToRank.size()+nRecvFrRank.size();
      
      /// Requests to be spawned
      MPI_Request reqs[nReqs];
      
      /// Number of degrees of freedom for extra components
      const auto nDof=
	mergedHostOutBuf.template getCompSize<MergedComp<DstRedComps>>()();
      
      /// Request under process
      MPI_Request* req=reqs;
      
      // Spawn the send requests
      for(BufComp sendOffset=0;
	  const auto& [dstRank,nEl] : nSendToRank)
	{
	  /// Pointer to location where the chunk for dstRank starts
	  const Fund& ptr=
	    invokeWithTypesOfTuple<MergedHostOutBufExtraComps>([&mergedHostOutBuf,
								&sendOffset]<DerivedFromComp...C>()->const Fund&
							       {
								 return mergedHostOutBuf(sendOffset,C(0)...);
							       });
	  //master_printf("Sending nel %d ndof %d %d bytes to rank %d\n",nEl,nDof,nEl*nDof*sizeof(Fund),dstRank());
	  MPI_Isend(&ptr,
		    nEl*nDof*sizeof(Fund),
		    MPI_CHAR,dstRank(),0,MPI_COMM_WORLD,req++);
	  sendOffset+=nEl;
	}
      
      for(BufComp recvOffset=0;
	  const auto& [rcvRank,nEl] : nRecvFrRank)
	{
	  /// Pointer to location where the chunk from recvRank starts
	  Fund& ptr=
	    invokeWithTypesOfTuple<MergedHostInBufExtraComps>([&mergedHostInBuf,
							       &recvOffset]<DerivedFromComp...C>()->Fund&
							      {
								return mergedHostInBuf(recvOffset,C(0)...);
							      });
	  //master_printf("Receiving %d bytes from rank %d ptr %p\n",nEl*nDof*sizeof(Fund),rcvRank(),&ptr);
	  MPI_Irecv(&ptr,
		    nEl*nDof*sizeof(Fund),
		    MPI_CHAR,rcvRank(),0,MPI_COMM_WORLD,req++);
	  recvOffset+=nEl;
	}
      
      // Waits all communications
      MPI_Waitall(nReqs,reqs,MPI_STATUS_IGNORE);
      
      /// Ensure that the incoming buffer can be accessed in the device
      decltype(auto) inBuf=
	hostInBuf.template copyToMemorySpaceIfNeeded<getMemoryType<execSpace>()>();
      
      // Fills the destination
      PAR_ON_EXEC_SPACE(execSpace,
			0,
			getInBufSize(),
			CAPTURE(TO_READ(inBuf),
				dstOfInBuf=dstOfInBuf.template getRefForExecSpace<execSpace>(),
				TO_WRITE(out)),
			iInBuf,
			{
			  out(dstOfInBuf(iInBuf))=inBuf(iInBuf);
			});
    }
    
    /// Communicate the expression, building the result
    template <DerivedFromNode SrcExpr,
	      DerivedFromNode DestExpr=DynamicTens<TupleReplaceType<typename SrcExpr::Comps,CSrc,CDst>,typename SrcExpr::Fund,getMemoryType<SrcExpr::execSpace>()>>
    DestExpr communicate(const SrcExpr& in) const
    {
      /// Create the result
      /// \todo possibly will fail if CSrc is not dynamic
      DestExpr out(std::tuple_cat(tupleFilterAllTypes<CompsList<CSrc>>(in.getDynamicSizes()),std::make_tuple(getNDst())));
      
      communicate(out,in);
      
      return out;
    }
    
    /// Type of the inverse communicator
    using Inverse=AllToAllComm<CSrc,CDst>;
    
    /// Returns an inverse communicator
    AllToAllComm<CSrc,CDst> inverse() const
    {
      /// Result
      Inverse res;
      
      /// Number of elments in the current source
      const CSrc nThisSrc=
	outBufOfSrc.template getCompSize<CSrc>();
      
      res.dstOfInBuf.allocate(nThisSrc.template castTo<typename Inverse::BufComp>());
      
      /// Gets fillable version of res.dstOfInBuf
      auto resDstOfInBuf=
	res.dstOfInBuf.getFillable();
      
      for(CSrc src=0;src<nThisSrc;src++)
       	resDstOfInBuf(outBufOfSrc(src).template castTo<typename Inverse::BufComp>())=src;
      
      /// Number of elements in the current destination
      const BufComp nThisBufOut=
	dstOfInBuf.template getCompSize<BufComp>();
      
      res.outBufOfSrc.allocate(nThisBufOut.template castTo<CDst>());
      
      /// Gets fillable version of res.outBufOsSrc
      auto resOutBufOfSrc=
	res.outBufOfSrc.getFillable();
      
      for(BufComp outBuf=0;outBuf<nThisBufOut;outBuf++)
	resOutBufOfSrc(dstOfInBuf(outBuf).template castTo<CDst>())=outBuf.template castTo<typename Inverse::BufComp>();
      
      res.nSendToRank=nRecvFrRank;
      res.nRecvFrRank=nSendToRank;
      
      res.verify();
      
      res.inited=true;
      
      return res;
    }
  };
  
  /// Product of two communicators
  ///
  /// The trick is to communicate the original rank and position
  /// twice, us it to define the inverse of the product of the
  /// application, and invert it
  template <DerivedFromComp CDst,
	    DerivedFromComp CTmp,
	    DerivedFromComp CSrc>
  auto operator*(const AllToAllComm<CDst,CTmp>& second,
		 const AllToAllComm<CTmp,CSrc>& first)
  {
    /// Number of elements in the source
    const CSrc nSrc=
      first.getNSrc();
    
    /// Needed to follow the application of first and second operator forward
    DynamicTens<CompsList<CSrc>,std::tuple<MpiRank,CSrc>,MemoryType::CPU> identity(nRanks,nSrc);
    
    for(CSrc src=0;src<nSrc;src++)
      identity(src)=std::make_tuple(thisRank,src);
    
    /// Gets the initial position of the ultimate destination of the product
    const auto initSrc=
      second.communicate(first.communicate(identity));
    
    /// Fills the inverse communicator
    AllToAllComm<CSrc,CDst> inverseRes(initSrc.template getCompSize<CSrc>(),
				       [&initSrc](const CDst& inverseDest)
    {
      return initSrc(inverseDest);
    });
    
    inverseRes.verify();
    
    return inverseRes.inverse();
  }
}

#endif
