#ifndef _HIT_HPP
#define _HIT_HPP

#include "prop.hpp"
#include "contr.hpp"
#include "routines/ios.hpp"

using namespace nissa;

namespace nissa
{
  struct HitLooper
  {
    /// Store all the dependencies
    std::map<std::string,std::set<std::string>> propDep;
    
    std::vector<std::string> oldDepInserted;
    
    /// List of computed propagators
    std::set<std::string> computedPropsList;
    
    /// List of computed 2pts
    std::set<int> computed2ptsList;
    
    /// List of computed handcuff sides
    std::set<std::string> computedHandcuffSidesList;
    
    /// Ordered dependecies
    std::vector<std::string> orderedDep;
    
    /// Number of passed dependencies
    int nPassedDep=0;
    
    /// Number of propagators which have been offloaded
    int nOffloaded=0;
    
    /// Number of propagators which have been recalled
    int nRecalled=0;
    
    /// Propagator status
    enum Status{NOT_COMPUTED,IN_MEMORY,OFFLOADED};
    
    /// Status of all propagators
    std::map<std::string,Status> propagatorsStatus;
    
    /// List of props which have been offloaded
    std::set<std::string> offloadedList;
    
    enum ORD{OFFLOAD,RECALL,DELETE};
    
    /// Offload, recall or delete individual propagators
    void offloadRecallDelete(const std::string& name,const ORD ord)
    {
      qprop_t& q=Q[name];
      
      if(ord==RECALL)
	{
	  if(propagatorsStatus[name]!=OFFLOADED)
	    CRASH("Asking to recall something not offloaded");
	  
	  offloadIfNeededToAccommodate(1);
	  q.alloc_storage();
	  propagatorsStatus[name]=IN_MEMORY;
	  nRecalled++;
	}
      
      double rwBeg=take_time();
      for(int id_so=0;id_so<nso_spi;id_so++)
	for(int ic_so=0;ic_so<nso_col;ic_so++)
	  {
	    int isou=so_sp_col_ind(id_so,ic_so);
	    const std::string path=combine("%s/prop%s_idso%d_icso%d",outfolder,name.c_str(),id_so,ic_so);
	    ReadWriteRealVector<spincolor,defaultSpaceTimeLayout,MemorySpace::CPU> rwTest(q[isou],path);
	    
	    switch(ord)
	      {
	      case OFFLOAD:
		rwTest.fastWrite();
		break;
	      case RECALL:
		rwTest.fastRead();
		break;
	      case DELETE:
		rwTest.cleanFiles();
		break;
	      }
	  }
      const char tag[3][15]={"Offloading","Recalling","Erasing"};
      MASTER_PRINTF("%s took: %lg s\n",tag[ord],take_time()-rwBeg);
      
      if(ord==OFFLOAD)
	{
	  if(propagatorsStatus[name]!=IN_MEMORY)
	    CRASH("Asking to offload something not in memory");
	  
	  const int64_t pre=required_memory;
	  q.free_storage();
	  const int64_t aft=required_memory;
	  MASTER_PRINTF("Freed after offloading, memory before: %ld bytes, after: %ld bytes\n",pre,aft);
	  propagatorsStatus[name]=OFFLOADED;
	  nOffloaded++;
	}
    }
    
    /// Erase all the propagators which are not needed
    void eraseUnnededProps()
    {
      MASTER_PRINTF("Checking what can be erased\n");
      
      /// List of propagators to be erased
      std::set<std::string> toErase;
      for(const auto& [n,s] : propagatorsStatus)
	{
	  if(propDep[n].empty())
	    {
	      MASTER_PRINTF("Erasing %s\n",n.c_str());
	      toErase.insert(n);
	    }
	  else
	    VERBOSITY_LV3_MASTER_PRINTF("keeping %s as it has %zu deps, first is: %s\n",n.c_str(),propDep[n].size(),propDep[n].begin()->c_str());
	}
      
      for(const auto& e : toErase)
	{
	  propagatorsStatus.erase(e);
	  Q[e].free_storage();
	  MASTER_PRINTF("%s freed\n",e.c_str());
	  
	  // Phyiscally remove from disk
	  if(auto o=offloadedList.find(e);o!=offloadedList.end())
	    {
	      MASTER_PRINTF("%s can be deleted from disk\n",e.c_str());
	      offloadedList.erase(o);
	      offloadRecallDelete(e,DELETE);
	    }
	}
      
      MASTER_PRINTF("n Live props: %zu\n",propagatorsStatus.size());
      MASTER_PRINTF("n offloaded props: %zu\n",offloadedList.size());
    }
    
    /// Offload some propagators, to accommodate nToAccommodate more
    void offloadIfNeededToAccommodate(const int& nToAccommodate)
    {
      if(nMaxPropsAllocated==0)
	return;
      
      /// Build the list of future dependencies
      std::vector<std::pair<int,std::string>> futureDeps;
      for(const auto& [s,n] : propagatorsStatus)
	if(n==1)
	  {
	    size_t firstUsage=nPassedDep;
	    while(firstUsage<orderedDep.size() and orderedDep[firstUsage]!=s)
	      firstUsage++;
	    
	    if(firstUsage==orderedDep.size())
	      CRASH("At nPassedDep %d unable to find the first usage for %s, why is it allocated?",nPassedDep,s.c_str());
	    else
	      futureDeps.emplace_back(firstUsage-nPassedDep,s);
	  }
      
      VERBOSITY_LV2_MASTER_PRINTF("futureDeps:\n");
      for([[maybe_unused]] const auto&  [delay,name] : futureDeps)
	{
	  VERBOSITY_LV2_MASTER_PRINTF(" %d %s\n",delay,name.c_str());
	}
      
      // \todo simplify
      
      int nLoaded=0;
      for(const auto& s : propagatorsStatus)
	nLoaded+=(s.second==1);
      VERBOSITY_LV2_MASTER_PRINTF("n Loaded props: %d\n",nLoaded);
      
      if((int)futureDeps.size()!=nLoaded)
	CRASH("unmatched loaded and future deps, %d %zu",nLoaded,futureDeps.size());
      
      MASTER_PRINTF("Needs to accommodate %d more props, %d are loaded and max is %d, needs to offload %d\n",
		    nToAccommodate,nLoaded,nMaxPropsAllocated,nLoaded+nToAccommodate-nMaxPropsAllocated);
      
      std::sort(futureDeps.begin(),futureDeps.end());
      for(size_t i=nMaxPropsAllocated-nToAccommodate;i<futureDeps.size();i++)
	{
	  const auto& [delay,name]=futureDeps[i];
	  if(delay==0)
	    CRASH("Unable to offload %s which will be needed immediately!",name.c_str());
	  MASTER_PRINTF("Offloading %s which will be used in %d\n",name.c_str(),delay);
	  offloadRecallDelete(name,OFFLOAD);
	  offloadedList.insert(name);
	}
    }
    
    /// Ensure that a prop is in memory
    void ensureInMemory(const std::string& name)
    {
      if(propagatorsStatus[name]==OFFLOADED)
	{
	  MASTER_PRINTF("Recalling %s\n",name.c_str());
	  offloadRecallDelete(name,RECALL);
	}
    }
    
    /// Generate a source, either a wall or a point in the origin
    void generate_original_source(qprop_t* sou,
				  const bool& skipOnly)
    {
      VERBOSITY_LV3_MASTER_PRINTF("creating the drawer\n");
      FieldRngOf<spincolor> drawer(field_rng_stream.getDrawer<spincolor>());
      if(stoch_source and skipOnly)
	{
	  VERBOSITY_LV3_MASTER_PRINTF("Going to destroy the drawer\n");
	  return;
	}
      
      //consistency check
      if(not stoch_source and (not diluted_spi_source or not diluted_col_source))
	CRASH("for a non-stochastic source, spin and color must be diluted");
      
      sou->sp[0]->initOn<defaultMemorySpace>([drawer,
					      sou](LxField<spincolor>& s)
      {
	PAR(0,locVol,
	    CAPTURE(tins=sou->tins,
		    drawer,
		    noise_type=sou->noise_type,
		    s=s.getWritable()),
	    ivol,
	    {
	      spincolor b{};
	      
	      if(stoch_source)
		{
		  if(tins==-1 or rel_coord_of_loclx(ivol,0)==tins)
		    {
		      spincolor drawn;
		      drawer.fillLocSite(drawn,ivol);
		      
		      for(int id=0;id<(diluted_spi_source?1:NDIRAC);id++)
			for(int ic=0;ic<(diluted_col_source?1:NCOL);ic++)
			  {
			    complex bi;
			    
			    complex_copy(bi,drawn[id][ic]);
			    
			    switch(noise_type)
			      {
			      case RND_ALL_PLUS_ONE:
				complex_put_to_real(bi,1);
				break;
			      case RND_ALL_MINUS_ONE:
				complex_put_to_real(bi,-1);
				break;
			      case RND_Z2:
				z2Transform(bi);
				break;
			      case RND_Z4:
				z4Transform(bi);
				break;
			      case RND_UNIF:
			      case RND_Z3:
			      case RND_GAUSS:
				CRASH("not implemented yet");
				break;
			      }
			    
			    complex_copy(b[id][ic],bi);
			  }
		    }
		}
	      else
		{
		  bool is_ori=true;
		  for(int mu=0;mu<NDIM;mu++)
		    is_ori&=(rel_coord_of_loclx(ivol,mu)==0);
		  
		  if(is_ori)
		    b[0][0][RE]=1.0;
		}
	      
	      spincolor_copy(s[ivol],b);
	    });
      });
      
      for(int i=1;i<nso_spi*nso_col;i++)
	{
	  const auto sc_so=sp_col_of_so_ind(i);
	  const int sp_so=sc_so.first,co_so=sc_so.second;
	  MASTER_PRINTF("Diluting on index %d={%d,%d}\n",i,sp_so,co_so);
	  
	  sou->sp[i]->initOn<defaultMemorySpace>([sp_so,co_so,&sou](LxField<spincolor>& d)
	  {
	    decltype(auto) s=sou->sp[0]->getSurelyReadableOn<defaultMemorySpace>();
	    
	    PAR(0,locVol,
		CAPTURE(TO_READ(s),
			TO_WRITE(d),
			sp_so,
			co_so),
		ivol,
		{
		  UNROLL_FOR_ALL_SPIN(sp_si)
		    UNROLL_FOR_ALL_COLS(co_si)
		    {
		      if((sp_si==sp_so or not diluted_spi_source) and (co_si==co_so or not diluted_col_source))
			complex_copy(d[ivol][sp_si][co_si],s[ivol][diluted_spi_source?0:sp_si][diluted_col_source?0:co_si]);
		      else
			complex_put_to_zero(d[ivol][sp_si][co_si]);
		    }
		});
	  });
	}
      
      //compute the norm2, set borders invalid
      double ori_source_norm2=0;
      for(int id_so=0;id_so<nso_spi;id_so++)
	for(int ic_so=0;ic_so<nso_col;ic_so++)
	  {
	    decltype(auto) s=
	      sou->sp[so_sp_col_ind(id_so,ic_so)]->getSurelyReadableOn<defaultMemorySpace>();
	    ori_source_norm2+=s.norm2();
	  }
      sou->ori_source_norm2=ori_source_norm2;
      
      // complex *n=nissa_malloc("n",locVol,complex);
      // spincolor *temp=nissa_malloc("temp",locVol+bord_vol,spincolor);
      // for(int id_so=0;id_so<nso_spi;id_so++)
      //   for(int ic_so=0;ic_so<nso_col;ic_so++)
      // 	{
      // 	  spincolor *s=sou->sp[so_sp_col_ind(id_so,ic_so)];
      // 	  MASTER_PRINTF("eta (0): %lg %lg\n",s[0][0][0][RE],s[0][0][0][IM]);
      // 	  fft4d((complex*)temp,(complex*)s,NDIRAC*NCOL,FFT_PLUS,FFT_NO_NORMALIZE);
      
      // 	  NISSA_PARALLEL_LOOP(ivol,0,locVol)
      // 	    {
      // 	      n[ivol][RE]=spincolor_norm2(temp[ivol]);
      // 	      n[ivol][IM]=0;
      // 	    }
      // 	  NISSA_PARALLEL_LOOP_END;
      
      // 	  fft4d(n,n,1,FFT_MINUS,FFT_NORMALIZE);
      
      // 	  const int64_t nPoints=((tins==-1)?glbVol:glbSpatVol);
      
      // 	  MASTER_PRINTF("eta+ eta (0): %lg , eta+ eta (1): %lg\n",n[0][RE],n[1][RE]);
      // 	  if(rank==0)
      // 	    n[0][RE]-=(double)NDIRAC*NCOL*nPoints/nso_spi/nso_col;
      // 	  MASTER_PRINTF("eta+ eta (0) after sub: %lg\n",n[0][RE]);
      
      // 	  NISSA_PARALLEL_LOOP(ivol,0,locVol)
      // 	    {
      // 	      n[ivol][RE]*=n[ivol][RE];
      // 	    }
      // 	  NISSA_PARALLEL_LOOP_END;
      
      // 	  complex res[1];
      // 	  glb_reduce(res,n,locVol);
      
      // 	  MASTER_PRINTF("Res: %lg\n",res[0][RE]);
      
      // 	  double exp=(double)NDIRAC*NCOL*sqr(nPoints)/(nso_spi*nso_col);
      // 	  MASTER_PRINTF("Exp: %lg\n",exp);
      
      // 	  double dev=res[0][RE]/exp-1;
      
      // 	  MASTER_PRINTF("Dev: %lg\n",dev);
      // 	}
      
      // nissa_free(temp);
      // nissa_free(n);
    }
    
    //Generate all the original sources
    void generate_original_sources(const int& ihit,
				   const bool& skipOnly=false)
    {
      for(size_t i=0;i<ori_source_name_list.size();i++)
	{
	  const std::string &name=
	    ori_source_name_list[i];
	  MASTER_PRINTF("Generating source \"%s\"\n",name.c_str());
	  
	  qprop_t *q=&Q[name];
	  q->alloc_storage();
	  generate_original_source(q,skipOnly);
	  
	  if(not skipOnly)
	    for(int id_so=0;id_so<nso_spi;id_so++)
	      for(int ic_so=0;ic_so<nso_col;ic_so++)
		{
		  //combine the filename
		  const std::string path=
		    combine("%s/hit%d_source%s_idso%d_icso%d",outfolder,ihit,name.c_str(),id_so,ic_so);
		  
		  const int isou=
		    so_sp_col_ind(id_so,ic_so);
		  
		  HostProp& sou=
		    (*q)[isou];
		  
		  ReadWriteRealVector<spincolor,defaultSpaceTimeLayout,MemorySpace::CPU> rw(sou,path);
		  
		  //if the prop exists read it
		  if(rw.canLoad())
		    {
		      MASTER_PRINTF("  loading the source dirac index %d, color %d\n",id_so,ic_so);
		      START_TIMING(read_prop_time,nread_prop);
		      rw.read();
		      STOP_TIMING(read_prop_time);
		    }
		  else
		    {
		      MASTER_PRINTF("  file %s not available, skipping loading\n",path.c_str());
		      
		      //and store if needed
		      if(q->store)
			{
			  MASTER_PRINTF("  writing the source dirac index %d, color %d\n",id_so,ic_so);
			  START_TIMING(store_prop_time,nstore_prop);
			  rw.write();
			  STOP_TIMING(store_prop_time);
			}
		    }
		}
	}
    }
    
    void start_hit(const int& iHit,
		   const bool& skip=false)
    {
      MASTER_PRINTF("\n=== Hit %d/%d ====\n",iHit+1,nHits);
      
      if(doNotAverageHits)
	clearCorrelations();
      
      for(int mu=0;mu<NDIM;mu++)
	{
	  using C=double[1];
	  C c;
	  field_rng_stream.drawScalar(c);
	  
	  if(iHit==0 or not doNotShiftSourceOfHits)
	    oriCoords[mu]=c[0]*glbSize[mu];
	}
      
      if(stoch_source)
	MASTER_PRINTF(" source time: %d\n",oriCoords[0]);
      else
	MASTER_PRINTF(" point source coords: %d %d %d %d\n",oriCoords[0],oriCoords[1],oriCoords[2],oriCoords[3]);
      
      if(need_photon)
	{
	  if(skip)
	    generate_photon_source(*photon_eta);
	  else
	    generate_photon_stochastic_propagator(iHit);
	}
      
      generate_original_sources(iHit,skip);
    }
    
    /// Perform one of step: dryRun, or eval
    void internalRun(const int& runStep,
		     const size_t& iHit)
    {
      computedPropsList.clear();
      computed2ptsList.clear();
      computedHandcuffSidesList.clear();
      offloadedList.clear();
      propagatorsStatus.clear();
      
      /// Keep backup of the propagators dependencies
      const std::map<std::string,std::set<std::string>> propDepBack=propDep;
      
      // Mark all the original sources as in memory and computed
      if(runStep==2)
	for(const auto& s : ori_source_name_list)
	  {
	    propagatorsStatus[s]=IN_MEMORY;
	    computedPropsList.insert(s);
	  }
      
      for(size_t i=0;i<qprop_name_list.size();i++)
	if(propDep[qprop_name_list[i]].size()==0)
	  MASTER_PRINTF("Skipping generation of prop %s as it has no dependencies\n",qprop_name_list[i].c_str());
	else
	  {
	    //get names
	    std::string name=qprop_name_list[i];
	    if(runStep==2)
	      MASTER_PRINTF("Computing prop %s\n",name.c_str());
	    
	    qprop_t &q=Q[name];
	    
	    for(const auto& [name,w] : Q[name].source_terms)
	      if(runStep==1)
		orderedDep.push_back(name);
	      else
		ensureInMemory(name);
	    
	    if(runStep==2)
	      {
		offloadIfNeededToAccommodate(1);
		generate_quark_propagator(name,q,iHit);
		propagatorsStatus[name]=IN_MEMORY;
	      }
	    
	    computedPropsList.insert(name);
	    
	    // Remove the dependencies from all sources
	    for(const auto& [s,w] : q.source_terms)
	      if(propDep[s].erase(name)!=1)
		CRASH("unable to remove the dependency %s of %s!",name.c_str(),s.c_str());
	      else
		VERBOSITY_LV2_MASTER_PRINTF("Erased dependency of %s from %s which has still %zu dependencies\n",name.c_str(),s.c_str(),propDep[s].size());
	    
	    if(runStep==2)
	      {
		nPassedDep+=Q[name].source_terms.size();
		eraseUnnededProps();
	      }
	    
	    /// Reorder to capture those that can recycle the props when needed to move to device
	    std::vector<size_t> canBeComputed2pts;
	    for(size_t icombo=0;icombo<mes2ptsContr.size();icombo++)
	      if(const std::string& a=mes2ptsContr[icombo].a,
		 b=mes2ptsContr[icombo].b;
		 computedPropsList.count(a) and
		 computedPropsList.count(b) and
		 computed2ptsList.count(icombo)==0)
		canBeComputed2pts.emplace_back(icombo);
	    
	    std::sort(canBeComputed2pts.begin(),
		      canBeComputed2pts.end(),
		      [](const size_t& i1,
			 const size_t& i2)
		      {
			auto getKey=
			  [](const size_t& i)
			  {
			    auto& m=mes2ptsContr[i];
			    return std::make_tuple(m.a,m.b,i);
			  };
			
			return getKey(i1)<getKey(i2);
		      });
	    
	    for(const size_t& iCombo : canBeComputed2pts)
	      {
		const std::string& a=mes2ptsContr[iCombo].a,
		  b=mes2ptsContr[iCombo].b;
		
		if(runStep==2)
		    MASTER_PRINTF("Can compute 2pts contraction: %s %s -> %s, %zu/%zu\n",
				  a.c_str(),b.c_str(),mes2ptsContr[iCombo].name.c_str(),computed2ptsList.size(),mes2ptsContr.size());
		
		  // Insert the correlation in the computed list
		  computed2ptsList.insert(iCombo);
		  
		  for(const std::string& p : {a,b})
		    if(runStep==1)
		      orderedDep.push_back(p);
		    else
		      ensureInMemory(p);
		  
		  // Compute the correlation function
		  if(runStep==2)
		    compute_mes2pt_contr(iCombo);
		  
		  // Remove the dependency from the 2pts of the props
		  for(const std::string& p : {a,b})
		    propDep[p].erase("2pts_"+std::to_string(iCombo));
		  
		  if(runStep==2)
		    {
		      eraseUnnededProps();
		      nPassedDep+=2;
		    }
		}
	    
	    MASTER_PRINTF("Removing the props from the 2pts contr list\n");
	    std::vector<std::string> toErase;
	    for(auto& [n,v] : mes2ptsPropsLib)
	      toErase.push_back(n);
	     for(const std::string& n : toErase)
	      removeMes2PtsProp(n);
	    
	    for(auto& [name,hs] : handcuffsSides)
	      if(const std::string& bw=hs.bw,
		 fw=hs.fw;
		 computedPropsList.count(bw) and
		 computedPropsList.count(fw) and
		 computedHandcuffSidesList.count(name)==0)
		{
		  if(runStep==2)
		    MASTER_PRINTF("Can compute handcuff contraction side: %s %s -> %s, %zu/%zu\n",
				  bw.c_str(),fw.c_str(),name.c_str(),computedHandcuffSidesList.size(),handcuffsSides.size());
		  
		  // Insert the correlation in the computed list
		  computedHandcuffSidesList.insert(name);
		  
		  for(const std::string& p : {bw,fw})
		    if(runStep==1)
		      orderedDep.push_back(p);
		    else
		      ensureInMemory(p);
		  
		  // Compute the correlation function
		  if(runStep==2)
		    {
		      computeHandcuffSide(hs);
		      if(hs.store==1)
			write_real_vector(combine("%s/hit%zu_handcuff_side_%s",outfolder,iHit,name.c_str()),hs.data,"handcuff");
		    }
		  
		  // Remove the dependency from the 2pts of the props
		  for(const std::string& p : {bw,fw})
		    propDep[p].erase("handcuffSide_"+name);
		  
		  if(runStep==2)
		    {
		      eraseUnnededProps();
		      nPassedDep+=2;
		    }
		}
	  }
      
      propDep=propDepBack;
      
      if(runStep==1)
	{
	  for([[maybe_unused]] const auto& s : orderedDep)
	    VERBOSITY_LV2_MASTER_PRINTF("  %s\n",s.c_str());
	  VERBOSITY_LV2_MASTER_PRINTF("Dependencies end\n");
	}
      else
	computeHandcuffsContr(false);
    }
    
    /// Run the simulation
    void run(const size_t& iHit)
    {
      for(int runStep=1;runStep<=2;runStep++)
	internalRun(runStep,iHit);
      
      std::ostringstream os;
      
      for(size_t iContr=0;iContr<mes2ptsContr.size();iContr++)
	if(const auto& m=mes2ptsContr[iContr];computed2ptsList.find(iContr)==computed2ptsList.end())
	  os<<" not computed corr "<<m.name<<" between "<<m.a<<" and "<<m.b<<std::endl;
      
      for(const auto& [name,handcuffSide] : handcuffsSides)
	if(computedHandcuffSidesList.find(name)==computedHandcuffSidesList.end())
	  os<<" not computed handcuff "<<name<<" between "<<handcuffSide.bw<<" and "<<handcuffSide.fw<<std::endl;
      
      if(os.str().size())
	CRASH("%s",os.str().c_str());
      
      computeHandcuffsContr(true);
    }
    
    /// Load the partial data, if exists
    int maybeLoadPartialData() const
    {
      if(handcuffs.size() or
	 bar2pts_alt_contr_size or
	 bar2pts_contr_size)
	CRASH("reimplement");
      
      int nHitsDone=0;
      if(fileExists(partialDataPath()))
	{
	  FILE* partialFile=open_file(partialDataPath(),"r");
	  if(is_master_rank())
	    {
	      if(fread(&nHitsDone,sizeof(nHitsDone),1,partialFile)!=1)
		CRASH("Failed to load nHitsDone");
	      
	      for(mes_contr_t& m : mes2ptsContr)
		if(fread(m.contr,sizeof(complex),m.contrSize,partialFile)!=m.contrSize)
		CRASH("Failed to load 2pts");
	    }
	  MPI_Bcast(&nHitsDone,1,MPI_INT,master_rank,MPI_COMM_WORLD);
	  for(mes_contr_t& m : mes2ptsContr)
	    MPI_Bcast(m.contr,m.contrSize,MPI_DOUBLE_COMPLEX,master_rank,MPI_COMM_WORLD);
	  
	  close_file(partialFile);
	}
      
      return nHitsDone;
    }
    
    /// Writes the partial data
    void writePartialData(const int& nHitsDone) const
    {
      FILE* partialFile=open_file(partialDataPath(),"w");
      if(is_master_rank())
	{
	  if(fwrite(&nHitsDone,sizeof(nHitsDone),1,partialFile)!=1)
	    CRASH("Failed to write nHitsDone");
	  
	  for(mes_contr_t& m : mes2ptsContr)
	    if(fwrite(m.contr,sizeof(complex),m.contrSize,partialFile)!=m.contrSize)
	    CRASH("Failed to write 2pts");
	  close_file(partialFile);
	}
    }
    
    /// Delete the partial data
    void deletePartialData() const
    {
      if(fileExists(partialDataPath()))
	if(is_master_rank())
	  remove(partialDataPath().c_str());
    }
    
    /// Constructor
    HitLooper()
    {
      // Insert all contractions for 2pts
      for(size_t icombo=0;icombo<mes2ptsContr.size();icombo++)
	for(const std::string& p : {mes2ptsContr[icombo].a,mes2ptsContr[icombo].b})
	  {
	    propDep[p].insert("2pts_"+std::to_string(icombo));
	    oldDepInserted.push_back(p);
	  }
      
      // Insert all contractions for handcuffs sides
      for(const auto& [name,handcuffsSide] : handcuffsSides)
	for(const std::string& p : {handcuffsSide.bw,handcuffsSide.fw})
	  {
	    propDep[p].insert("handcuffSide_"+name);
	    oldDepInserted.push_back(p);
	  }
      
      /// Iterate until no more dependencies are found
      while(oldDepInserted.size())
	{
	  std::vector<std::string> newDepInserted;
	  for(const std::string& p : oldDepInserted)
	    for(const auto& [n,w] : Q[p].source_terms)
	      {
		newDepInserted.push_back(n);
		propDep[n].insert(p);
	      }
	  oldDepInserted=newDepInserted;
	}
    }
  };
}

#endif
