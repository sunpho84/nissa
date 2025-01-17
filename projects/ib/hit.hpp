#ifndef _HIT_HPP
#define _HIT_HPP

#include "prop.hpp"
#include "contr.hpp"

using namespace nissa;

struct HitLooper
{
  /// Store all the dependencies
  std::map<std::string,std::set<std::string>> propDep;
  
  std::vector<std::string> oldDepInserted;
  
  /// List of computed propagators
  std::set<std::string> computedPropsList;
  
  /// List of computed 2pts
  std::set<int> computed2ptsList;
  
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
	  ReadWriteRealVector<spincolor> rwTest(q[isou],path);
	  
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
  void offloadIfNeededToAccommodate(const int nToAccommodate)
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
    FieldRngOf<spincolor> drawer(field_rng_stream.getDrawer<spincolor>());
    if(stoch_source and skipOnly)
      return;
    
    //consistency check
    if(not stoch_source and (not diluted_spi_source or not diluted_col_source))
      CRASH("for a non-stochastic source, spin and color must be diluted");
    
    PAR(0,locVol,
	CAPTURE(tins=sou->tins,
		drawer,
		noise_type=sou->noise_type,
		s=sou->sp[0]->getWritable()),
	ivol,
	{
	  spincolor b{};
	  
	  if(stoch_source)
	    {
	      if(tins==-1 or rel_coord_of_loclx(ivol,0)==0)
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
		rel_coord_of_loclx(ivol,mu);
	      
	      if(is_ori)
		b[0][0][RE]=1.0;
	    }
	  
	  spincolor_copy(s[ivol],b);
	});
    
    for(int i=1;i<nso_spi*nso_col;i++)
      {
	PAR(0,locVol,
	  CAPTURE(s=sou->sp[0]->getReadable(),
		  d=sou->sp[i]->getWritable(),
		  i),
	  ivol,
	  {
	    const auto [sp_so,co_so]=sp_col_of_so_ind(i);
	    
	    const int sp_si=diluted_spi_source?0:sp_so;
	    const int co_si=diluted_col_source?0:co_so;
	    complex_copy(d[ivol][sp_so][co_so],s[ivol][sp_si][co_si]);
	  });
      }
    
    //compute the norm2, set borders invalid
    double ori_source_norm2=0;
    for(int id_so=0;id_so<nso_spi;id_so++)
      for(int ic_so=0;ic_so<nso_col;ic_so++)
	{
	  LxField<spincolor>& s=*sou->sp[so_sp_col_ind(id_so,ic_so)];
	  s.invalidateHalo();
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
  void generate_original_sources(int ihit,bool skipOnly)
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
		
		LxField<spincolor>& sou=
		  (*q)[isou];
		
		ReadWriteRealVector<spincolor> rw(sou,path);
		
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
  
  void start_hit(int ihit,bool skip=false)
  {
    MASTER_PRINTF("\n=== Hit %d/%d ====\n",ihit+1,nhits);
    
    if(doNotAverageHits)
      clearCorrelations();
    
    for(int mu=0;mu<NDIM;mu++)
      {
	using C=double[1];
	C c;
	field_rng_stream.drawScalar(c);
	source_coord[mu]=c[0]*glbSize[mu];
      }
    
    if(stoch_source) MASTER_PRINTF(" source time: %d\n",source_coord[0]);
    else             MASTER_PRINTF(" point source coords: %d %d %d %d\n",source_coord[0],source_coord[1],source_coord[2],source_coord[3]);
    if(need_photon)
      {
	if(skip)
	  generate_photon_source(*photon_eta);
	else
	  generate_photon_stochastic_propagator(ihit);
      }
    generate_original_sources(ihit,skip);
  }
  
  /// Perform one of step: dryRun, or eval
  void internalRun(const int runStep,const size_t iHit)
  {
    computedPropsList.clear();
    computed2ptsList.clear();
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
	  
	  for(size_t icombo=0;icombo<mes2pts_contr_map.size();icombo++)
	    if(const std::string& a=mes2pts_contr_map[icombo].a,
	       b=mes2pts_contr_map[icombo].b;
	       computedPropsList.count(a) and
	       computedPropsList.count(b) and
	       computed2ptsList.count(icombo)==0)
	      {
		if(runStep==2)
		  MASTER_PRINTF("Can compute contraction: %s %s -> %s, %zu/%zu\n",
			      a.c_str(),b.c_str(),mes2pts_contr_map[icombo].name.c_str(),computed2ptsList.size(),mes2pts_contr_map.size());
		
		// Insert the correlation in the computed list
		computed2ptsList.insert(icombo);
		
		for(const std::string& p : {a,b})
		  if(runStep==1)
		    orderedDep.push_back(p);
		  else
		    ensureInMemory(p);
		
		// Compute the correlation function
		if(runStep==2)
		  compute_mes2pt_contr(icombo);
		
		// Remove the dependency from the 2pts of the props
		for(const std::string& p : {a,b})
		  propDep[p].erase("2pts_"+std::to_string(icombo));
		
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
  }
  
  void run(const size_t iHit)
  {
    for(int runStep=1;runStep<=2;runStep++)
      internalRun(runStep,iHit);
    
    std::ostringstream os;
    for(size_t iContr=0;iContr<mes2pts_contr_map.size();iContr++)
      if(const auto& m=mes2pts_contr_map[iContr];computed2ptsList.find(iContr)==computed2ptsList.end())
	os<<" not computed corr "<<m.name<<" between "<<m.a<<" and "<<m.b<<std::endl;
    
    if(os.str().size())
      CRASH("%s",os.str().c_str());
  }
  
  /// Constructor
  HitLooper()
  {
    // Insert all contractions
    for(size_t icombo=0;icombo<mes2pts_contr_map.size();icombo++)
      for(const std::string& p : {mes2pts_contr_map[icombo].a,mes2pts_contr_map[icombo].b})
	{
	  propDep[p].insert("2pts_"+std::to_string(icombo));
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

#endif
