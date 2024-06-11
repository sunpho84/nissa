#ifndef _HIT_HPP
#define _HIT_HPP

#include <ib/prop.hpp>
#include <ib/contr.hpp>

using namespace nissa;

struct HitLooper
{
  /// Store all the dependencies
  std::map<std::string,std::set<std::string>> propDep;
  
  std::vector<std::string> oldDepInserted;
  
  /// List of computed propagators
  std::set<std::string> computedProps;
  
  /// List of computed 2pts
  std::set<int> computed2pts;
  
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
  std::map<std::string,Status> status;
  
  /// List of props which have been offloaded
  std::set<std::string> offloadedList;
  
  enum ORD{OFFLOAD,RECALL,DELETE};
  
  /// Offload, recall or delete individual propagators
  void offloadRecallDelete(const std::string& name,const ORD ord)
  {
    qprop_t& q=Q[name];
    
    if(ord==RECALL)
      {
	if(status[name]!=OFFLOADED)
	  crash("Asking to recall something not offloaded");
	
	offloadIfNeededToAccommodate(1);
	q.alloc_storage();
	status[name]=IN_MEMORY;
	nRecalled++;
      }
    
    double rwBeg=take_time();
    for(int id_so=0;id_so<nso_spi;id_so++)
      for(int ic_so=0;ic_so<nso_col;ic_so++)
	{
	  int isou=so_sp_col_ind(id_so,ic_so);
	  spincolor *sol=nullptr;
	  if(ord!=DELETE)
	    sol=q[isou];
	  const std::string path=combine("%s/prop%s_idso%d_icso%d",outfolder,name.c_str(),id_so,ic_so);
	  ReadWriteRealVector<spincolor> rwTest(sol,path,true);
	  
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
    master_printf("%s took: %lg s\n",tag[ord],take_time()-rwBeg);
    
    if(ord==OFFLOAD)
      {
	if(status[name]!=IN_MEMORY)
	  crash("Asking to offload something not in memory");
	
        const int64_t pre=required_memory;
        q.free_storage();
        const int64_t aft=required_memory;
	master_printf("Freed after offloading, memory before: %ld bytes, after: %ld bytes\n",pre,aft);
	status[name]=OFFLOADED;
	nOffloaded++;
      }
  }
  
  /// Erase all the propagators which are not needed
  void eraseUnnededProps()
  {
    master_printf("Checking what can be erased\n");
    
    /// List of propagators to be erased
    std::set<std::string> toErase;
    for(const auto& [n,s] : status)
      {
	if(propDep[n].empty())
	  {
	    master_printf("Erasing %s\n",n.c_str());
	    toErase.insert(n);
	  }
	else
	  verbosity_lv3_master_printf("keeping %s as it has %zu deps, first is: %s\n",n.c_str(),propDep[n].size(),propDep[n].begin()->c_str());
      }
    
    for(const auto& e : toErase)
      {
	status.erase(e);
	Q[e].free_storage();
	master_printf("%s freed\n",e.c_str());
	
	// Phyiscally remove from disk
	if(auto o=offloadedList.find(e);o!=offloadedList.end())
	  {
	    master_printf("%s can be deleted from disk\n",e.c_str());
	    offloadedList.erase(o);
	    offloadRecallDelete(e,DELETE);
	  }
      }
    
    master_printf("n Live props: %zu\n",status.size());
    master_printf("n offloaded props: %zu\n",offloadedList.size());
  }
  
  /// Offload some propagators, to accommodate nToAccommodate more
  void offloadIfNeededToAccommodate(const int nToAccommodate)
  {
    if(nMaxPropsAllocated==0)
      return;
    
    /// Build the list of future dependencies
    std::vector<std::pair<int,std::string>> futureDeps;
    for(const auto& [s,n] : status)
      if(n==1)
	{
	  size_t firstUsage=nPassedDep;
	  while(firstUsage<orderedDep.size() and orderedDep[firstUsage]!=s)
	    firstUsage++;
	  
	  if(firstUsage==orderedDep.size())
	    crash("At nPassedDep %d unable to find the first usage for %s, why is it allocated?",nPassedDep,s.c_str());
	  else
	    futureDeps.emplace_back(firstUsage-nPassedDep,s);
	}
    
    verbosity_lv2_master_printf("futureDeps:\n");
    for(const auto& [delay,name] : futureDeps)
      verbosity_lv2_master_printf(" %zu %s\n",delay,name.c_str());
    
    // \todo simplify
    
    int nLoaded=0;
    for(const auto& s : status)
      nLoaded+=(s.second==1);
    verbosity_lv2_master_printf("n Loaded props: %zu\n",nLoaded);
    
    if((int)futureDeps.size()!=nLoaded)
      crash("unmatched loaded and future deps, %d %zu",nLoaded,futureDeps.size());
    
    master_printf("Needs to accommodate %d more props, %d are loaded and max is %d, needs to offload %d\n",
		  nToAccommodate,nLoaded,nMaxPropsAllocated,nLoaded+nToAccommodate-nMaxPropsAllocated);
    
    std::sort(futureDeps.begin(),futureDeps.end());
    for(size_t i=nMaxPropsAllocated-nToAccommodate;i<futureDeps.size();i++)
      {
	const auto& [delay,name]=futureDeps[i];
	if(delay==0)
	  crash("Unable to offload %s which will be needed immediately!",name.c_str());
	master_printf("Offloading %s which will be used in %d\n",name.c_str(),delay);
	offloadRecallDelete(name,OFFLOAD);
	offloadedList.insert(name);
      }
  }
  
  /// Ensure that a prop is in memory
  void ensureInMemory(const std::string& name)
  {
    if(status[name]==OFFLOADED)
      {
	master_printf("Recalling %s\n",name.c_str());
	offloadRecallDelete(name,RECALL);
      }
  }
  
  //generate a source, wither a wall or a point in the origin
  void generate_original_source(qprop_t* sou,bool skipOnly)
  {
    const rnd_t noise_type=sou->noise_type;
    
    std::unique_ptr<FieldRngOf<spincolor>> drawer;
    if(stoch_source and use_new_generator)
      {
	drawer=std::make_unique<FieldRngOf<spincolor>>(field_rng_stream.getDrawer<spincolor>());
	if(noise_type!=RND_Z4 and noise_type!=RND_Z2)
	  crash("Noise type different from Z4 or Z2 not implemented yet");
	
	if(skipOnly)
	  return;
      }
    
    //consistency check
    if(not stoch_source and (not diluted_spi_source or not diluted_col_source)) crash("for a non-stochastic source, spin and color must be diluted");
    
    //reset all to begin
    for(int i=0;i<nso_spi*nso_col;i++) vector_reset(sou->sp[i]);
    
    spincolor **sou_proxy=nissa_malloc("sou_proxy",nso_spi*nso_col,spincolor*);
    for(int id_so=0;id_so<nso_spi;id_so++)
      for(int ic_so=0;ic_so<nso_col;ic_so++)
	sou_proxy[so_sp_col_ind(id_so,ic_so)]=sou->sp[so_sp_col_ind(id_so,ic_so)];
    
    const int tins=sou->tins;
    
    //NISSA_PARALLEL_LOOP(ivol,0,locVol)
    NISSA_LOC_VOL_LOOP(ivol)
      {
	spincolor c;
	spincolor_put_to_zero(c);
	
	//compute relative coords
	bool is_spat_orig=true;
	coords_t rel_c;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    rel_c[mu]=rel_coord_of_loclx(ivol,mu);
	    if(mu) is_spat_orig&=(rel_c[mu]==0);
	  }
	
	//dilute in space
	int mask=1;
	for(int mu=0;mu<NDIM;mu++) mask&=(rel_c[mu]%diluted_spat_source==0);
	
	//fill colour and spin index 0
	for(int id_si=0;id_si<(diluted_spi_source?1:NDIRAC);id_si++)
	  for(int ic_si=0;ic_si<(diluted_col_source?1:NCOL);ic_si++)
	    {
	      if(stoch_source and mask and (tins==-1 or rel_c[0]==tins))
		{
		  if(use_new_generator)
		    {
		      drawer->fillLocSite(c,ivol);
		      for(int id=0;id<NDIRAC;id++)
			for(int ic=0;ic<NCOL;ic++)
			  {
			    if(noise_type==RND_Z4)
			      z4Transform(c[id][ic]);
			    else
			      z2Transform(c[id][ic]);
			  }
		    }
		  else
		    comp_get_rnd(c[id_si][ic_si],&(loc_rnd_gen[ivol]),noise_type);
		}
	      if(not stoch_source and is_spat_orig and (tins==-1 or rel_c[0]==tins)) complex_put_to_real(c[id_si][ic_si],1);
	    }
	
	//fill other spin indices
	for(int id_so=0;id_so<nso_spi;id_so++)
	  for(int ic_so=0;ic_so<nso_col;ic_so++)
	    for(int id_si=0;id_si<NDIRAC;id_si++)
	      for(int ic_si=0;ic_si<NCOL;ic_si++)
		  if((!diluted_spi_source or (id_so==id_si)) and (!diluted_col_source or (ic_so==ic_si)))
		    complex_copy(sou_proxy[so_sp_col_ind(id_so,ic_so)][ivol][id_si][ic_si],c[diluted_spi_source?0:id_si][diluted_col_source?0:ic_si]);
      }
    //NISSA_PARALLEL_LOOP_END;
    
    //compute the norm2, set borders invalid
    double ori_source_norm2=0;
    for(int id_so=0;id_so<nso_spi;id_so++)
      for(int ic_so=0;ic_so<nso_col;ic_so++)
	{
	  spincolor *s=sou->sp[so_sp_col_ind(id_so,ic_so)];
	  set_borders_invalid(s);
	  ori_source_norm2+=double_vector_glb_norm2(s,locVol);
	}
    if(IS_MASTER_THREAD) sou->ori_source_norm2=ori_source_norm2;
    
    // complex *n=nissa_malloc("n",locVol,complex);
    // spincolor *temp=nissa_malloc("temp",locVol+bord_vol,spincolor);
    // for(int id_so=0;id_so<nso_spi;id_so++)
    //   for(int ic_so=0;ic_so<nso_col;ic_so++)
    // 	{
    // 	  spincolor *s=sou->sp[so_sp_col_ind(id_so,ic_so)];
    // 	  master_printf("eta (0): %lg %lg\n",s[0][0][0][RE],s[0][0][0][IM]);
    // 	  fft4d((complex*)temp,(complex*)s,NDIRAC*NCOL,FFT_PLUS,FFT_NO_NORMALIZE);
	  
    // 	  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    // 	    {
    // 	      n[ivol][RE]=spincolor_norm2(temp[ivol]);
    // 	      n[ivol][IM]=0;
    // 	    }
    // 	  NISSA_PARALLEL_LOOP_END;
	  
    // 	  fft4d(n,n,1,FFT_MINUS,FFT_NORMALIZE);
	  
    // 	  const int64_t nPoints=((tins==-1)?glbVol:glbSpatVol);
	  
    // 	  master_printf("eta+ eta (0): %lg , eta+ eta (1): %lg\n",n[0][RE],n[1][RE]);
    // 	  if(rank==0)
    // 	    n[0][RE]-=(double)NDIRAC*NCOL*nPoints/nso_spi/nso_col;
    // 	  master_printf("eta+ eta (0) after sub: %lg\n",n[0][RE]);
	  
    // 	  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    // 	    {
    // 	      n[ivol][RE]*=n[ivol][RE];
    // 	    }
    // 	  NISSA_PARALLEL_LOOP_END;
	  
    // 	  complex res[1];
    // 	  glb_reduce(res,n,locVol);
	  
    // 	  master_printf("Res: %lg\n",res[0][RE]);
	  
    // 	  double exp=(double)NDIRAC*NCOL*sqr(nPoints)/(nso_spi*nso_col);
    // 	  master_printf("Exp: %lg\n",exp);
	  
    // 	  double dev=res[0][RE]/exp-1;
	  
    // 	  master_printf("Dev: %lg\n",dev);
    // 	}
    
    // nissa_free(temp);
    // nissa_free(n);
    
    nissa_free(sou_proxy);
  }
  
  //Generate all the original sources
  void generate_original_sources(int ihit,bool skipOnly)
  {
    for(size_t i=0;i<ori_source_name_list.size();i++)
      {
	std::string &name=ori_source_name_list[i];
	master_printf("Generating source \"%s\"\n",name.c_str());
	qprop_t *q=&Q[name];
	q->alloc_storage();
	generate_original_source(q,skipOnly);
	
	if(not skipOnly)
	  for(int id_so=0;id_so<nso_spi;id_so++)
	    for(int ic_so=0;ic_so<nso_col;ic_so++)
	      {
		//combine the filename
		const std::string path=combine("%s/hit%d_source%s_idso%d_icso%d",outfolder,ihit,name.c_str(),id_so,ic_so);
		
		int isou=so_sp_col_ind(id_so,ic_so);
		spincolor *sou=(*q)[isou];
		
		ReadWriteRealVector<spincolor> rw(sou,path);
		
		//if the prop exists read it
		if(rw.canLoad())
		  {
		    master_printf("  loading the source dirac index %d, color %d\n",id_so,ic_so);
		    START_TIMING(read_prop_time,nread_prop);
		    rw.read();
		    STOP_TIMING(read_prop_time);
		  }
		else
		  {
		    master_printf("  file %s not available, skipping loading\n",path.c_str());
		    
		    //and store if needed
		    if(q->store)
		      {
			master_printf("  writing the source dirac index %d, color %d\n",id_so,ic_so);
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
    master_printf("\n=== Hit %d/%d ====\n",ihit+1,nhits);
    
    if(use_new_generator)
      {
	for(int mu=0;mu<NDIM;mu++)
	  {
	    using C=double[1];
	    C c;
	    field_rng_stream.drawScalar(c);
	    source_coord[mu]=c[0]*glbSize[mu];
	  }
      }
    else
      source_coord=generate_random_coord();
    
    if(stoch_source) master_printf(" source time: %d\n",source_coord[0]);
    else             master_printf(" point source coords: %d %d %d %d\n",source_coord[0],source_coord[1],source_coord[2],source_coord[3]);
    if(need_photon)
      {
	if(skip)
	  generate_photon_source(photon_eta);
	else
	  generate_photon_stochastic_propagator(ihit);
      }
    generate_original_sources(ihit,skip);
  }
  
  /// Perform one of step: dryRun, or eval
  void internalRun(const int runStep,const size_t iHit)
  {
    computedProps.clear();
    computed2pts.clear();
    offloadedList.clear();
    status.clear();
    
    /// Keep backup of the propagators dependencies
    const std::map<std::string,std::set<std::string>> propDepBack=propDep;
    
    // Mark all the original sources as in memory and computed
    if(runStep==2)
      for(const auto& s : ori_source_name_list)
	{
	  status[s]=IN_MEMORY;
	  computedProps.insert(s);
	}
    
    for(size_t i=0;i<qprop_name_list.size();i++)
      if(propDep[qprop_name_list[i]].size()==0)
	master_printf("Skipping generation of prop %s as it has no dependencies\n",qprop_name_list[i].c_str());
      else
	{
	  //get names
	  std::string name=qprop_name_list[i];
	  if(runStep==2)
	    master_printf("Computing prop %s\n",name.c_str());
	  
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
	      status[name]=IN_MEMORY;
	    }
	  
	  computedProps.insert(name);
	  
	  // Remove the dependencies from all sources
	  for(const auto& [s,w] : q.source_terms)
	    if(propDep[s].erase(name)!=1)
	      crash("unable to remove the dependency %s of %s!",name.c_str(),s.c_str());
	    else
	      verbosity_lv2_master_printf("Erased dependency of %s from %s which has still %zu dependencies\n",name.c_str(),s.c_str(),propDep[s].size());
	  
	  if(runStep==2)
	    {
	      nPassedDep+=Q[name].source_terms.size();
	      eraseUnnededProps();
	    }
	  
	  for(size_t icombo=0;icombo<mes2pts_contr_map.size();icombo++)
	    if(const std::string& a=mes2pts_contr_map[icombo].a,
	       b=mes2pts_contr_map[icombo].b;
	       computedProps.count(a) and
	       computedProps.count(b) and
	       computed2pts.count(icombo)==0)
	      {
		if(runStep==2)
		  master_printf("Can compute contraction: %s %s -> %s, %zu/%zu\n",
			      a.c_str(),b.c_str(),mes2pts_contr_map[icombo].name.c_str(),computed2pts.size(),mes2pts_contr_map.size());
		
		// Insert the correlation in the computed list
		computed2pts.insert(icombo);
		
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
	for(const auto& s : orderedDep)
	  verbosity_lv2_master_printf("  %s\n",s.c_str());
	verbosity_lv2_master_printf("Dependencies end\n");
      }
  }
  
  void run(const size_t iHit)
  {
    for(int runStep=1;runStep<=2;runStep++)
      internalRun(runStep,iHit);
    
    std::ostringstream os;
    for(size_t iContr=0;iContr<mes2pts_contr_map.size();iContr++)
      if(const auto& m=mes2pts_contr_map[iContr];computed2pts.find(iContr)==computed2pts.end())
	os<<" not computed corr "<<m.name<<" between "<<m.a<<" and "<<m.b<<std::endl;
    
    if(os.str().size())
      crash("%s",os.str().c_str());
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
