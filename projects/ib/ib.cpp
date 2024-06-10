#include <nissa.hpp>

#include "conf.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

#include <sstream>
#include <complex>

using namespace nissa;

///////////////////////////////// initialise the library, read input file, allocate /////////////////////////////////////

void init_simulation(int narg,char **arg)
{
  //check argument
  if(narg<2) crash("Use: %s input_file [stop_path]|periodic/antiperiodic|store/load_photons",arg[0]);
  
  const char *path=arg[1];
  
  //parse the rest of the args
  for(int iarg=2;iarg<narg;iarg++)
    {
      bool parsed=false;
      master_printf("parsing argument %d: '%s'\n",iarg,arg[iarg]);
      
      //check if we passed "periodic"
      if(not parsed and not strcasecmp(arg[iarg],"periodic"))
	{
	  temporal_bc=PERIODIC_BC;
	  master_printf(" Setting temporal bc to %lg = 'periodic'\n",temporal_bc);
	  parsed=true;
	}
      
      //check if we passed "antiperiodic"
      if(not parsed and not strcasecmp(arg[iarg],"antiperiodic"))
	{
	  temporal_bc=ANTIPERIODIC_BC;
	  master_printf(" Setting temporal bc to %lg = 'antiperiodic'\n",temporal_bc);
	  parsed=true;
	}
      
      //check if we passed "store_photons"
      if(not parsed and not strcasecmp(arg[iarg],"store_photons"))
	{
	  master_printf(" Store the photon fields\n");
	  store_photons=true;
	  parsed=true;
	}
      
      //check if we passed "load_photons"
      if(not parsed and not strcasecmp(arg[iarg],"load_photons"))
	{
	  master_printf(" Load (if present) the photon fields\n");
	  load_photons=true;
	  parsed=true;
	}
      
      //otherwise take this as a suffix for stop, running and finished filenames
      if(not parsed)
	{
	  std::string app(arg[iarg]);
	  stop_path+="_"+app;
	  running_filename+="_"+app;
	  finished_filename+="_"+app;
	  master_printf("Adding to stop,finished,and running filenames the suffix: '%s'\n",arg[iarg]);
	  master_printf("Stop filename: '%s'\n",stop_path.c_str());
	  master_printf("Running filename: '%s'\n",running_filename.c_str());
	  master_printf("Finished filename: '%s'\n",finished_filename.c_str());
	  parsed=true;
	}
    }
  
  //open input file
  open_input(path);
  
  //init the grid
  read_init_grid();
  
  //Wall time
  read_str_double("WallTime",&wall_time);
  
  //Sources
  read_seed_start_random();
  read_stoch_source();
  read_set_dilutions_if_stoch_source(stoch_source);
  read_nhits();
  int nsources;
  read_str_int("NSources",&nsources);
  ori_source_name_list.resize(nsources*ncopies);
  //discard header
  expect_str("Name");
  if(stoch_source) expect_str("NoiseType");
  expect_str("Tins");
  expect_str("Store");
  //loop over sources
  for(int isource=0;isource<nsources;isource++)
    {
      //name
      char name[1024];
      read_str(name,1024);
      //noise type and tins
      rnd_t noise_type=RND_ALL_PLUS_ONE;
      int tins=ALL_TIMES;
      if(stoch_source)
	{
	  char str_noise_type[20];
	  read_str(str_noise_type,20);
	  noise_type=convert_str_to_rnd_t(str_noise_type);
	}
      read_int(&tins);
      //store
      int store_source;
      read_int(&store_source);
      
      //add
      for(int icopy=0;icopy<ncopies;icopy++)
	{
	  char suffix[128]="";
	  if(ncopies>1) sprintf(suffix,"_copy%d",icopy);
	  
	  char fullName[1024+129];
	  sprintf(fullName,"%s%s",name,suffix);
	  ori_source_name_list[isource+nsources*icopy]=fullName;
	  Q[fullName].init_as_source(noise_type,tins,0,store_source);
	}
    }
  
  read_twisted_run();
  read_clover_run();
  
  //NProps
  int nprops;
  read_str_int("NProps",&nprops);
  qprop_name_list.resize(nprops*ncopies);
  //Discard header
  expect_str("Name");
  expect_str("Ins");
  expect_str("SourceName");
  expect_str("Tins");
  expect_str("Kappa");
  if(twisted_run)
    {
      expect_str("Mass");
      expect_str("R");
    }
  expect_str("Charge");
  read_iso_theta();
  expect_str("Residue");
  expect_str("Store");
  for(int iq=0;iq<nprops;iq++)
    {
      //name
      char name[1024];
      read_str(name,1024);
      master_printf("Read variable 'Name' with value: %s\n",name);
      
      //ins name
      char ins[INS_TAG_MAX_LENGTH+1];
      read_str(ins,INS_TAG_MAX_LENGTH);
      master_printf("Read variable 'Ins' with value: %s\n",ins);
      
      //source_name
      std::vector<source_term_t> source_terms;
      char source_name[1024];
      read_str(source_name,1024);
      master_printf("Read variable 'SourceName' with value: %s\n",source_name);
      
      bool multi_source=(strcasecmp(source_name,"LINCOMB")==0);
      
      int nsources;
      if(multi_source)
	{
	  read_int(&nsources);
	  master_printf("Read variable 'NSources' with value: %d\n",nsources);
	}
      else
	nsources=1;
      
      for(int isource=0;isource<nsources;isource++)
	{
	  std::pair<double,double> weight={1.0,0.0};
	  if(multi_source)
	    {
	      read_str(source_name,1024);
	      
	      master_printf("Read variable 'Sourcename' with value: %s\n",source_name);
	      
	      //read weight as string
	      char sweight[128];
	      read_str(sweight,128);
	      
	      //convert to double
	      std::istringstream is(sweight);
	      std::complex<double> cweight;
	      is>>cweight;
	      
	      weight.first=cweight.real();
	      weight.second=cweight.imag();
	    }
	  
	  source_terms.push_back(std::make_pair(source_name,weight));
	}
      
      //insertion time
      int tins;
      read_int(&tins);
      master_printf("Read variable 'Tins' with value: %d\n",tins);
      
      double kappa=0.125,mass=0.0,charge=0,residue=1e-16;
      momentum_t theta;
      char ext_field_path[32]="";
      theta[0]=temporal_bc;
      for(int mu=1;mu<NDIM;mu++) theta[mu]=0;
      int r=0,store_prop=0;
      
      bool decripted=false;
      
      if(strcasecmp(ins,ins_tag[PROP])==0 or strcasecmp(ins,ins_tag[DIROP])==0)
	{
	  decripted=true;
	  
	  read_double(&kappa);
	  master_printf("Read variable 'Kappa' with value: %lg\n",kappa);
	  if(twisted_run)
	    {
	      read_double(&mass);
	      master_printf("Read variable 'Mass' with value: %lg\n",mass);
	      read_int(&r);
	      master_printf("Read variable 'R' with value: %d\n",r);
	      
	      //include tau in the mass
	      mass*=tau3[r];
	    }
	  read_double(&charge);
	  master_printf("Read variable 'Charge' with value: %lg\n",charge);
	  read_theta(theta);
	  if(strcasecmp(ins,ins_tag[DIROP])!=0)
	    {
	      read_double(&residue);
	      master_printf("Read variable 'Residue' with value: %lg\n",residue);
	    }
	}
      
      //read phasing
      if(strcasecmp(ins,ins_tag[PHASING])==0)
	{
	  decripted=true;
	  
	  read_theta(theta);
	}
      
      bool vph=false;
      for(const auto& possIns : {VPHOTON0,VPHOTON1,VPHOTON2,VPHOTON3,
				 VBHOTON0,VBHOTON1,VBHOTON2,VBHOTON3})
	vph|=(strcasecmp(ins,ins_tag[possIns])==0);
      
      bool ph=false;
      for(const auto& possIns : {CVEC0,CVEC1,CVEC2,CVEC3})
	ph|=(strcasecmp(ins,ins_tag[possIns])==0);
      
      if(vph)
	{
	  decripted=true;
	  
	  read_double(&mass);
	  master_printf("Read variable 'Mass' with value: %lg\n",mass);
	  
	  
	  read_int(&r);
	  master_printf("Read variable 'R' with value: %d\n",r);
	  
	  read_theta(theta);
	}
      
      //read smearing
      if(strcasecmp(ins,ins_tag[SMEARING])==0 or strcasecmp(ins,ins_tag[WFLOW])==0 or strcasecmp(ins,ins_tag[BACK_WFLOW])==0)
	{
	  decripted=true;
	  
	  read_double(&kappa);
	  master_printf("Read variable 'Kappa' with value: %lg\n",kappa);
	  
	  read_int(&r);
	  master_printf("Read variable 'R' with value: %d\n",r);
	  
	  read_theta(theta);
	}
      
      //read smearing
      double kappa1=0.0,kappa2=0.0,kappa3=0.0;
      if(strcasecmp(ins,ins_tag[ANYSM])==0)
	{
	  decripted=true;
	  
	  read_double(&kappa1);
	  master_printf("Read variable 'Kappa1' with value: %lg\n",kappa1);
	  
	  read_double(&kappa2);
	  master_printf("Read variable 'Kappa2' with value: %lg\n",kappa2);
	  
	  read_double(&kappa3);
	  master_printf("Read variable 'Kappa3' with value: %lg\n",kappa3);
	  
	  read_int(&r);
	  master_printf("Read variable 'R' with value: %d\n",r);
	  
	  read_theta(theta);
	}
      double kappa_asymm[4]={0.0,kappa1,kappa2,kappa3};
      
      //everything else
      if(not decripted)
	{
	  //external source
	  if(strcasecmp(ins,ins_tag[EXT_FIELD])==0)
	    {
	      read_str(ext_field_path,32);
	      master_printf("Read variable 'ext_field_path' with value: %s\n",ext_field_path);
	    }
	  
	  if(twisted_run)
	    {
	      read_int(&r);
	      master_printf("Read variable 'R' with value: %d\n",r);
	    }
	  read_double(&charge);
	  master_printf("Read variable 'Charge' with value: %lg\n",charge);
	  
	  if(ph)
	    read_theta(theta);
	}
      
      read_int(&store_prop);
      master_printf("Read variable 'Store' with value: %d\n",store_prop);
      
      for(int icopy=0;icopy<ncopies;icopy++)
	{
	  char suffix[128]="";
	  if(ncopies>1) sprintf(suffix,"_copy%d",icopy);
	  
	  std::vector<source_term_t> source_full_terms=source_terms;
	  for(auto& [name,weight] : source_full_terms)
	    {
	      name+=suffix;
	      if(Q.find(name)==Q.end()) crash("unable to find source %s",name.c_str());
	    }
	  
	  char fullName[1024+129];
	  sprintf(fullName,"%s%s",name,suffix);
	  if(Q.find(fullName)!=Q.end()) crash("name \'%s\' already included",fullName);
	  
	  Q[fullName].init_as_propagator(ins_from_tag(ins),source_full_terms,tins,residue,kappa,kappa_asymm,mass,ext_field_path,r,charge,theta,store_prop);
	  qprop_name_list[icopy+ncopies*iq]=fullName;
	}
    }
  
  read_photon_pars();
  
  read_free_theory_flag();
  read_random_gauge_transform();
  read_Landau_gauge_fix();
  read_store_conf();
  
  read_loc_hadr_curr();
  read_loc_muon_curr();
  
  //mesons
  read_mes2pts_contr_pars();
  
  //meslept
  read_meslep_contr_pars();
  
  //barions
  set_Cg5();
  read_bar2pts_contr_pars();
  
  //meson handcuffs
  read_handcuffs_contr_pars();
  
  //read the fast Fourier transform parameters
  read_fft_prop_pars();
  
  //smearing
  read_ape_smearing_pars();
  
  read_ngauge_conf();
  
  ///////////////////// finished reading apart from conf list ///////////////
  
  if(clover_run)
    {
      Cl=nissa_malloc("Cl",locVol,clover_term_t);
      invCl=nissa_malloc("invCl",locVol,inv_clover_term_t);
    }
  
  allocate_loop_source();
  allocate_photon_fields();
  
  loc_contr=nissa_malloc("loc_contr",locVol,complex);
  
  allocate_mes2pts_contr();
  allocate_handcuffs_contr();
  
  nmeslep_corr=nquark_lep_combos*nindep_meslep_weak;
  meslep_hadr_part=nissa_malloc("hadr",locVol,spinspin);
  meslep_contr=nissa_malloc("meslep_contr",glbSize[0]*nindep_meslep_weak*nmeslep_proj*nmeslep_corr,complex);
  
  allocate_bar2pts_contr();
  
  allocate_L_prop();
  
  lock_file.init();
}

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
	computedProps.insert(name);
	
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
  void run(const int runStep,const size_t iHit)
  {
    status.clear();
    
    /// Keep backup of the propagators dependencies
    const std::map<std::string,std::set<std::string>> propDepBack=propDep;
    
    // Mark all the original sources as in memory
    if(runStep==2)
      for(const auto& s : ori_source_name_list)
	status[s]=IN_MEMORY;
    
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
    
    computedProps.clear();
    computed2pts.clear();
    propDep=propDepBack;
    offloadedList.clear();
    
    if(runStep==1)
      {
	for(const auto& s : orderedDep)
	  verbosity_lv2_master_printf("  %s\n",s.c_str());
	verbosity_lv2_master_printf("Dependencies end\n");
      }
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

//carry out a single hit
void hit_loop(int iHit)
{
  HitLooper hitLooper;
  hitLooper.start_hit(iHit);
  
  for(int runStep=1;runStep<=2;runStep++)
    hitLooper.run(runStep,iHit);
  
  master_printf("NOffloaded: %d\n",hitLooper.nOffloaded);
  master_printf("NRecalled: %d\n",hitLooper.nRecalled);
}

//close deallocating everything
void close()
{
  print_statistics();
  
  Q.clear();
  
  free_photon_fields();
  free_loop_source();
  free_L_prop();
  
  free_mes2pts_contr();
  free_handcuffs_contr();
  
  nissa_free(loc_contr);
  
  nissa_free(meslep_hadr_part);
  nissa_free(meslep_contr);
  nissa_free(lep_contr_iq1);
  nissa_free(lep_contr_iq2);
  nissa_free(leps);
  nissa_free(lep_energy);
  nissa_free(neu_energy);
  
  for(auto &f : fft_filterer)
    if(f.nfft_filtered) f.fft_filter_remap.destroy();
  
  if(clover_run)
    {
      nissa_free(Cl);
      nissa_free(invCl);
    }
  free_bar2pts_contr();
  
  free_confs();
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //init simulation according to input file
  init_simulation(narg,arg);
  
  constexpr char NMAX_PROPS_ALLOCATED_STR[]="NMAX_PROPS_ALLOCATED";
  if(const char* nMaxAllocatedStr=getenv(NMAX_PROPS_ALLOCATED_STR))
    {
      nMaxPropsAllocated=atoi(nMaxAllocatedStr);
      master_printf("NMAX_PROPS_ALLOCATED=%d\n",nMaxPropsAllocated);
    }
  else
    {
      master_printf("No maximum number of propagators to be allocated passed\n");
      master_printf("Optionally specify the maximal number of propagators to be allocated by exporting %s\n",NMAX_PROPS_ALLOCATED_STR);
    }
  
  //loop over the configs
  int iconf=0;
  while(read_conf_parameters(iconf,finish_file_present))
    {
      for(int ihit=0;ihit<nhits;ihit++)
	{
	  hit_loop(ihit);
	  compute_contractions(); //not working, here only to emit errors
	  propagators_fft(ihit); // same
	}
      
      free_confs();
      print_contractions();
      
      mark_finished();
    }
  
  //close the simulation
  tot_prog_time+=take_time();
  close();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
