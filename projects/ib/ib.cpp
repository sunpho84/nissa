#include <nissa.hpp>

#include "conf.hpp"
#include "hit.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

#include <sstream>
#include <complex>

using namespace nissa;

///////////////////////////////// initialise the library, read input file, allocate /////////////////////////////////////

extern int runningUpdateTime;

void init_simulation(int narg,char **arg)
{
  //check argument
  if(narg<2)
    CRASH("Use: %s input_file [stop_path]|periodic/antiperiodic|store/load_photons",arg[0]);
  
  const char RUNNING_UPDATE_TIME_STRING[]=
    "RUNNING_UPDATE_TIME";
  if(const char* p=getenv(RUNNING_UPDATE_TIME_STRING))
    {
      runningUpdateTime=atoi(p);
      MASTER_PRINTF("Time between running file update set to: %d s\n",runningUpdateTime);
    }
  else
    {
      MASTER_PRINTF("Time between running file update set to default value: %d s\n",runningUpdateTime);
      MASTER_PRINTF("To change it, export the envirnoment variable %s\n",RUNNING_UPDATE_TIME_STRING);
    }
  
  const char *path=arg[1];
  
  //parse the rest of the args
  for(int iarg=2;iarg<narg;iarg++)
    {
      bool parsed=false;
      MASTER_PRINTF("parsing argument %d: '%s'\n",iarg,arg[iarg]);
      
      //check if we passed "periodic"
      if(not parsed and not strcasecmp(arg[iarg],"periodic"))
	{
	  temporal_bc=PERIODIC_BC;
	  MASTER_PRINTF(" Setting temporal bc to %lg = 'periodic'\n",temporal_bc);
	  parsed=true;
	}
      
      //check if we passed "antiperiodic"
      if(not parsed and not strcasecmp(arg[iarg],"antiperiodic"))
	{
	  temporal_bc=ANTIPERIODIC_BC;
	  MASTER_PRINTF(" Setting temporal bc to %lg = 'antiperiodic'\n",temporal_bc);
	  parsed=true;
	}
      
      //check if we passed "store_photons"
      if(not parsed and not strcasecmp(arg[iarg],"store_photons"))
	{
	  MASTER_PRINTF(" Store the photon fields\n");
	  store_photons=true;
	  parsed=true;
	}
      
      //check if we passed "load_photons"
      if(not parsed and not strcasecmp(arg[iarg],"load_photons"))
	{
	  MASTER_PRINTF(" Load (if present) the photon fields\n");
	  load_photons=true;
	  parsed=true;
	}
      
      //otherwise take this as a suffix for stop, running and finished filenames
      if(not parsed)
	{
	  std::string app(arg[iarg]);
	  stopPath+="_"+app;
	  runningFilename+="_"+app;
	  finishedFilename+="_"+app;
	  MASTER_PRINTF("Adding to stop,finished,and running filenames the suffix: '%s'\n",arg[iarg]);
	  MASTER_PRINTF("Stop path: '%s'\n",stopPath.c_str());
	  MASTER_PRINTF("Running filename: '%s'\n",runningFilename.c_str());
	  MASTER_PRINTF("Finished filename: '%s'\n",finishedFilename.c_str());
	  MASTER_PRINTF("Partial data filename: '%s'\n",partialDataFilename.c_str());
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
  ori_source_name_list.resize(nsources*nCopies);
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
      for(int icopy=0;icopy<nCopies;icopy++)
	{
	  char suffix[128]="";
	  if(nCopies>1) sprintf(suffix,"_copy%d",icopy);
	  
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
  qprop_name_list.resize(nprops*nCopies);
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
      MASTER_PRINTF("Read variable 'Name' with value: %s\n",name);
      
      //ins name
      char ins[INS_TAG_MAX_LENGTH+1];
      read_str(ins,INS_TAG_MAX_LENGTH);
      MASTER_PRINTF("Read variable 'Ins' with value: %s\n",ins);
      
      //source_name
      std::vector<source_term_t> source_terms;
      char source_name[1024];
      read_str(source_name,1024);
      MASTER_PRINTF("Read variable 'SourceName' with value: %s\n",source_name);
      
      bool multi_source=(strcasecmp(source_name,"LINCOMB")==0);
      
      int nsources;
      if(multi_source)
	{
	  read_int(&nsources);
	  MASTER_PRINTF("Read variable 'NSources' with value: %d\n",nsources);
	}
      else
	nsources=1;
      
      for(int isource=0;isource<nsources;isource++)
	{
	  std::pair<double,double> weight={1.0,0.0};
	  if(multi_source)
	    {
	      read_str(source_name,1024);
	      
	      MASTER_PRINTF("Read variable 'Sourcename' with value: %s\n",source_name);
	      
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
      int tins=-1;
      if(strcasecmp(ins,ins_tag[DEL_POS])!=0)
	{
	  read_int(&tins);
	  MASTER_PRINTF("Read variable 'Tins' with value: %d\n",tins);
	}
      
      double kappa=0.125,mass=0.0,charge=0,residue=1e-16;
      Momentum theta;
      char ext_field_path[32]="";
      theta[0]=temporal_bc;
      for(int mu=1;mu<NDIM;mu++) theta[mu]=0;
      int r=0,store_prop=0;
      
      bool decripted=false;
      
      if(strcasecmp(ins,ins_tag[PROP])==0 or
	 strcasecmp(ins,ins_tag[DIROP])==0 or
	 strcasecmp(ins,ins_tag[LEP_LOOP])==0)
	{
	  decripted=true;
	  
	  read_double(&kappa);
	  MASTER_PRINTF("Read variable 'Kappa' with value: %lg\n",kappa);
	  if(twisted_run)
	    {
	      read_double(&mass);
	      MASTER_PRINTF("Read variable 'Mass' with value: %lg\n",mass);
	      
	      read_int(&r);
	      MASTER_PRINTF("Read variable 'R' with value: %d\n",r);
	      
	      //include tau in the mass
	      mass*=tau3[r];
	    }
	  read_double(&charge);
	  MASTER_PRINTF("Read variable 'Charge' with value: %lg\n",charge);
	  read_theta(theta);
	  if(strcasecmp(ins,ins_tag[DIROP])!=0)
	    {
	      read_double(&residue);
	      MASTER_PRINTF("Read variable 'Residue' with value: %lg\n",residue);
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
      for(const auto& possIns : {CVEC0,CVECBW0,CVECFW0,
				 CVEC1,
				 CVEC2,
				 CVEC3})
	ph|=(strcasecmp(ins,ins_tag[possIns])==0);
      
      if(vph)
	{
	  decripted=true;
	  
	  read_double(&mass);
	  MASTER_PRINTF("Read variable 'Mass' with value: %lg\n",mass);
	  
	  
	  read_int(&r);
	  MASTER_PRINTF("Read variable 'R' with value: %d\n",r);
	  
	  read_theta(theta);
	}
      
      //read smearing
      if(strcasecmp(ins,ins_tag[SMEARING])==0 or strcasecmp(ins,ins_tag[WFLOW])==0 or strcasecmp(ins,ins_tag[BACK_WFLOW])==0)
	{
	  decripted=true;
	  
	  read_double(&kappa);
	  MASTER_PRINTF("Read variable 'Kappa' with value: %lg\n",kappa);
	  
	  read_int(&r);
	  MASTER_PRINTF("Read variable 'R' with value: %d\n",r);
	  
	  read_theta(theta);
	}
      
      //read smearing
      double kappa1=0.0,kappa2=0.0,kappa3=0.0;
      if(strcasecmp(ins,ins_tag[ANYSM])==0)
	{
	  decripted=true;
	  
	  read_double(&kappa1);
	  MASTER_PRINTF("Read variable 'Kappa1' with value: %lg\n",kappa1);
	  
	  read_double(&kappa2);
	  MASTER_PRINTF("Read variable 'Kappa2' with value: %lg\n",kappa2);
	  
	  read_double(&kappa3);
	  MASTER_PRINTF("Read variable 'Kappa3' with value: %lg\n",kappa3);
	  
	  read_int(&r);
	  MASTER_PRINTF("Read variable 'R' with value: %d\n",r);
	  
	  read_theta(theta);
	}
      double kappa_asymm[4]={0.0,kappa1,kappa2,kappa3};
      
      if(strcasecmp(ins,ins_tag[DEL_POS])==0)
	{
	  Coords c;
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      read_int(&c[mu]);
	      if(c[mu]<0)
		CRASH("dir %d has been chosen negative value %d",mu,c[mu]);
	      if(c[mu]>=glbSize[mu])
		CRASH("dir %d has been chosen larger than glb size %d",mu,glbSize[mu]);
	    }
	  
	  r=glblxOfCoord(c);
	  MASTER_PRINTF("Choosing coords {%d,%d,%d,%d} corresponding to site %d\n",c[0],c[1],c[2],c[3],r);
	  
	  decripted=true;
	}
      
      if(strcasecmp(ins,ins_tag[DEL_SPIN])==0)
	{
	  read_int(&r);
	  MASTER_PRINTF("Choosing spin %d\n",r);
	  decripted=true;
	}
      
      if(strcasecmp(ins,ins_tag[DEL_COL])==0)
	{
	  read_int(&r);
	  MASTER_PRINTF("Choosing color %d\n",r);
	  decripted=true;
	}
      
      //everything else
      if(not decripted)
	{
	  //external source
	  if(strcasecmp(ins,ins_tag[EXT_FIELD])==0)
	    {
	      read_str(ext_field_path,32);
	      MASTER_PRINTF("Read variable 'ext_field_path' with value: %s\n",ext_field_path);
	    }
	  
	  if(twisted_run)
	    {
	      read_int(&r);
	      MASTER_PRINTF("Read variable 'R' with value: %d\n",r);
	    }
	  read_double(&charge);
	  MASTER_PRINTF("Read variable 'Charge' with value: %lg\n",charge);
	  
	  if(ph)
	    read_theta(theta);
	}
      
      read_int(&store_prop);
      MASTER_PRINTF("Read variable 'Store' with value: %d\n",store_prop);
      
      for(int icopy=0;icopy<nCopies;icopy++)
	{
	  char suffix[128]="";
	  if(nCopies>1) sprintf(suffix,"_copy%d",icopy);
	  
	  std::vector<source_term_t> source_full_terms=source_terms;
	  for(auto& [name,weight] : source_full_terms)
	    {
	      name+=suffix;
	      if(Q.find(name)==Q.end()) CRASH("unable to find source %s",name.c_str());
	    }
	  
	  char fullName[1024+129];
	  sprintf(fullName,"%s%s",name,suffix);
	  if(Q.find(fullName)!=Q.end()) CRASH("name \'%s\' already included",fullName);
	  
	  Q[fullName].init_as_propagator(ins_from_tag(ins),source_full_terms,tins,residue,kappa,kappa_asymm,mass,ext_field_path,r,charge,theta,store_prop);
	  qprop_name_list[icopy+nCopies*iq]=fullName;
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
  
  constexpr bool decryptNameTest=false;
  if(decryptNameTest)
    for(const mes_contr_t& m : mes2ptsContr)
      {
	int bwOrder=0;
	std::string line;
	for(const std::string& u : {m.a,m.b})
	  {
	    std::string cur=u;
	    while(not Q[cur].is_source)
	      {
		std::string wh;
		switch(Q[cur].insertion)
		  {
		  case PROP:
		    wh="-";
		    break;
		  case SCALAR:
		    wh="G[0]";
		    break;
		  case PSEUDO:
		    wh="G[5]";
		    break;
		  case GAMMA:
		    wh="G["+std::to_string(Q[cur].r)+"]";
		    break;
		  default:
		    break;
		  }
		
		if(const int t=Q[cur].tins;t!=-1)
		  wh+="|t="+std::to_string(t)+"|";
		
		if(bwOrder==0)
		  line=wh+line;
		else
		  line=line+wh;
		
		cur=Q[cur].source_terms.front().first;
	      }
	    
	    line+="     ";
	    bwOrder++;
	  }
	MASTER_PRINTF("Decrypting %s = %s\n",m.name.c_str(),line.c_str());
      }
  
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
      Cl=new LxField<clover_term_t>("Cl");
      invCl=new LxField<inv_clover_term_t>("invCl");
    }
  
  allocate_loop_source();
  allocate_photon_fields();
  
  loc_contr=new LxField<complex>("loc_contr");
  
  for(mes_contr_t& m : mes2ptsContr)
    m.alloc();
  
  nmeslep_corr=nquark_lep_combos*nindep_meslep_weak;
  meslep_hadr_part=nissa_malloc("hadr",locVol,spinspin);
  meslep_contr=nissa_malloc("meslep_contr",glbSize[0]*nindep_meslep_weak*nmeslep_proj*nmeslep_corr,complex);
  
  allocate_bar2pts_contr();
  
  allocate_L_prop();
  
  lock_file.init();
}

//close deallocating everything
void close()
{
  print_statistics();
  
  Q.clear();
  
  cachedLepLoops.clear();
  
  free_photon_fields();
  free_loop_source();
  free_L_prop();
  
  for(mes_contr_t& m : mes2ptsContr)
    m.unalloc();
  freeHandcuffsContr();
  
  delete loc_contr;
  
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
      delete Cl;
      delete invCl;
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
  
  constexpr char PRESERVE_PARTIAL_DATA_STR[]="PRESERVE_PARTIAL_DATA";
  if(const char* preservePartialDataStr=getenv(PRESERVE_PARTIAL_DATA_STR))
    {
      preservePartialData=atoi(preservePartialDataStr);
      MASTER_PRINTF("Preserve partial data: %d\n",preservePartialData);
    }
  else
    MASTER_PRINTF("Optionally preserve the partial data by exporting %s\n",PRESERVE_PARTIAL_DATA_STR);
  
  constexpr char NMAX_PROPS_ALLOCATED_STR[]="NMAX_PROPS_ALLOCATED";
  if(const char* nMaxAllocatedStr=getenv(NMAX_PROPS_ALLOCATED_STR))
    {
      nMaxPropsAllocated=atoi(nMaxAllocatedStr);
      MASTER_PRINTF("NMAX_PROPS_ALLOCATED=%d\n",nMaxPropsAllocated);
    }
  else
    {
      MASTER_PRINTF("No maximum number of propagators to be allocated passed\n");
      MASTER_PRINTF("Optionally specify the maximal number of propagators to be allocated by exporting %s\n",NMAX_PROPS_ALLOCATED_STR);
    }
  
  constexpr char DO_NOT_AVERAGE_HITS_STR[]="DO_NOT_AVERAGE_HITS";
  doNotAverageHits=getenv(DO_NOT_AVERAGE_HITS_STR)!=nullptr;
  if(doNotAverageHits)
    MASTER_PRINTF("%s exported, not averaging hits\n",DO_NOT_AVERAGE_HITS_STR);
  else
    MASTER_PRINTF("Averaging hits, export %s if needed otherwise\n",DO_NOT_AVERAGE_HITS_STR);
  
  //loop over the configs
  int iconf=0;
  while(read_conf_parameters(iconf) and not fileExists(finishedPath()))
    {
      HitLooper hitLooper;
      const int nHitsDoneSoFar=
	hitLooper.maybeLoadPartialData();
      
      if(nHitsDoneSoFar)
	MASTER_PRINTF("Found partial file with %d hits\n",nHitsDoneSoFar);
      
      for(int iHit=0;iHit<nHits;iHit++)
	{
	  const bool skip=
	    iHit<nHitsDoneSoFar;
	  
	  hitLooper.start_hit(iHit,skip);
	  if(skip)
	    MASTER_PRINTF("Skipping\n");
	  else
	    {
	      hitLooper.run(iHit);
	      
	      compute_contractions(); //not working, here only to emit errors
	      propagators_fft(iHit); // same
	      
	      if(doNotAverageHits)
		print_contractions(iHit);
	      
	      hitLooper.writePartialData(iHit+1);
	    }
	}
      
      MASTER_PRINTF("NOffloaded: %d\n",hitLooper.nOffloaded);
      MASTER_PRINTF("NRecalled: %d\n",hitLooper.nRecalled);
      
      free_confs();
      
      if(not doNotAverageHits)
	print_contractions();
      
      finalizeConf(hitLooper);
    }
  
  //close the simulation
  tot_prog_time+=take_time();
  close();
}

int main(int narg,char **arg)
{
  initNissa(narg,arg);
  in_main(narg,arg);
  closeNissa();
  
  return 0;
}
