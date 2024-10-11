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
      int tins=-1;
      if(strcasecmp(ins,ins_tag[DEL_POS])!=0)
	{
	  read_int(&tins);
	  master_printf("Read variable 'Tins' with value: %d\n",tins);
	}
      
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
      
      if(strcasecmp(ins,ins_tag[DEL_POS])==0)
	{
	  coords_t c;
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      read_int(&c[mu]);
	      if(c[mu]<0)
		crash("dir %d has been chosen negative value %d",mu,c[mu]);
	      if(c[mu]>=glbSize[mu])
		crash("dir %d has been chosen larger than glb size %d",mu,glbSize[mu]);
	    }
	  
	  r=glblx_of_coord(c);
	  master_printf("Choosing coords {%d,%d,%d,%d} corresponding to site %d\n",c[0],c[1],c[2],c[3],r);
	  
	  decripted=true;
	}
      
      if(strcasecmp(ins,ins_tag[DEL_SPIN])==0)
	{
	  read_int(&r);
	  master_printf("Choosing spin %d\n",r);
	  decripted=true;
	}
      
      if(strcasecmp(ins,ins_tag[DEL_COL])==0)
	{
	  read_int(&r);
	  master_printf("Choosing color %d\n",r);
	  decripted=true;
	}
      
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
  
  constexpr bool decryptNameTest=false;
  if(decryptNameTest)
    for(const auto& [name,bw,fw] : mes2pts_contr_map)
      {
	int bwOrder=0;
	std::string line;
	for(const std::string& u : {bw,fw})
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
	master_printf("Decrypting %s = %s\n",name.c_str(),line.c_str());
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
  
  constexpr char DO_NOT_AVERAGE_HITS_STR[]="DO_NOT_AVERAGE_HITS";
  doNotAverageHits=getenv(DO_NOT_AVERAGE_HITS_STR)!=nullptr;
  if(doNotAverageHits)
    master_printf("%s exported, not averaging hits\n",DO_NOT_AVERAGE_HITS_STR);
  else
    master_printf("Averaging hits, export %s if needed otherwise\n",DO_NOT_AVERAGE_HITS_STR);
  
  //loop over the configs
  int iconf=0;
  while(read_conf_parameters(iconf,finish_file_present))
    {
      HitLooper hitLooper;
      
      for(int iHit=0;iHit<nhits;iHit++)
	{
	  hitLooper.start_hit(iHit);
	  hitLooper.run(iHit);
	  
	  compute_contractions(); //not working, here only to emit errors
	  propagators_fft(iHit); // same
	  
	  if(doNotAverageHits)
	    print_contractions();
	}
      
      master_printf("NOffloaded: %d\n",hitLooper.nOffloaded);
      master_printf("NRecalled: %d\n",hitLooper.nRecalled);
      
      free_confs();
      
      if(not doNotAverageHits)
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
