#include <nissa.hpp>

using namespace nissa;

//convention on gospel
const int follow_chris=0,follow_nazario=1;

//define types of quark propagator used
const int nins_kind=5;
enum insertion_t{                    ORIGINAL,  SCALAR,  PSEUDO,  PHOTON,   TADPOLE};//,  VECTOR};
const char ins_name[nins_kind][20]={"ORIGINAL","SCALAR","PSEUDO","PHOTON", "TADPOLE"};//, "VECTOR"};
const int nqprop_kind=6;
enum qprop_t{                           PROP_0,  PROP_S,  PROP_P,  PROP_T,  PROP_PHOTON,  PROP_PHOTON2};//,  PROP_VECTOR};
const char prop_name[nqprop_kind][20]={"PROP_0","PROP_S","PROP_P","PROP_T","PROP_PHOTON","PROP_PHOTON2"};//,"PROP_VECTOR"};

//map the source, the destination and the insertion for each propagator
const qprop_t prop_map[nqprop_kind]=         {PROP_0,   PROP_S, PROP_P, PROP_T,   PROP_PHOTON,  PROP_PHOTON2};//,  PROP_VECTOR};
const insertion_t insertion_map[nqprop_kind]={ORIGINAL, SCALAR, PSEUDO, TADPOLE,  PHOTON,       PHOTON};//,        VECTOR};
const qprop_t source_map[nqprop_kind]=       {PROP_0,   PROP_0, PROP_0, PROP_0,   PROP_0,       PROP_PHOTON};//,   PROP_0};
const char prop_abbr[]=                       "0"       "S"     "P"     "T"       "L"           "M";//            "V";

//sign of the lepton momentum
const int norie=2;
const int sign_orie[2]={-1,+1};

const int nins=3;
const int nlins=2;
const int nrev=2;

int nanalyzed_conf=0;
double tot_prog_time=0,wall_time;

//read the conf and setup it
void setup_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta,const char *conf_path,int rnd_gauge_transform,int free_theory)
{
  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  if(!free_theory)
    {
      read_ildg_gauge_conf(conf,conf_path);
      master_printf("plaq: %+016.016g\n",global_plaquette_lx_conf(conf));
    }
  else generate_cold_lx_conf(conf);
  
  //if asked, randomly transform the configurations
  if(rnd_gauge_transform) perform_random_gauge_transform(conf,conf);
  
  //put anti-periodic boundary condition for the fermionic propagator
  old_theta[0]=old_theta[1]=old_theta[2]=old_theta[3]=0;
  put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
  adapt_theta(conf,old_theta,put_theta,0,0);
}

//used to shift the configuration
void index_shift(int &irank_out,int &ivol_out,int ivol_in,void *pars)
{
  int *source_coord=(int*)pars;
  coords co;
  for(int nu=0;nu<NDIM;nu++) co[nu]=(glbCoordOfLoclx[ivol_in][nu]+source_coord[nu])%glbSize[nu];
  get_loclx_and_rank_of_coord(&ivol_out,&irank_out,co);
}

//perform a random shift
void random_shift_gauge_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta)
{
  //remove phase
  put_theta[0]=0;put_theta[1]=put_theta[2]=put_theta[3]=0;
  adapt_theta(conf,old_theta,put_theta,0,0);
  
  //source coord
  coords shift_coord;
  generate_random_coord(shift_coord);
  
  //shift the configuration
  double shift_time=-take_time();
  vector_remap_t shifter(locVol,index_shift,(void*)shift_coord);
  shifter.remap(conf,conf,sizeof(quad_su3));
  shift_time+=take_time();
  master_printf("Shifted of %d %d %d %d in %lg sec, plaquette after shift: %+016.016lg\n",shift_coord[0],shift_coord[1],shift_coord[2],shift_coord[3],shift_time,global_plaquette_lx_conf(conf));
  
  //put back the phase
  put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
  adapt_theta(conf,old_theta,put_theta,0,0);
}

//check if the time is enough
int check_remaining_time()
{
  int enough_time;
  
  //check remaining time
  double temp_time=take_time()+tot_prog_time;
  double ave_time=temp_time/nanalyzed_conf;
  double left_time=wall_time-temp_time;
  enough_time=left_time>(ave_time*1.1);
  
  master_printf("Remaining time: %lg sec\n",left_time);
  master_printf("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
  if(enough_time) master_printf("Continuing with next conf!\n");
  else master_printf("Not enough time, exiting!\n");
  
  return enough_time;
}
int map_mu[4]={4,1,2,3};

#ifdef POINT_SOURCE_VERSION
 #define PROP_TYPE su3spinspin
#else
 #define PROP_TYPE colorspinspin
#endif

/* the loop is normalised such that the physical rate at leading order
  is obtained multiplying the loop by Gf^2 fpi^2 * phi2 (space phase
  factor) which is (1-rl^2)/(16 pi mpi) where rl=ml/mpi, whereas the
  interference is obtained by the full hadrolepton correlation
  multiplied by 4*mpi*fpi*Gf^2*phi2 */

/////////////////////////////////////// data //////////////////////////////

int loc_pion_curr;
int loc_muon_curr;

int ninv_tot=0,nhadr_contr_tot=0,nlept_contr_tot=0,nsource_tot=0,nphoton_prop_tot=0;
double inv_time=0,hadr_contr_time=0,lept_contr_time=0,print_time=0;
double source_time=0,photon_prop_time=0,lepton_prop_time=0;

int follow_chris_or_nazario;
int free_theory,rnd_gauge_transform;
int ngauge_conf;
char conf_path[1024],outfolder[1024];
quad_su3 *conf;

tm_basis_t base;
int pure_wilson;
double kappa;
double put_theta[4],old_theta[4]={0,0,0,0};

PROP_TYPE *source,*original_source;
int seed,noise_type;

int nqmass,nr,nsources;
double *qmass,*qkappa,*residue;
namespace{
  PROP_TYPE **Q;
}
spincolor *temp_source;
spincolor *temp_solution;

gauge_info photon;
double tadpole[4];
spin1field *photon_field;

int hadr_corr_length;
complex *hadr_corr;
const int nhadr_contr=16+12;
const int ig_hadr_so[nhadr_contr]={ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5  , 1, 2, 3, 10,11,12, 10,11,12,13,14,15};
const int ig_hadr_si[nhadr_contr]={ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15  , 10,11,12,1, 2, 3,  10,11,12,13,14,15};
complex *glb_corr,*loc_corr;

//list the 8 matrices to insert for the weak current
const int nweak_ins=17;
const int nweak_ind=9;
//const int nhadrolept_proj=4,hadrolept_projs[nhadrolept_proj]={9,4,5,0};
const int nhadrolept_proj=1,hadrolept_projs[nhadrolept_proj]={4};
const int list_weak_insq[nweak_ins]=     {1,2,3,4, 6,7,8,9,  1,2,3,4, 6,7,8,9, 5};
const int list_weak_insl[nweak_ins]=     {1,2,3,4, 6,7,8,9,  6,7,8,9, 1,2,3,4, 5};
const int list_weak_ind_contr[nweak_ins]={0,0,0,1, 2,2,2,3,  4,4,4,5, 6,6,6,7, 8};
const char list_weak_ind_nameq[nweak_ind][3]={"VK","V0","AK","A0","VK","V0","AK","A0","P5"};
const char list_weak_ind_namel[nweak_ind][3]={"VK","V0","AK","A0","AK","A0","VK","V0","V0"};
int nind;
spinspin *hadr;
complex *hadrolept_corr;

//hadron contractions
const int ncombo_hadr_corr=6;
const qprop_t prop1_hadr_map[ncombo_hadr_corr]={PROP_0,PROP_0,PROP_0,PROP_0,PROP_0      ,PROP_PHOTON};//,PROP_0};
const qprop_t prop2_hadr_map[ncombo_hadr_corr]={PROP_0,PROP_S,PROP_P,PROP_T,PROP_PHOTON2,PROP_PHOTON};//,PROP_VECTOR};

//parameters of the leptons
int nleptons;
int *lep_corr_iq1;
int *lep_corr_iq2;
tm_quark_info *leps;
double *lep_energy,*neu_energy;
spinspin **L,*temp_lep;

//prototype
void generate_original_source();
void generate_photon_stochastic_propagator();

////////////////////////////////////////////// get propagator and prop. info /////////////////////////////////////////////

//return appropriate propagator
int nqprop,nlprop;
int iqprop(int imass,int ip,int r)
{return r+nr*(imass+nqmass*ip);}
int ilprop(int ilepton,int ilins,int orie,int r)
{return r+nr*(ilins+nlins*(orie+norie*ilepton));}

//return appropriately modified info
tm_quark_info get_lepton_info(int ilepton,int orie,int r)
{
  tm_quark_info le=leps[ilepton];
  le.r=r;
  for(int i=1;i<NDIM;i++) le.bc[i]*=sign_orie[orie];
  
  return le;
}

///////////////////////////////// initialise the library, read input file, allocate /////////////////////////////////////

void init_simulation(char *path)
{
  //open input file
  open_input(path);
  
  //init the grid
  {
    int L,T;
    read_str_int("L",&L);
    read_str_int("T",&T);
    //Init the MPI grid
    init_grid(T,L);
  }
  
  //Wall time
  read_str_double("WallTime",&wall_time);
  //Pure Wilson
  read_str_int("PureWilson",&pure_wilson);
  if(pure_wilson)
    {
      base=WILSON_BASE;
      nr=1;
      read_list_of_double_pairs("QKappaResidues",&nqmass,&qkappa,&residue);
    }
  else
    {
      base=MAX_TWIST_BASE;
      //Kappa
      read_str_double("Kappa",&kappa);
      //One or two r
      read_str_int("NR",&nr);
      //Masses and residue
      read_list_of_double_pairs("QMassResidues",&nqmass,&qmass,&residue);
    }
  
  //Leptons
  read_str_int("LeptonicCorrs",&nleptons);
  lep_corr_iq1=nissa_malloc("lep_corr_iq1",nleptons,int);
  lep_corr_iq2=nissa_malloc("lep_corr_iq2",nleptons,int);
  leps=nissa_malloc("leps",nleptons,tm_quark_info);
  lep_energy=nissa_malloc("lep_energy",nleptons,double);
  neu_energy=nissa_malloc("neu_energy",nleptons,double);
  if(!pure_wilson) expect_str("Q1Q2LepmassMesmass");
  else             expect_str("Q1Q2LepkappaMesmass");
  for(int il=0;il<nleptons;il++)
    {
      //read quarks identfiying the mesons
      read_int(lep_corr_iq1+il);
      read_int(lep_corr_iq2+il);
      
      //if not pure wilson read mass
      if(pure_wilson) leps[il].mass=0;
      else            read_double(&leps[il].mass);
      
      //antiperiodic
      leps[il].bc[0]=1;
      
      //maximal twist (if tm), otherwise read kappa
      if(pure_wilson) read_double(&leps[il].kappa);
      else            leps[il].kappa=0.125;
      leps[il].r=0;
      
      //read the mass of the meson (that must have been determined outside)
      double mes_mass;
      read_double(&mes_mass);
      
      //set initial value of bc and check kinematic
      for(int i=1;i<NDIM;i++) leps[il].bc[i]=0;
      if(tm_quark_energy(leps[il],0)>=mes_mass) crash("initial state is lighter (%lg) than final state at rest (%lg)!",mes_mass,tm_quark_energy(leps[il],0));
      
      //compute meson momentum and bc
      double err;
      do
      	{
      	  //compute the error
	  double lep_energy=tm_quark_energy(leps[il],0);
	  double neu_energy=naive_massless_quark_energy(leps[il].bc,0);
      	  err=lep_energy+neu_energy-mes_mass;
      	  //compute the derivative
      	  double eps=1e-8;
      	  for(int i=1;i<NDIM;i++) leps[il].bc[i]+=eps;
      	  double der=(tm_quark_energy(leps[il],0)+naive_massless_quark_energy(leps[il].bc,0)-mes_mass-err)/eps;
      	  for(int i=1;i<NDIM;i++) leps[il].bc[i]-=eps+err/der;
	  
      	  master_printf("lep_e: %+010.10lg, neu_e: %+010.10lg, mes_mass: %lg, error: %lg, der: %lg\n",lep_energy,neu_energy,mes_mass,err,der);
      	}
      while(fabs(err)>1e-14);
      
      //write down energy
      lep_energy[il]=tm_quark_energy(leps[il],0);
      neu_energy[il]=naive_massless_quark_energy(leps[il].bc,0);
      master_printf(" ilepton %d, lepton energy: %lg, neutrino energy: %lg\n",il,lep_energy[il],neu_energy[il]);
      master_printf(" lep+neut energy: %lg\n",lep_energy[il]+neu_energy[il]);
      master_printf(" bc: %+016.016lg\n",leps[il].bc[1]);
    }
  
  //Zero mode subtraction
  char zero_mode_sub_str[100];
  read_str_str("ZeroModeSubtraction",zero_mode_sub_str,100);
  
  if(strncasecmp(zero_mode_sub_str,"PECIONA",100)==0) photon.zms=PECIONA;
  else
    if(strncasecmp(zero_mode_sub_str,"UNNO_ALEMANNA",100)==0) photon.zms=UNNO_ALEMANNA;
    else
      if(strncasecmp(zero_mode_sub_str,"ONLY_100",100)==0) photon.zms=ONLY_100;
      else crash("Unkwnown zero mode subtraction: %s",zero_mode_sub_str);
  
  //gauge for photon propagator
  char photon_gauge_str[100];
  read_str_str("PhotonGauge",photon_gauge_str,100);
  if(strncasecmp(photon_gauge_str,"FEYNMAN",100)==0) photon.alpha=FEYNMAN_ALPHA;
  else
    if(strncasecmp(photon_gauge_str,"LANDAU",100)==0) photon.alpha=LANDAU_ALPHA;
    else
      if(strncasecmp(photon_gauge_str,"LANDAU",100)==0) read_str_double("Alpha",&photon.alpha);
      else crash("Unkwnown photon gauge: %s",photon_gauge_str);
  
  //discretization for photon propagator
  char photon_discrete_str[100];
  read_str_str("PhotonDiscretization",photon_discrete_str,100);
  if(strncasecmp(photon_discrete_str,"WILSON",100)==0) photon.c1=WILSON_C1;
  else
    if(strncasecmp(photon_discrete_str,"TLSYM",100)==0) photon.c1=TLSYM_C1;
    else crash("Unkwnown photon discretization: %s",photon_discrete_str);
  
  //initialize the random generator with the read seed
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  //noise type
  read_str_int("NoiseType",&noise_type);
  
  //flag to simulate in the free theory
  read_str_int("FreeTheory",&free_theory);
  
  //flag to make the muon with or without the external line
  read_str_int("FollowChrisOrNazario",&follow_chris_or_nazario);
  
  //perform a random gauge transformation
  read_str_int("RandomGaugeTransform",&rnd_gauge_transform);
  
  //local current on muon or pion
  read_str_int("LocPionCurr",&loc_pion_curr);
  read_str_int("LocMuonCurr",&loc_muon_curr);
  
  //number of sources
  read_str_int("NSources",&nsources);
  
  //number of configurations
  read_str_int("NGaugeConf",&ngauge_conf);
  
  ///////////////////// finihed reading apart from conf list ///////////////
  
  //compute the tadpole summing all momentum
  compute_tadpole(tadpole,photon);
  
  //Allocate
  nqprop=iqprop(nqmass-1,nqprop_kind-1,nr-1)+1;
  nlprop=ilprop(nleptons-1,nlins-1,norie-1,nr-1)+1;
  
  //allocate temporary vectors
  temp_source=nissa_malloc("temp_source",locVol,spincolor);
  temp_solution=nissa_malloc("temp_solution",locVol+bord_vol,spincolor);
  hadr_corr_length=glbSize[0]*nhadr_contr*ncombo_hadr_corr*nqmass*nqmass*nr;
  hadr_corr=nissa_malloc("hadr_corr",hadr_corr_length,complex);
  glb_corr=nissa_malloc("glb_corr",glbSize[0]*nhadr_contr,complex);
  loc_corr=nissa_malloc("loc_corr",glbSize[0]*nhadr_contr,complex);
  nind=nleptons*nweak_ind*norie*nr*nins;
  hadr=nissa_malloc("hadr",locVol,spinspin);
  hadrolept_corr=nissa_malloc("hadrolept_corr",glbSize[0]*nweak_ind*nhadrolept_proj*nind,complex);
  original_source=nissa_malloc("source",locVol,PROP_TYPE);
  source=nissa_malloc("source",locVol,PROP_TYPE);
  photon_field=nissa_malloc("photon_phield",locVol+bord_vol,spin1field);
  Q=nissa_malloc("Q*",nqprop,PROP_TYPE*);
  for(int iprop=0;iprop<nqprop;iprop++) Q[iprop]=nissa_malloc("Q",locVol+bord_vol,PROP_TYPE);
  L=nissa_malloc("L*",nlprop,spinspin*);
  for(int iprop=0;iprop<nlprop;iprop++) L[iprop]=nissa_malloc("L",locVol+bord_vol,spinspin);
  temp_lep=nissa_malloc("temp_lep",locVol+bord_vol,spinspin);
  conf=nissa_malloc("conf",locVol+bord_vol,quad_su3);
}

//find a new conf
int read_conf_parameters(int &iconf)
{
  int ok_conf;
  
  do
    {
      //Gauge path
      read_str(conf_path,1024);
      
      //Out folder
      read_str(outfolder,1024);
      
      //Check if the conf has been finished or is already running
      master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
      char fin_file[1024],run_file[1024];
      safe_snprintf(fin_file,1024,"%s/%s",outfolder, finished_filename.c_str());
      safe_snprintf(run_file,1024,"%s/%s",outfolder, running_filename.c_str());
      ok_conf=!(file_exists(fin_file)) && !(file_exists(run_file));
      
      //if not finished
      if(ok_conf)
	{
	  master_printf(" Configuration \"%s\" not yet analyzed, starting",conf_path);
	  if(!dir_exists(outfolder))
	    {
	      int ris=create_dir(outfolder);
	      if(ris==0) master_printf(" Output path \"%s\" not present, created.\n",outfolder);
	      else
		crash(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
	    }
	  file_touch(run_file);
	}
      else
	{
	  master_printf(" In output path \"%s\" terminating file already present: configuration \"%s\" already analyzed, skipping.\n",outfolder,conf_path);
	  for(int isource=0;isource<nsources;isource++)
	    {
	      coords coord;
	      generate_random_coord(coord);
	      generate_stochastic_tlSym_gauge_propagator_source(photon_field);
	      generate_original_source();
	    }
	}
      iconf++;
    }
  while(!ok_conf && iconf<ngauge_conf);
  
  master_printf("\n");
  
  return ok_conf;
}

//init a new conf
void start_new_conf()
{
  setup_conf(conf,old_theta,put_theta,conf_path,rnd_gauge_transform,free_theory);
  
  //reset correlations
  vector_reset(hadr_corr);
  vector_reset(hadrolept_corr);
}

//generate a wall-source for stochastic QCD propagator
void generate_original_source()
{
  //reset the real source position
  coords origin_coord;
  for(int mu=0;mu<NDIM;mu++) origin_coord[mu]=0;
  
#ifdef POINT_SOURCE_VERSION
  master_printf("Source position: t=%d x=%d y=%d z=%d\n",origin_coord[0],origin_coord[1],origin_coord[2],origin_coord[3]);
  generate_delta_source(original_source,origin_coord);
#else
  enum rnd_t rnd_type_map[6]={RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_Z2,RND_Z2,RND_Z4,RND_GAUSS};
  generate_spindiluted_source(original_source,rnd_type_map[noise_type],origin_coord[0]);
#endif
}

//////////////////////////////////////// quark propagators /////////////////////////////////////////////////

//insert the photon on the source side
void insert_external_loc_source(PROP_TYPE *out,spin1field *curr,bool *dirs,PROP_TYPE *in,int t)
{ 
  
  if(in==out) crash("in==out");
  
  vector_reset(out);
  
  for(int mu=0;mu<NDIM;mu++)
    if(dirs[mu])
      NISSA_PARALLEL_LOOP(ivol,0,locVol)
	if(t==-1||glbCoordOfLoclx[ivol][0]==t)
	  {
	    PROP_TYPE temp1,temp2;
	    NAME2(unsafe_dirac_prod,PROP_TYPE)(temp1,base_gamma+map_mu[mu],in[ivol]);
	    NAME3(unsafe,PROP_TYPE,prod_complex)(temp2,temp1,curr[ivol][mu]);
	    NAME2(PROP_TYPE,summ_the_prod_idouble)(out[ivol],temp2,1);
	  }
  NISSA_PARALLEL_LOOP_END;
  
  set_borders_invalid(out);
}
//insert the photon on the source
void insert_external_loc_source(PROP_TYPE *out,spin1field *curr,PROP_TYPE *in,int t)
{insert_external_loc_source(out,curr,all_dirs,in,t);}

void insert_external_source(PROP_TYPE *out,spin1field *curr,PROP_TYPE *ori,int t,int r,int loc)
{
  if(loc) insert_external_loc_source(source,curr,ori,t);
  else
    if(!pure_wilson) insert_tm_external_source(source,conf,curr,ori,r,all_dirs,t);
    else             insert_Wilson_external_source(source,conf,curr,ori,all_dirs,t);
}

//generate a sequential source
void generate_source(insertion_t inser,int r,PROP_TYPE *ori,int t=-1)
{
  source_time-=take_time();
  
  switch(inser)
    {
    case ORIGINAL:prop_multiply_with_gamma(source,0,original_source);break;
    case SCALAR:prop_multiply_with_gamma(source,0,ori);break;
    case PSEUDO:prop_multiply_with_gamma(source,5,ori);break;
    case PHOTON:insert_external_source(source,photon_field,ori,t,r,loc_pion_curr);break;
    case TADPOLE:
      if(!pure_wilson) insert_tm_tadpole(source,conf,ori,r,tadpole,-1);
      else             insert_Wilson_tadpole(source,conf,ori,tadpole,-1);
      break;
      //case VECTOR:insert_external_source(source,NULL,ori,t,r,loc_pion_curr);break;
    }
  
  source_time+=take_time();
  nsource_tot++;
}

//invert on top of a source, putting all needed for the appropriate quark
void get_qprop(PROP_TYPE *out,PROP_TYPE *in,int imass,bool r)
{
  //these are the ways in which Dirac operator rotates - propagator is opposite, see below
#ifdef POINT_SOURCE_VERSION
  for(int ic=0;ic<NCOL;ic++)
#endif
    for(int id=0;id<4;id++)
      { 
	//read the source out
#ifdef POINT_SOURCE_VERSION
	get_spincolor_from_su3spinspin(temp_source,in,id,ic);
#else
	get_spincolor_from_colorspinspin(temp_source,in,id);
#endif
	
	//rotate the source index - the propagator rotate AS the sign of mass term
	if(!pure_wilson) safe_dirac_prod_spincolor(temp_source,(tau3[r]==-1)?&Pminus:&Pplus,temp_source);
	
	//invert
	inv_time-=take_time();
	if(!pure_wilson) inv_tmD_cg_eoprec(temp_solution,NULL,conf,kappa,tau3[r]*qmass[imass],100000,residue[imass],temp_source);
	else             inv_tmD_cg_eoprec(temp_solution,NULL,conf,qkappa[imass],0,100000,residue[imass],temp_source);
	ninv_tot++;inv_time+=take_time();
	
	//rotate the sink index
	if(!pure_wilson) safe_dirac_prod_spincolor(temp_solution,(tau3[r]==-1)?&Pminus:&Pplus,temp_solution);
	
	//put the output on place
#ifdef POINT_SOURCE_VERSION
	master_printf("  finished the inversion dirac index %d, color %d\n",id,ic);
	put_spincolor_into_su3spinspin(out,temp_solution,id,ic);
#else
	master_printf("  finished the inversion dirac index %d\n",id);
	put_spincolor_into_colorspinspin(out,temp_solution,id);
#endif
      }
}

//generate all the quark propagators
void generate_quark_propagators()
{
  for(int ip=0;ip<nqprop_kind;ip++)
    {
      master_printf("Generating propagtor of type %s inserting %s on source %s\n",prop_name[prop_map[ip]],ins_name[insertion_map[ip]],prop_name[source_map[ip]]);
      for(int imass=0;imass<nqmass;imass++)
	for(int r=0;r<nr;r++)
	  {
	    if(!pure_wilson) master_printf(" mass[%d]=%lg, r=%d\n",imass,qmass[imass],r);
	    else             master_printf(" kappa[%d]=%lg\n",imass,qkappa[imass]);
	    generate_source(insertion_map[ip],r,Q[iqprop(imass,source_map[ip],r)]);
	    get_qprop(Q[iqprop(imass,prop_map[ip],r)],source,imass,r);
	  }
    }
}

/////////////////////////////////////////////// photon propagators ///////////////////////////////////////////

//wrapper to generate a stochastic propagator
void generate_photon_stochastic_propagator()
{
  photon_prop_time-=take_time();
  generate_stochastic_tlSym_gauge_propagator_source(photon_field);
  
  //do also the alternative version
  multiply_by_sqrt_tlSym_gauge_propagator(photon_field,photon_field,photon);
  
  photon_prop_time+=take_time();
  nphoton_prop_tot++;
}

/////////////////////////////////////////////// lepton propagators ///////////////////////////////////////////

//compute phase exponent for space part: vec{p}*\vec{x}
double get_space_arg(const int ivol,const momentum_t bc)
{
  double arg=0;
  for(int mu=1;mu<NDIM;mu++)
    {
      double step=bc[mu]*M_PI/glbSize[mu];
      arg+=step*glbCoordOfLoclx[ivol][mu];
    }
  return arg;
}

//compute the phase for lepton on its sink
void get_lepton_sink_phase_factor(complex out,int ivol,int ilepton,tm_quark_info le)
{
  //compute space and time factor
  double arg=get_space_arg(ivol,le.bc);
  int t=glbCoordOfLoclx[ivol][0];
  if(follow_chris_or_nazario==follow_nazario && t>=glbSize[0]/2) t=glbSize[0]-t;
  double ext=exp(t*lep_energy[ilepton]);
  
  //compute full exponential (notice the factor -1)
  out[RE]=cos(-arg)*ext;
  out[IM]=sin(-arg)*ext;
}

//compute the phase for antineutrino - the orientation is that of the muon (as above)
void get_antineutrino_source_phase_factor(complex out,const int ivol,const int ilepton,const momentum_t bc)
{
  //compute space and time factor
  double arg=get_space_arg(ivol,bc);
  int t=glbCoordOfLoclx[ivol][0];
  if(follow_chris_or_nazario==follow_nazario && t>=glbSize[0]/2) t=glbSize[0]-t;
  double ext=exp(t*neu_energy[ilepton]);
  
  //compute full exponential (notice the factor +1)
  out[RE]=cos(+arg)*ext;
  out[IM]=sin(+arg)*ext;
}

//set everything to a phase factor
void set_to_lepton_sink_phase_factor(spinspin *prop,int ilepton,tm_quark_info &le)
{
  
  vector_reset(prop);
  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    {
      complex ph;
      get_lepton_sink_phase_factor(ph,ivol,ilepton,le);
      spinspin_put_to_diag(prop[ivol],ph);
    }
  NISSA_PARALLEL_LOOP_END;
  set_borders_invalid(prop);
}

//insert the photon on the source side
void insert_photon_on_the_source(spinspin* prop,spin1field* A,bool* dirs,tm_quark_info le,int twall)
{
  
  //select A
  communicate_lx_spin1field_borders(A);
  
  //copy on the temporary and communicate borders
  vector_copy(temp_lep,prop);
  communicate_lx_spinspin_borders(temp_lep);
  vector_reset(prop);
  
  if(!loc_muon_curr)
    {
      dirac_matr GAMMA;
      if(pure_wilson) dirac_prod_double(&GAMMA,base_gamma+0,1);
      else dirac_prod_idouble(&GAMMA,base_gamma+5,-tau3[le.r]);
      
      //phases
      quad_u1 phases;
      for(int mu=0;mu<NDIM;mu++)
	{
	  phases[mu][0]=cos(le.bc[mu]*M_PI);
	  phases[mu][1]=sin(le.bc[mu]*M_PI);
	}
      
      //prepare each propagator for a single lepton
      //by computing i(phi(x-mu)A_mu(x-mu)(-i t3 g5-gmu)/2-phi(x+mu)A_mu(x)(-i t3 g5+gmu)/2)=
      //(ph0 A_mu(x-mu)g[r][0][mu]-ph0 A_mu(x)g[r][1][mu])=
      for(int mu=0;mu<NDIM;mu++)
	if(dirs[mu])
	  NISSA_PARALLEL_LOOP(ivol,0,locVol)
	    if(twall==-1||glbCoordOfLoclx[ivol][0]==twall)
	      {
		//find neighbors
		int ifw=loclxNeighup[ivol][mu];
		int ibw=loclxNeighdw[ivol][mu];
		
		//compute phase factor
		spinspin ph_bw,ph_fw;
		
		//transport down and up
		if(glbCoordOfLoclx[ivol][mu]==glbSize[mu]-1) unsafe_spinspin_prod_complex_conj2(ph_fw,temp_lep[ifw],phases[mu]);
		else spinspin_copy(ph_fw,temp_lep[ifw]);
		if(glbCoordOfLoclx[ivol][mu]==0) unsafe_spinspin_prod_complex(ph_bw,temp_lep[ibw],phases[mu]);
		else spinspin_copy(ph_bw,temp_lep[ibw]);
		
		//fix coefficients - i is inserted here!
		//also dir selection is made here
		spinspin_prodassign_idouble(ph_fw,-0.5*dirs[mu]);
		spinspin_prodassign_idouble(ph_bw,+0.5*dirs[mu]);
		
		//fix insertion of the current
		safe_spinspin_prod_complex(ph_fw,ph_fw,A[ivol][mu]);
		safe_spinspin_prod_complex(ph_bw,ph_bw,A[ibw][mu]);
		
		//summ and subtract the two
		spinspin fw_M_bw,fw_P_bw;
		spinspin_subt(fw_M_bw,ph_fw,ph_bw);
		spinspin_summ(fw_P_bw,ph_fw,ph_bw);
		
		//put GAMMA on the summ
		spinspin temp_P;
		unsafe_spinspin_prod_dirac(temp_P,fw_P_bw,&GAMMA);
		spinspin_summassign(prop[ivol],temp_P);
		
		//put gmu on the diff
		spinspin temp_M;
		unsafe_spinspin_prod_dirac(temp_M,fw_M_bw,base_gamma+map_mu[mu]);
		spinspin_summassign(prop[ivol],temp_M);
	      }
      NISSA_PARALLEL_LOOP_END;
    }
  else
    {
      for(int mu=0;mu<NDIM;mu++)
	if(dirs[mu])
	  NISSA_PARALLEL_LOOP(ivol,0,locVol)
	    if(twall==-1||glbCoordOfLoclx[ivol][0]==twall)
	      {
		spinspin temp1,temp2;
		unsafe_spinspin_prod_dirac(temp1,temp_lep[ivol],base_gamma+map_mu[mu]);
		unsafe_spinspin_prod_complex(temp2,temp1,A[ivol][mu]);
		spinspin_summ_the_prod_idouble(prop[ivol],temp2,1);
	      }
     NISSA_PARALLEL_LOOP_END; 
    }
  
  set_borders_invalid(prop);
}

void insert_photon_on_the_source(spinspin *prop,bool *dirs,tm_quark_info le,int twall)
{
  if(!loc_muon_curr) master_printf("Inserting photon point-split on time %d\n",twall);
  else master_printf("Inserting photon locally on time %d\n");
  insert_photon_on_the_source(prop,photon_field,dirs,le,twall);
}

//insert the photon on the source
void insert_photon_on_the_source(spinspin *prop,tm_quark_info &le,int twall)
{insert_photon_on_the_source(prop,all_dirs,le,twall);}

//generate all the lepton propagators, pointing outward
//the computations is done by:
// 1)putting the correct phase in x space, given by exp(E_mu*t-i*vec(p)*vec(x))
// 2)multiplying it by the conserved current inserting eta or phi
// 3)going to momentum space
// 4)multiplying by the lepton propagator in momentum space
// 5)coming back to x space
void generate_lepton_propagators()
{
  
  if(IS_MASTER_THREAD) lepton_prop_time-=take_time();
  
  for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int ilins=0;ilins<nlins;ilins++)
      for(int ori=0;ori<norie;ori++)
	for(int r=0;r<nr;r++)
	  {
	    //set the properties of the meson
	    //time boundaries are anti-periodic, space are as for external line
	    tm_quark_info le=get_lepton_info(ilepton,ori,r);
	    
	    //select the propagator
	    int iprop=ilprop(ilepton,ilins,ori,r);
	    spinspin *prop=L[iprop];
	    
	    //put it to a phase
	    set_to_lepton_sink_phase_factor(prop,ilepton,le);
	    
	    //if we are doing Nazario's way (with the external line) add it
	    if(follow_chris_or_nazario==follow_nazario)
	      {
		//select only the wall
		int tmiddle=glbSize[0]/2;
		select_propagator_timeslice(prop,prop,tmiddle);
		multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le,base,true);
	      }
	    
	    //insert or not photon
	    if(ilins)
	      {
		//insert photon and prolong
		insert_photon_on_the_source(prop,le,-1); //all times
		multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le,base,true);
	      }
	  }
  
  if(IS_MASTER_THREAD) lepton_prop_time+=take_time();
}

////////////////////////////////////////// purely hadronic correlators ///////////////////////////////////////////

//compute all the hadronic correlations
void compute_hadronic_correlations()
{
  master_printf("Computing hadronic correlation functions\n");
  
  hadr_contr_time-=take_time();
  for(int icombo=0;icombo<ncombo_hadr_corr;icombo++)
    for(int imass=0;imass<nqmass;imass++)
      for(int jmass=0;jmass<nqmass;jmass++)
	for(int r=0;r<nr;r++)
	  {
	    //compute the correlation function
	    int ip1=iqprop(imass,prop1_hadr_map[icombo],r);
	    int ip2=iqprop(jmass,prop2_hadr_map[icombo],r);
	    std::swap(ip1,ip2);
	    meson_two_points_Wilson_prop(glb_corr,loc_corr,ig_hadr_so,Q[ip1],ig_hadr_si,Q[ip2],nhadr_contr);
	    nhadr_contr_tot+=nhadr_contr;
	    
	    //save to the total stack
	    for(int ihadr_contr=0;ihadr_contr<nhadr_contr;ihadr_contr++)
	      for(int t=0;t<glbSize[0];t++)
		{
		  int i=t+glbSize[0]*(ihadr_contr+nhadr_contr*(r+nr*(jmass+nqmass*(imass+nqmass*icombo))));
		  complex_summassign(hadr_corr[i],glb_corr[t+glbSize[0]*ihadr_contr]);
		}
	  }
  hadr_contr_time+=take_time();
}

/////////////////////////////////////////////// hadroleptonic correlators //////////////////////////////////////////

//compute the hadronic part of the lepton correlation function
//as usual, FIRST propagator is reverted
void hadronic_part_leptonic_correlation(spinspin* hadr,PROP_TYPE* S1,PROP_TYPE* S2)
{
  
  vector_reset(hadr);
  
  //it's just the matter of inserting gamma5*gamma5=identity between S1^dag and S2
  //on the sink gamma5 must still be inserted!
  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    for(int ic_si=0;ic_si<NCOL;ic_si++)
#ifdef POINT_SOURCE_VERSION
      for(int ic_so=0;ic_so<NCOL;ic_so++)
#endif
	for(int id_si1=0;id_si1<4;id_si1++)
	  for(int id_si2=0;id_si2<4;id_si2++)
	    for(int id_so=0;id_so<4;id_so++)
	      complex_summ_the_conj1_prod
		(hadr[ivol][id_si2][id_si1], //this way when taking the trace with dirac matrix, that is acting on S2, as it should
#ifdef POINT_SOURCE_VERSION
		 S1[ivol][ic_si][ic_so][id_si1][id_so],S2[ivol][ic_si][ic_so][id_si2][id_so])
#else
                 S1[ivol][ic_si][id_si1][id_so],S2[ivol][ic_si][id_si2][id_so])
#endif
		 ;
   NISSA_PARALLEL_LOOP_END;
   THREAD_BARRIER();
}

//compute the leptonic part of the correlation function
void attach_leptonic_correlation(spinspin* hadr,int iprop,int ilepton,int orie,int rl,int ext_ind)
{
  
  vector_reset(loc_corr);
  
  //get the lepton info and prop
  tm_quark_info le=get_lepton_info(ilepton,orie,rl);
  spinspin *lept=L[iprop];
  
  //get the projectors
  spinspin promu[2],pronu[2];
  twisted_on_shell_operator_of_imom(promu[0],le,0,false,-1,base);
  if(follow_chris_or_nazario==follow_nazario) twisted_on_shell_operator_of_imom(promu[1],le,0,false,+1,base);
  else twisted_on_shell_operator_of_imom(promu[1],le,0,false,-1,base);
  naive_massless_on_shell_operator_of_imom(pronu[0],le.bc,0,-1);
  if(follow_chris_or_nazario==follow_nazario) naive_massless_on_shell_operator_of_imom(pronu[1],le.bc,0,+1);
  else naive_massless_on_shell_operator_of_imom(pronu[1],le.bc,0,-1);
  if(follow_chris_or_nazario==follow_chris)
    for(int i=0;i<2;i++)
      safe_spinspin_prod_dirac(promu[i],promu[i],base_gamma+map_mu[0]);
  
  //compute the right part of the leptonic loop: G0 G^dag
  dirac_matr hadrolept_proj_gamma[nhadrolept_proj];
  for(int ig_proj=0;ig_proj<nhadrolept_proj;ig_proj++)
    {
      int ig=hadrolept_projs[ig_proj];
      dirac_matr temp_gamma;
      dirac_herm(&temp_gamma,base_gamma+ig);
      dirac_prod(hadrolept_proj_gamma+ig_proj,base_gamma+map_mu[0],&temp_gamma);
    }
  //insert gamma5 on the sink-hadron-gamma: S1^dag G5 GW S2 (G5 G5) - no dagger, no commutator because it's on the LO leptonic part
  dirac_matr weak_ins_hadr_gamma[nweak_ins];
  for(int ins=0;ins<nweak_ins;ins++) dirac_prod(weak_ins_hadr_gamma+ins,base_gamma+5,base_gamma+list_weak_insq[ins]);
  
  //define the combined weak projectors (see below)
  dirac_matr neutr_1m_g5_proj;
  dirac_subt(&neutr_1m_g5_proj,base_gamma+0,base_gamma+5);
  
  for(int ins=0;ins<nweak_ins;ins++)
    {
      //define a local storage
      spinspin hl_loc_corr[locSize[0]];
      for(int i=0;i<locSize[0];i++) spinspin_put_to_zero(hl_loc_corr[i]);
      
      NISSA_PARALLEL_LOOP(ivol,0,locVol)
	{
	  int t=locCoordOfLoclx[ivol][0];
	  
	  //multiply lepton side on the right (source) side
	  spinspin la;
	  unsafe_spinspin_prod_dirac(la,lept[ivol],base_gamma+list_weak_insl[ins]);
	  
	  //include 4*(1-5)/2/2=(1-5) coming from the two neturino projector+(1-g5) weak lepton structure
	  //the second /2 comes from sqr(1/sqrt(2)) of 1502.00257
	  spinspin l;
	  unsafe_spinspin_prod_dirac(l,la,&neutr_1m_g5_proj);
	  
	  //get the neutrino phase (multiply hadron side) - notice that the sign of momentum is internally reversed
	  complex ph;
	  get_antineutrino_source_phase_factor(ph,ivol,ilepton,le.bc);
	  
	  //trace hadron side
	  complex h;
	  trace_spinspin_with_dirac(h,hadr[ivol],weak_ins_hadr_gamma+ins);
	  
	  //combine hl
	  complex_prodassign(h,ph);
	  spinspin_summ_the_complex_prod(hl_loc_corr[t],l,h);
	}
      NISSA_PARALLEL_LOOP_END;
      glb_threads_reduce_double_vect((double*)hl_loc_corr,locSize[0]*sizeof(spinspin)/sizeof(double));
      
      //save projection on LO
      for(int ig_proj=0;ig_proj<nhadrolept_proj;ig_proj++)
	NISSA_PARALLEL_LOOP(loc_t,0,locSize[0])
	  {
	    int glb_t=loc_t+rank_coord[0]*locSize[0];
	    int ilnp=(glb_t>=glbSize[0]/2); //select the lepton/neutrino projector
	    
	    spinspin td;
	    unsafe_spinspin_prod_spinspin(td,hl_loc_corr[loc_t],pronu[ilnp]);
	    spinspin dtd;
	    unsafe_spinspin_prod_spinspin(dtd,promu[ilnp],td);
	    complex hl;
	    trace_spinspin_with_dirac(hl,dtd,hadrolept_proj_gamma+ig_proj);
	    
	    //summ the average
	    int i=glb_t+glbSize[0]*(ig_proj+nhadrolept_proj*(list_weak_ind_contr[ins]+nweak_ind*ext_ind));
	    complex_summ_the_prod_double(hadrolept_corr[i],hl,1.0/glbSpatVol); //here to remove the statistical average on xw
	  }
      NISSA_PARALLEL_LOOP_END;
      if(IS_MASTER_THREAD) nlept_contr_tot+=nhadrolept_proj;
      THREAD_BARRIER();
    }
}

//return the index of the combination of r, orientation, etc
int hadrolept_corrpack_ind(int rl,int orie,int r2,int irev,int qins,int ilepton)
{return rl+nr*(orie+norie*(r2+nr*(irev+nrev*(qins+nins*ilepton))));}

//compute the total hadroleptonic correlation functions
void compute_hadroleptonic_correlations()
{
  master_printf("Computing leptonic correlation functions\n");
  lept_contr_time-=take_time();
  
 for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int qins=0;qins<nins;qins++)
      for(int irev=0;irev<nrev;irev++)
	for(int r2=0;r2<nr;r2++)
	  {
	    //takes the index of the quarks
	    int iq1=lep_corr_iq1[ilepton];
	    int iq2=lep_corr_iq2[ilepton];
	      
	    //takes the propagators
	    qprop_t PROP1_TYPE=PROP_0,PROP2_TYPE=PROP_0;
	    if(qins==1) PROP1_TYPE=PROP_PHOTON;
	    if(qins==2) PROP2_TYPE=PROP_PHOTON;
	    
	    int ip1=iqprop(iq1,PROP1_TYPE,r2);
	    int ip2=iqprop(iq2,PROP2_TYPE,r2);
	    
	    if(irev==1) std::swap(ip1,ip2); //select the propagator to revert
	    
	    //compute the hadronic part
	    hadronic_part_leptonic_correlation(hadr,Q[ip1],Q[ip2]);
	    
	    for(int orie=0;orie<norie;orie++)
	      for(int rl=0;rl<nr;rl++)
		{
		  //contract with lepton
		  int ilins=(qins!=0);
		  int iprop=ilprop(ilepton,ilins,orie,rl);
		  int ind=hadrolept_corrpack_ind(rl,orie,r2,irev,qins,ilepton);
		  attach_leptonic_correlation(hadr,iprop,ilepton,orie,rl,ind);
		}
	  }
  
  lept_contr_time+=take_time();
}

//print out correlations
void print_correlations()
{
  print_time-=take_time();
  
  //open file and reduce
  FILE *fout=open_file(combine("%s/corr_hl",outfolder).c_str(),"w");
  glb_nodes_reduce_complex_vect(hadrolept_corr,glbSize[0]*nweak_ind*nhadrolept_proj*nind);
  
  //write down
  for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int qins=0;qins<nins;qins++)
      for(int irev=0;irev<nrev;irev++)
	for(int r2=0;r2<nr;r2++)
	  {
	    //takes the index of the quarks
	    int iq1=lep_corr_iq1[ilepton];
	    int iq2=lep_corr_iq2[ilepton];
	    if(irev==1) std::swap(iq1,iq2);
	    for(int orie=0;orie<norie;orie++)
	      for(int rl=0;rl<nr;rl++)
		{
		  int corrpack_ind=hadrolept_corrpack_ind(rl,orie,r2,irev,qins,ilepton);
		  
		  if(!pure_wilson) master_fprintf(fout," # mq1=%lg mq2=%lg qins=%d qrev=%d rq1=%d rq2=%d lep_orie=%+d rl=%d\n\n",
						  qmass[iq1],qmass[iq2],qins,irev+1,!r2,r2,sign_orie[orie],rl);
		  else             master_fprintf(fout," # kappaq1=%lg kappaq2=%lg qins=%d qrev=%d lep_orie=%+d\n\n",
						  qkappa[iq1],qkappa[iq2],qins,irev+1,sign_orie[orie]);
		  for(int ind=0;ind<nweak_ind;ind++)
		    for(int ig_proj=0;ig_proj<nhadrolept_proj;ig_proj++)
		      {
			master_fprintf(fout," # qins=%s lins=%s proj=%s\n\n",list_weak_ind_nameq[ind],list_weak_ind_namel[ind],gtag[hadrolept_projs[ig_proj]]);
			for(int t=0;t<glbSize[0];t++)
			  {
			    int i=t+glbSize[0]*(ig_proj+nhadrolept_proj*(ind+nweak_ind*corrpack_ind));
			    master_fprintf(fout,"%+16.16lg %+16.16lg\n",hadrolept_corr[i][RE]/(glb_spat_vol*nsources),hadrolept_corr[i][IM]/(glb_spat_vol*nsources));
			  }
			master_fprintf(fout,"\n");
		      }
		}
	  }
  close_file(fout);
  
  /////////////////////////////////// purely hadronic part ////////////////////////////////////////////
  
  //normalise
  double n=1.0/nsources;
  for(int i=0;i<hadr_corr_length;i++) complex_prodassign_double(hadr_corr[i],n);
  
  int ind=0;
  for(int icombo=0;icombo<ncombo_hadr_corr;icombo++)
    {
      fout=open_file(combine("%s/corr_%c%c",outfolder,prop_abbr[prop1_hadr_map[icombo]],prop_abbr[prop2_hadr_map[icombo]]).c_str(),"w");
      
      for(int imass=0;imass<nqmass;imass++)
	for(int jmass=0;jmass<nqmass;jmass++)
	  for(int r=0;r<nr;r++)
	    {
	      if(!pure_wilson) master_fprintf(fout," # m1(rev)=%lg m2(ins)=%lg r=%d\n",qmass[imass],qmass[jmass],r);
	      else             master_fprintf(fout," # kappa1(rev)=%lg kappa2(ins)=%lg\n",qkappa[imass],qkappa[jmass]);
	      print_contractions_to_file(fout,nhadr_contr,ig_hadr_so,ig_hadr_si,hadr_corr+ind*glbSize[0],0,"",1.0/glb_spat_vol);
	      master_fprintf(fout,"\n");
	      ind+=nhadr_contr;
	    }
      
      //close the file
      close_file(fout);
    }
  
  print_time+=take_time();
}

//close deallocating everything
void close()
{
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g, of which:\n",tot_prog_time);
  master_printf(" - %02.2f%s to prepare %d photon stochastic propagators (%2.2gs avg)\n",photon_prop_time/tot_prog_time*100,"%",nphoton_prop_tot,photon_prop_time/nphoton_prop_tot);
  master_printf(" - %02.2f%s to prepare %d lepton propagators (%2.2gs avg)\n",lepton_prop_time/tot_prog_time*100,"%",nlprop,lepton_prop_time/nlprop);
  master_printf(" - %02.2f%s to prepare %d generalized sources (%2.2gs avg)\n",source_time/tot_prog_time*100,"%",nsource_tot,source_time/nsource_tot);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf("    of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d hadronic contractions (%2.2gs avg)\n",hadr_contr_time/tot_prog_time*100,"%",nhadr_contr_tot,hadr_contr_time/nhadr_contr_tot);
  master_printf(" - %02.2f%s to perform %d leptonic contractions (%2.2gs avg)\n",lept_contr_time/tot_prog_time*100,"%",nlept_contr_tot,lept_contr_time/nlept_contr_tot);
  master_printf(" - %02.2f%s to print hadro-leptonic contractions\n",print_time/tot_prog_time*100,"%");
  
  nissa_free(photon_field);
  nissa_free(source);
  nissa_free(original_source);
  for(int iprop=0;iprop<nqprop;iprop++) nissa_free(Q[iprop]);
  nissa_free(Q);
  for(int iprop=0;iprop<nlprop;iprop++) nissa_free(L[iprop]);
  nissa_free(L);
  nissa_free(temp_lep);
  nissa_free(conf);
  nissa_free(hadr_corr);
  nissa_free(glb_corr);
  nissa_free(loc_corr);
  nissa_free(hadr);
  nissa_free(hadrolept_corr);
  nissa_free(temp_source);
  nissa_free(temp_solution);
  nissa_free(lep_corr_iq1);
  nissa_free(lep_corr_iq2);
  nissa_free(leps);
  nissa_free(lep_energy);
  nissa_free(neu_energy);
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  //loop over the configs
  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf && enough_time && !file_exists("stop") && read_conf_parameters(iconf))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      for(int isource=0;isource<nsources;isource++)
	{
	  master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
	  
	  //init
	  random_shift_gauge_conf(conf,old_theta,put_theta);
	  generate_photon_stochastic_propagator();
	  generate_original_source();
	  
	  generate_lepton_propagators();
	  generate_quark_propagators();
	  
	  compute_hadroleptonic_correlations();
	  compute_hadronic_correlations();
	}
      
      //print out correlations
      print_correlations();
      
      //pass to the next conf if there is enough time
      char fin_file[1024];
      if(snprintf(fin_file,1024,"%s/finished",outfolder)<0) crash("creating path for finished file");
      file_touch(fin_file);
      
      nanalyzed_conf++;
      enough_time=check_remaining_time();
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
