#include <nissa.hpp>

#ifdef POINT_SOURCE_VERSION
 #define PROP_TYPE su3spinspin
#else
 #define PROP_TYPE colorspinspin
#endif

using namespace nissa;

/////////////////////////////////////// data //////////////////////////////

int ninv_tot=0,nhadr_contr_tot=0,nlept_contr_tot=0,nsource_tot=0,nphoton_prop_tot=0;
double inv_time=0,hadr_contr_time=0,lept_contr_time=0,print_time=0;
double tot_prog_time=0,source_time=0,photon_prop_time=0,lepton_prop_time=0;

int wall_time;
int free_theory,rnd_gauge_transform;
int ngauge_conf,nanalyzed_conf=0;
char conf_path[1024],outfolder[1024];
quad_su3 *conf;

double kappa;
double put_theta[4],old_theta[4]={0,0,0,0};

coords source_coord;
PROP_TYPE *source,*original_source;
int seed,noise_type;

int nqmass,nr,nsources;
double *qmass,*residue;
PROP_TYPE **Q;

spincolor *temp_source;
spincolor *temp_solution;
  
gauge_info photon;
double tadpole[4];
spin1field *photon_phi,*photon_eta;

int hadr_corr_length;
complex *hadr_ave;
complex *hadr_err;
complex *glb_corr,*loc_corr;

//sign of the muon momentum
const int sign_orie[2]={-1,+1};

//list the 8 matrices to insert for the weak current
const int nweak_ins=16;
const int nvitt_g_proj=2,vitt_g_projs[nvitt_g_proj]={4,9};
int list_weak_insq[nweak_ins]={1,2,3,4, 6,7,8,9,  1,2,3,4, 6,7,8,9};
int list_weak_insl[nweak_ins]={1,2,3,4, 6,7,8,9,  6,7,8,9, 1,2,3,4};
int nind;
spinspin *hadr;
complex *glb_weak_vitt_corr;
complex *glb_weak_vitt_ave;
complex *glb_weak_vitt_err;
complex *glb_weak_proj_corr;

//compute the eigenvalues of (1-+g0)/2
double W=1/sqrt(2);
spin ompg0_eig[2][2]={{{{+W, 0},{ 0, 0},{+W, 0},{ 0, 0}},
		       {{ 0, 0},{+W, 0},{ 0, 0},{+W, 0}}},
		      {{{+W, 0},{ 0, 0},{-W, 0},{ 0, 0}},
		       {{ 0, 0},{+W, 0},{ 0, 0},{-W, 0}}}};


//#define DEBUG
//#define NOINSQ
//#define ONLYLEP

//define types of quark propagator used
const int nins_kind=6;
enum insertion_t{                    ORIGINAL,  SCALAR,  PSEUDO,  STOCH_PHI,  STOCH_ETA,  TADPOLE};
const char ins_name[nins_kind][20]={"ORIGINAL","SCALAR","PSEUDO","STOCH_PHI","STOCH_ETA","TADPOLE"};
const int nqprop_kind=7;
enum qprop_t{                           PROP_0,  PROP_S,  PROP_P,  PROP_PHI,  PROP_ETA,  PROP_PHIETA,  PROP_T};
const char prop_name[nqprop_kind][20]={"PROP_0","PROP_S","PROP_P","PROP_PHI","PROP_ETA","PROP_PHIETA","PROP_T"};
#ifndef NOINSQ
const qprop_t PROP_PHI_ETA[2]={PROP_PHI,PROP_ETA};
#endif
//map the source, the destination and the insertion for each propagator
const qprop_t prop_map[nqprop_kind]=         {PROP_0,   PROP_S, PROP_P, PROP_PHI,  PROP_ETA,  PROP_PHIETA, PROP_T};
const insertion_t insertion_map[nqprop_kind]={ORIGINAL, SCALAR, PSEUDO, STOCH_PHI, STOCH_ETA, STOCH_ETA,   TADPOLE};
const qprop_t source_map[nqprop_kind]=       {PROP_0,   PROP_0, PROP_0, PROP_0,    PROP_0,    PROP_PHI,    PROP_0};
const char prop_abbr[]=                       "0"       "S"     "P"     "A"        "B"        "X"          "T";

//hadron contractions
const int ncombo_hadr_corr=9;
const qprop_t prop1_hadr_map[ncombo_hadr_corr]={PROP_0,PROP_0,PROP_0,PROP_0,     PROP_0,PROP_PHI};
const qprop_t prop2_hadr_map[ncombo_hadr_corr]={PROP_0,PROP_S,PROP_P,PROP_PHIETA,PROP_T,PROP_ETA};

//parameters of the leptons
int nleptons;
int *lep_corr_iq1;
int *lep_corr_iq2;
int lepton_mom_sign[2]={-1,+1};
tm_quark_info *leps;
double *lep_energy,*neu_energy;
spinspin **L,*temp_lep;

////////////////////////////////////////////// get propagator and prop. info /////////////////////////////////////////////

//return appropriate propagator
int nqprop,nlprop;
int iqprop(int imass,qprop_t ip,int r)
{return r+nr*(imass+nqmass*ip);}
int ilprop(int ilepton,int orie,int phi_eta,int r)
{return r+nr*(phi_eta+2*(orie+2*ilepton));}

//return appropriately modified info
tm_quark_info get_lepton_info(int ilepton,int orie,int r)
{
  tm_quark_info le=leps[ilepton];
  le.r=r;
  for(int i=1;i<4;i++) le.bc[i]*=sign_orie[orie];
  
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
  read_str_int("WallTime",&wall_time);
  //Kappa
  read_str_double("Kappa",&kappa);
  //One or two r
  read_str_int("NR",&nr);
  //Masses and residue
  read_list_of_double_pairs("QMassResidues",&nqmass,&qmass,&residue);
  
  //Leptons
  read_str_int("LeptonicCorrs",&nleptons);
  lep_corr_iq1=nissa_malloc("lep_corr_iq1",nleptons,int);
  lep_corr_iq2=nissa_malloc("lep_corr_iq2",nleptons,int);
  leps=nissa_malloc("leps",nleptons,tm_quark_info);
  lep_energy=nissa_malloc("lep_energy",nleptons,double);
  neu_energy=nissa_malloc("neu_energy",nleptons,double);
  expect_str("Q1Q2LepmassMesmass");
  for(int il=0;il<nleptons;il++)
    {
      //read quarks identfiying the mesons, and lepton mass
      read_int(lep_corr_iq1+il);
      read_int(lep_corr_iq2+il);
      read_double(&leps[il].mass);
      
      //maximal twist and antiperiodic
      leps[il].bc[0]=1;
      leps[il].kappa=0.125;
      leps[il].r=0;
      
      //read the mass of the meson (that must have been determined outside)
      double mes_mass;
      read_double(&mes_mass);
      
      double free_quark_ener[2];
      for(int iq=0;iq<2;iq++)
	{
	  tm_quark_info q;
	  q.kappa=0.125;
	  q.bc[0]=1;
	  for(int i=1;i<4;i++) q.bc[i]=0;
	  q.mass=qmass[((iq==0)?lep_corr_iq1:lep_corr_iq2)[il]];
	  free_quark_ener[iq]=tm_quark_energy(q,0);
	  master_printf(" supposed free quark energy[%d]: %+016.016lg\n",iq,free_quark_ener[iq]);
	}
      double free_mes_ener=free_quark_ener[0]+free_quark_ener[1];
      master_printf(" supposed free meson energy: %+016.016lg\n",free_mes_ener);
      
      //compute meson momentum and bc
      for(int i=1;i<4;i++) leps[il].bc[i]=(sqr(mes_mass)-sqr(leps[il].mass))/(2*mes_mass)/sqrt(3)/M_PI*glb_size[i];
      double err;
      do
      	{
      	  //compute the error
	  double lep_energy=tm_quark_energy(leps[il],0);
	  double neu_energy=naive_massless_quark_energy(leps[il].bc,0);
      	  err=lep_energy+neu_energy-mes_mass;
      	  //compute the derivative
      	  double eps=1e-8;
      	  for(int i=1;i<4;i++) leps[il].bc[i]+=eps;
      	  double der=(tm_quark_energy(leps[il],0)+naive_massless_quark_energy(leps[il].bc,0)-mes_mass-err)/eps;
      	  for(int i=1;i<4;i++) leps[il].bc[i]-=eps+err/der;
	  
      	  master_printf("rank %d lep_e: %+010.10lg, neu_e: %+010.10lg, mes_mass: %lg, error: %lg, der: %lg\n",rank,lep_energy,neu_energy,mes_mass,err,der);
      	}
      while(fabs(err)>1e-14);
      
      //write down energy
      lep_energy[il]=tm_quark_energy(leps[il],0);
      neu_energy[il]=naive_massless_quark_energy(leps[il].bc,0);
      master_printf(" ilepton %d, lepton mass %lg, lepton energy: %lg, neutrino energy: %lg\n",il,leps[il].mass,lep_energy[il],neu_energy[il]);
      master_printf(" lep+neut energy: %lg\n",lep_energy[il]+neu_energy[il]);
    }
  
  //Zero mode subtraction
  char zero_mode_sub_str[100];
  read_str_str("ZeroModeSubtraction",zero_mode_sub_str,100);
  
  if(strncasecmp(zero_mode_sub_str,"PECIONA",100)==0) photon.zms=PECIONA;
  else
    if(strncasecmp(zero_mode_sub_str,"UNNO_ALEMANNA",100)==0) photon.zms=UNNO_ALEMANNA;
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
  
  //perform a random gauge transformation
  read_str_int("RandomGaugeTransform",&rnd_gauge_transform);
  
  //number of sources
  read_str_int("NSources",&nsources);
  
  //number of configurations
  read_str_int("NGaugeConf",&ngauge_conf);
  
  ///////////////////// finihed reading apart from conf list ///////////////
  
  //compute the tadpole summing all momentum
  compute_tadpole(tadpole,photon);
  
  //Allocate
  nqprop=iqprop(nqmass-1,prop_map[nqprop_kind-1],nr-1)+1;
  nlprop=ilprop(nleptons-1,1,1,nr-1)+1;
  
  //allocate temporary vectors
  temp_source=nissa_malloc("temp_source",loc_vol,spincolor);
  temp_solution=nissa_malloc("temp_solution",loc_vol,spincolor);
  hadr_corr_length=glb_size[0]*ncombo_hadr_corr*nqmass*nqmass*nr;
  hadr_ave=nissa_malloc("hadr_ave",hadr_corr_length,complex);
  hadr_err=nissa_malloc("hadr_err",hadr_corr_length,complex);
  glb_corr=nissa_malloc("glb_corr",glb_size[0],complex);
  loc_corr=nissa_malloc("loc_corr",glb_size[0],complex);
  nind=nleptons*nweak_ins*2*2*nr;
  hadr=nissa_malloc("hadr",loc_vol,spinspin);
  glb_weak_proj_corr=nissa_malloc("glb_weak_proj_corr",glb_size[0]*nweak_ins*4*nind,complex);
  glb_weak_vitt_corr=nissa_malloc("glb_weak_vitt_corr",glb_size[0]*nweak_ins*nvitt_g_proj*nind,complex);
  glb_weak_vitt_ave=nissa_malloc("glb_weak_vitt_ave",glb_size[0]*nweak_ins*nvitt_g_proj*nind,complex); 
  glb_weak_vitt_err=nissa_malloc("glb_weak_vitt_err",glb_size[0]*nweak_ins*nvitt_g_proj*nind,complex); 
  original_source=nissa_malloc("source",loc_vol,PROP_TYPE);
  source=nissa_malloc("source",loc_vol,PROP_TYPE);
  photon_eta=nissa_malloc("photon_eta",loc_vol+bord_vol,spin1field);
  photon_phi=nissa_malloc("photon_phi",loc_vol+bord_vol,spin1field);
  Q=nissa_malloc("Q*",nqprop,PROP_TYPE*);
  for(int iprop=0;iprop<nqprop;iprop++) Q[iprop]=nissa_malloc("Q",loc_vol+bord_vol,PROP_TYPE);
  L=nissa_malloc("L*",nlprop,spinspin*);
  for(int iprop=0;iprop<nlprop;iprop++) L[iprop]=nissa_malloc("L",loc_vol+bord_vol,spinspin);
  temp_lep=nissa_malloc("temp_lep",loc_vol+bord_vol,spinspin);
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
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
      sprintf(fin_file,"%s/finished",outfolder);
      sprintf(run_file,"%s/running",outfolder);
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
	master_printf(" In output path \"%s\" terminating file already present: configuration \"%s\" already analyzed, skipping.\n",outfolder,conf_path);
      iconf++;
    }
  while(!ok_conf && iconf<ngauge_conf);
  
  master_printf("\n");
  
  return ok_conf;
}

//read the conf and setup it
void setup_conf()
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
  
  //reset correlations
  vector_reset(hadr_ave);
  vector_reset(hadr_err);
  vector_reset(glb_weak_vitt_ave);
  vector_reset(glb_weak_vitt_err);
  vector_reset(glb_weak_proj_corr);
}

//generate a wall-source for stochastic QCD propagator
void generate_original_source()
{
  //source coord
  //coords M={glb_size[0]/2,glb_size[1],glb_size[2],glb_size[3]}; //temporarily commented to fix the source in 0
  for(int mu=0;mu<4;mu++) source_coord[mu]=0;//(int)(rnd_get_unif(&glb_rnd_gen,0,1)*M[mu]);
  
#ifdef POINT_SOURCE_VERSION
  master_printf("Source position: t=%d x=%d y=%d z=%d\n",source_coord[0],source_coord[1],source_coord[2],source_coord[3]);
  generate_delta_source(original_source,source_coord);
#else
  master_printf("Source position: t=%d\n",source_coord[0]);
  generate_spindiluted_source(original_source,rnd_type_map[noise_type],source_coord[0]);
#endif
}

//////////////////////////////////////// quark propagators /////////////////////////////////////////////////

//generate a sequential source
void generate_source(insertion_t inser,int r,PROP_TYPE *ori)
{
  source_time-=take_time();

#ifdef DEBUG
  int dirs[4]={0,1,0,0};
#endif
  switch(inser)
    {
    case ORIGINAL:prop_multiply_with_gamma(source,0,original_source);break;
    case SCALAR:prop_multiply_with_gamma(source,0,ori);break;
    case PSEUDO:prop_multiply_with_gamma(source,5,ori);break;
#ifdef DEBUG
    case STOCH_PHI:insert_conserved_current(source,conf,ori,r,dirs);break;
    case STOCH_ETA:insert_conserved_current(source,conf,ori,r,dirs);break;
#else
    case STOCH_PHI:insert_external_source(source,conf,photon_phi,ori,r);break;
    case STOCH_ETA:insert_external_source(source,conf,photon_eta,ori,r);break;
#endif
    case TADPOLE:insert_tadpole(source,conf,ori,r,tadpole);break;
    }
  
  source_time+=take_time();
  nsource_tot++;
}

//invert on top of a source, putting all needed for the appropriate quark
void get_qprop(PROP_TYPE *out,PROP_TYPE *in,int imass,bool r,int rotate=true)
{
  //these are the way in which Dirac operator rotate - propagator is opposite, see below
  
#ifdef POINT_SOURCE_VERSION
  for(int ic=0;ic<3;ic++)
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
	if(rotate) safe_dirac_prod_spincolor(temp_source,(tau3[r]==-1)?&Pminus:&Pplus,temp_source);
	
	//invert
	inv_time-=take_time();
	inv_tmD_cg_eoprec_eos(temp_solution,NULL,conf,kappa,tau3[r]*qmass[imass],100000,residue[imass],temp_source);
	ninv_tot++;inv_time+=take_time();
	
	//rotate the sink index
	if(rotate) safe_dirac_prod_spincolor(temp_solution,(tau3[r]==-1)?&Pminus:&Pplus,temp_solution);      
	
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
	    master_printf(" mass[%d]=%lg, r=%d\n",imass,qmass[imass],r);
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
  generate_stochastic_tlSym_gauge_propagator(photon_phi,photon_eta,photon);
  photon_prop_time+=take_time();
  nphoton_prop_tot++;
}

/////////////////////////////////////////////// lepton propagators ///////////////////////////////////////////

//compute phase exponent for space part: vec{p}*\vec{x}
double get_space_arg(int ivol,momentum_t bc)
{
  double arg=0;
  for(int mu=1;mu<4;mu++)
    {
      double step=bc[mu]*M_PI/glb_size[mu];
      arg+=step*glb_coord_of_loclx[ivol][mu];
    }
  return arg;
}

//compute the phase for lepton on its sink
void get_lepton_sink_phase_factor(complex out,int ivol,int ilepton,tm_quark_info le)
{
  //compute space and time factor
  double arg=get_space_arg(ivol,le.bc);
  int t=(glb_coord_of_loclx[ivol][0]-source_coord[0]+glb_size[0])%glb_size[0];
  double ext=exp(t*lep_energy[ilepton]);
  
  //compute full exponential (notice the factor -1)
  out[RE]=cos(-arg)*ext;
  out[IM]=sin(-arg)*ext;
}

//compute the phase for antineutrino - the orientation is that of the muon (as above)
void get_antineutrino_source_phase_factor(complex out,int ivol,int ilepton,momentum_t bc)
{
  //compute space and time factor
  double arg=get_space_arg(ivol,bc);
  int t=(glb_coord_of_loclx[ivol][0]-source_coord[0]+glb_size[0])%glb_size[0];
  //t*=(1-2*(t>=glb_size[0]/2));
  if(t>=glb_size[0]/2) t=glb_size[0]-t;
  double ext=exp(t*neu_energy[ilepton]);
  
  //compute full exponential (notice the factor +1)
  out[RE]=cos(+arg)*ext;
  out[IM]=sin(+arg)*ext;
}

//set everything to a phase factor
void set_to_lepton_sink_phase_factor(spinspin *prop,int ilepton,tm_quark_info &le,int twall)
{
  GET_THREAD_ID();
  
  vector_reset(prop);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    if(twall==-1||glb_coord_of_loclx[ivol][0]==twall)
      {
	complex ph;
	get_lepton_sink_phase_factor(ph,ivol,ilepton,le);
	spinspin_put_to_diag(prop[ivol],ph);
      }
  set_borders_invalid(prop);
}

void trace_test_lep_prop_source(complex *c,spinspin *prop,tm_quark_info le,int ilepton)
{
  GET_THREAD_ID();

  if(IS_MASTER_THREAD)
    {
      //trace
      int imom=0;
      int par_apar=0;
      memset(c,0,sizeof(complex)*glb_size[0]);

      //fetch the energy and momentum ot the muon
      //double E=lep_energy[ilepton];

      //compute steps
      momentum_t steps;
      for(int mu=1;mu<4;mu++)
	steps[mu]=+le.bc[mu]*M_PI/glb_size[mu];
      
      spinspin pr;
      twisted_particle_anti_particle_projector_of_imom(pr,le,imom,par_apar);

      NISSA_LOC_VOL_LOOP(ivol)
	{
	  //compute phase exponent
	  double arg=0;
	  for(int mu=1;mu<4;mu++)
	    arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	  
	  //compute pure time factor
	  double ext=1; //ANNA
	  //double ext=exp(-(glb_size[0]-glb_coord_of_loclx[ivol][0])*E);
	  
	  //compute the phase
	  complex ph={ext*cos(arg),ext*sin(arg)};
	  
	  spinspin o;
	  unsafe_spinspin_prod_spinspin(o,pr,prop[ivol]);
	  complex tr;
	  trace_spinspin(tr,o);
	  complex_summ_the_prod(c[glb_coord_of_loclx[ivol][0]],tr,ph);
	}
      MPI_Allreduce(MPI_IN_PLACE,c,2*glb_size[0],MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
  THREAD_BARRIER();
}

//insert the photon on the source side
void insert_photon_on_the_source(spinspin *prop,int ilepton,int phi_eta,coords dirs,tm_quark_info &le)
{ 
  GET_THREAD_ID();

  //select A
  spin1field *A=(phi_eta==0)?photon_phi:photon_eta;
  communicate_lx_spin1field_borders(A);
  
  //phases
  complex phases[4];
  for(int mu=0;mu<NDIM;mu++)
    {
      phases[mu][0]=cos(le.bc[mu]*M_PI);
      phases[mu][1]=sin(le.bc[mu]*M_PI);
    }
  
  //copy on the temporary and communicate borders
  vector_copy(temp_lep,prop);
  communicate_lx_spinspin_borders(temp_lep);
  vector_reset(prop);
  
  //prepare each propagator for a single lepton
  //by computing i(phi(x-mu)A_mu(x-mu)(-i t3 g5-gmu)/2-phi(x+mu)A_mu(x)(-i t3 g5+gmu)/2)=
  //(ph0 A_mu(x-mu)g[r][0][mu]-ph0 A_mu(x)g[r][1][mu])=
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      {
	//find neighbors
	int ifw=loclx_neighup[ivol][mu];
	int ibw=loclx_neighdw[ivol][mu];
	
	//compute phase factor
	spinspin ph_bw,ph_fw;

	//transport down and up
	if(glb_coord_of_loclx[ivol][mu]==glb_size[mu]-1) unsafe_spinspin_prod_complex_conj2(ph_fw,temp_lep[ifw],phases[mu]);
	else spinspin_copy(ph_fw,temp_lep[ifw]);
	if(glb_coord_of_loclx[ivol][mu]==0) unsafe_spinspin_prod_complex(ph_bw,temp_lep[ibw],phases[mu]);
	else spinspin_copy(ph_bw,temp_lep[ibw]);
	
	//fix coefficients - i is inserted here!
	//also dir selection is made here
#ifndef DEBUG
	spinspin_prodassign_idouble(ph_fw,-0.5*dirs[mu]);
	spinspin_prodassign_idouble(ph_bw,+0.5*dirs[mu]);
	
	//fix insertion of the current
	safe_spinspin_prod_complex(ph_fw,ph_fw,A[ivol][mu]);
	safe_spinspin_prod_complex(ph_bw,ph_bw,A[ibw][mu]);
#else
	spinspin_prodassign_double(ph_fw,-0.5*dirs[mu]);
	spinspin_prodassign_double(ph_bw,+0.5*dirs[mu]);
#endif
	
	//summ and subtract the two
	spinspin bw_M_fw,bw_P_fw;
	spinspin_subt(bw_M_fw,ph_bw,ph_fw);
	spinspin_summ(bw_P_fw,ph_bw,ph_fw);
		  
	//put -i g5 t3 on the summ
	spinspin temp_P;
	spinspin_prodassign_idouble(bw_P_fw,-tau3[le.r]);
	unsafe_spinspin_prod_dirac(temp_P,bw_P_fw,base_gamma+5);
	spinspin_summassign(prop[ivol],temp_P);
	
	//put gmu on the diff
	spinspin temp_M;
	unsafe_spinspin_prod_dirac(temp_M,bw_M_fw,base_gamma+map_mu[mu]);
	spinspin_summassign(prop[ivol],temp_M);
      }
  set_borders_invalid(prop);
}
//insert the photon on the source
void insert_photon_on_the_source(spinspin *prop,int ilepton,int phi_eta,tm_quark_info &le)
{
  coords dirs={1,1,1,1};
  insert_photon_on_the_source(prop,ilepton,phi_eta,dirs,le);
}

//generate all the lepton propagators, pointing outward
//the computations is done by:
// 1)putting the correct phase in x space, given by exp(E_mu*t-i*vec(p)*vec(x))
// 2)multiplying it by the conserved current inserting eta or phi
// 3)going to momentum space
// 4)multiplying by the lepton propagator in momentum space
// 5)coming back to x space
THREADABLE_FUNCTION_0ARG(generate_lepton_propagators)
{
  GET_THREAD_ID();

  if(IS_MASTER_THREAD) lepton_prop_time-=take_time();
  master_printf("Generating lepton propagators\n");

  for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int orie=0;orie<2;orie++)
      for(int phi_eta=0;phi_eta<2;phi_eta++)
	for(int r=0;r<nr;r++)
	  {
	    //set the properties of the meson
	    //time boundaries are anti-periodic, space are as for external line
	    tm_quark_info le=get_lepton_info(ilepton,orie,r);
	    
	    //select the propagator
	    int iprop=ilprop(ilepton,orie,phi_eta,r);
	    spinspin *prop=L[iprop];
	    	    
	    //put it to a phase
	    int twall=((glb_size[0]/2+0+source_coord[0])%glb_size[0]);
	    set_to_lepton_sink_phase_factor(prop,ilepton,le,twall);
	    
	    //multiply and the insert the current in between, on the source side
	    multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le);
#ifndef DEBUG
	    insert_photon_on_the_source(prop,ilepton,phi_eta,le);
#else
	    int dir0[4]={1,0,0,0};
	    insert_photon_on_the_source(prop,ilepton,phi_eta,dir0,le);
#endif
	    multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le);
	  }
  
  if(IS_MASTER_THREAD) lepton_prop_time+=take_time();
}
THREADABLE_FUNCTION_END

////////////////////////////////////////// purely hadronic correlators ///////////////////////////////////////////

//compute all the hadronic correlations
void compute_hadronic_correlation()
{
  master_printf("Computing hadronic correlation functions\n");
  
  hadr_contr_time-=take_time();
  int ind=0;
  for(int icombo=0;icombo<ncombo_hadr_corr;icombo++)
    for(int imass=0;imass<nqmass;imass++)
      for(int jmass=0;jmass<nqmass;jmass++)
	for(int r=0;r<nr;r++)
	  {
	    //compute the correlation function
	    int ig_so=5,ig_si=5;
	    meson_two_points_Wilson_prop(glb_corr,loc_corr,&ig_so,Q[iqprop(imass,prop1_hadr_map[icombo],r)],&ig_si,Q[iqprop(jmass,prop2_hadr_map[icombo],r)],1);
	    
	    //save to the total stack
	    for(int t=0;t<glb_size[0];t++)
	      {
		int i=(t-source_coord[0]+glb_size[0])%glb_size[0]+glb_size[0]*ind;
		complex_summassign(hadr_ave[i],glb_corr[t]);
		complex temp={sqr(glb_corr[t][RE]),sqr(glb_corr[t][IM])};
		complex_summassign(hadr_err[i],temp);
	      }
	    ind++;
	  }
  hadr_contr_time+=take_time();
  nhadr_contr_tot+=ind;
}

//compute all the hadronic correlations
void print_hadronic_correlations()
{
  //normalise
  double n=1.0/nsources,n2=1.0/std::max(1,nsources-1);
  for(int i=0;i<hadr_corr_length;i++)
    {
      //normalise
      complex_prodassign_double(hadr_ave[i],n);
      complex_prodassign_double(hadr_err[i],n);
      
      //compute the true error
      complex temp={sqr(hadr_ave[i][RE]),sqr(hadr_ave[i][IM])};
      complex_subtassign(hadr_err[i],temp);
      for(int ri=0;ri<2;ri++) hadr_err[i][ri]=sqrt(fabs(hadr_err[i][ri])*n2);
    }
  
  int ind=0;
  for(int icombo=0;icombo<ncombo_hadr_corr;icombo++)
    {
      FILE *fout[2];
      const char suff[2][10]={"corr","err"};
      for(int i=0;i<2;i++) fout[i]=open_file(combine("%s/%s_%c%c",outfolder,suff[i],prop_abbr[prop1_hadr_map[icombo]],prop_abbr[prop2_hadr_map[icombo]]).c_str(),"w");
      
      for(int imass=0;imass<nqmass;imass++)
	for(int jmass=0;jmass<nqmass;jmass++)
	  for(int r=0;r<nr;r++)
	    {
	      //print out
	      int ig_so=5,ig_si=5;

	      complex *hadr[2]={hadr_ave,hadr_err};
	      for(int i=0;i<2;i++)
		{
		  master_fprintf(fout[i]," # m1(rev)=%lg m2(ins)=%lg r=%d\n",qmass[imass],qmass[jmass],r);
		  print_contractions_to_file(fout[i],1,&ig_so,&ig_si,hadr[i]+ind*glb_size[0],0,"",1.0);
		  master_fprintf(fout[i],"\n");
		}
	      ind++;
	    }
      
      //close the file
      for(int i=0;i<2;i++) close_file(fout[i]);
    }
}

/////////////////////////////////////////////// hadroleptonic correlators //////////////////////////////////////////

//compute the hadronic part of the lepton correlation function
//as usual, FIRST propagator is reverted
THREADABLE_FUNCTION_3ARG(hadronic_part_leptonic_correlation, spinspin*,hadr, PROP_TYPE*,S1, PROP_TYPE*,S2)
{
  GET_THREAD_ID();
  
  vector_reset(hadr);
  
  //it's just the matter of inserting gamma5*gamma5=identity between S1^dag and S2
  //on the sink gamma5 must still be inserted!
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int ic_si=0;ic_si<3;ic_si++)
#ifdef POINT_SOURCE_VERSION
      for(int ic_so=0;ic_so<3;ic_so++)
#endif
	for(int id_si1=0;id_si1<4;id_si1++)
	  for(int id_si2=0;id_si2<4;id_si2++)
	    for(int id_so=0;id_so<4;id_so++)
	      complex_summ_the_conj1_prod
		(hadr[ivol][id_si2][id_si1], //this way when taking the trace with dirac matrix, that is acting on S2, as it should
#ifdef POINT_SOURCE_VERSION
		 S1[ivol][ic_si][ic_so][id_si1][id_so],S2[ivol][ic_so][ic_si][id_si2][id_so])
#else
                 S1[ivol][ic_si][id_si1][id_so],S2[ivol][ic_si][id_si2][id_so])
#endif
		 ;
  THREAD_BARRIER();
}
THREADABLE_FUNCTION_END

//compute the polarisation vector
void get_polvect(spin *u,spin *vbar,tm_quark_info &le)
{
  for(int s=0;s<2;s++)
    {
      twisted_wavefunction_of_imom(u[s],le,0,0,s);
      spin v;
      naive_massless_wavefunction_of_imom(v,le.bc,0,1,s);
      spin vdag;
      spin_conj(vdag,v);
      unsafe_spin_prod_dirac(vbar[s],vdag,base_gamma+map_mu[0]);
    }
}

//compute the polarisation vector barred
void get_polvect_bar(spin *ubar,spin *v,tm_quark_info &le)
{
  for(int s=0;s<2;s++)
    {
      naive_massless_wavefunction_of_imom(v[s],le.bc,0,1,s);
      spin u;
      twisted_wavefunction_of_imom(u,le,0,0,s);
      spin udag;
      spin_conj(udag,u);
      unsafe_spin_prod_dirac(ubar[s],udag,base_gamma+map_mu[0]);
    }
}

//compute the leptonic part of the correlation function
THREADABLE_FUNCTION_6ARG(attach_leptonic_correlation, spinspin*,hadr, int,iprop, int,ilepton, int,orie, int,rl, int,ind)
{
  GET_THREAD_ID();
  
  vector_reset(loc_corr);
  
  //get the lepton info and prop
  tm_quark_info le=get_lepton_info(ilepton,orie,rl);
  spinspin *lept=L[iprop];
  
  //gets the polarisation vector
  spin u[2],vbar[2];
  get_polvect(u,vbar,le);
  spin v[2],ubar[2];
  get_polvect_bar(ubar,v,le);
  
  //gets the projectors for Vittorio
  spinspin promu[2],pronu[2];
  twisted_on_shell_operator_of_imom(promu[0],le,0,false,-1);
  twisted_on_shell_operator_of_imom(promu[1],le,0,false,+1);
  naive_massless_on_shell_operator_of_imom(pronu[0],le.bc,0,-1);
  naive_massless_on_shell_operator_of_imom(pronu[1],le.bc,0,+1);
  
  //compute the right part of the leptonic loop: G0 G^dag
  dirac_matr vitt_proj_gamma[nvitt_g_proj];
  for(int ig_proj=0;ig_proj<nvitt_g_proj;ig_proj++)
    {
      int ig=vitt_g_projs[ig_proj];
      dirac_matr temp_gamma;
      dirac_herm(&temp_gamma,base_gamma+ig);
      dirac_prod(vitt_proj_gamma+ig_proj,base_gamma+map_mu[0],&temp_gamma);
    }
  //insert gamma5 on the sink-hadron-gamma: S1^dag G5 GW S2 (G5 G5) - will dag and commutator with g0 come into?
  dirac_matr weak_ins_hadr_gamma[nweak_ins];
  for(int ins=0;ins<nweak_ins;ins++) dirac_prod(weak_ins_hadr_gamma+ins,base_gamma+5,base_gamma+list_weak_insl[ins]);
  
  for(int ins=0;ins<nweak_ins;ins++)
    {
      //define a local storage
      spinspin hl_loc_corr[loc_size[0]];
      for(int i=0;i<loc_size[0];i++) spinspin_put_to_zero(hl_loc_corr[i]);
      
      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	{
	  int t=loc_coord_of_loclx[ivol][0];
	  
	  //multiply lepton side on the right (source) side
	  spinspin l;
	  unsafe_spinspin_prod_dirac(l,lept[ivol],weak_ins_hadr_gamma+ins);
	  
	  //trace hadron side
	  complex h;
	  trace_spinspin_with_dirac(h,hadr[ivol],base_gamma+list_weak_insq[ins]);
	  
	  //get the neutrino phase (multiply hadron side) - notice that the sign of momentum is internally reversed
	  complex ph;
	  get_antineutrino_source_phase_factor(ph,ivol,ilepton,le.bc);
	  complex_prodassign(h,ph);
#ifndef ONLYLEP
	  spinspin_summ_the_complex_prod(hl_loc_corr[t],l,h);
#else
	  spinspin_summ_the_complex_prod(hl_loc_corr[t],lept[ivol],ph);
#endif
	}
      glb_threads_reduce_double_vect((double*)hl_loc_corr,loc_size[0]*sizeof(spinspin)/sizeof(double));
      
      //change sign when crossing 0 for averaging corr function properly
      for(int loc_t=0;loc_t<loc_size[0];loc_t++)
	{
	  int glb_t=loc_t+glb_coord_of_loclx[0][0];
	  int sign=1-2*(glb_t<source_coord[0]);
	  spinspin_prodassign_double(hl_loc_corr[loc_t],sign);
	}
      
      //save the decomposition on spinor
      for(int smu=0;smu<2;smu++)
	for(int snu=0;snu<2;snu++)
	  NISSA_PARALLEL_LOOP(loc_t,0,loc_size[0])
	    {
	      int glb_t=(loc_t+glb_coord_of_loclx[0][0]-source_coord[0]+glb_size[0])%glb_size[0];
	      complex hl={0,0};
	      //u multiplies the sink, vbar the source
	      spin Sv;
	      unsafe_spinspin_prod_spin(Sv,hl_loc_corr[loc_t],v[snu]);
	      for(int id_si=0;id_si<4;id_si++) complex_summ_the_prod(hl,ubar[smu][id_si],Sv[id_si]);
	      complex_summassign(glb_weak_proj_corr[glb_t+glb_size[0]*(smu+2*(snu+2*(ins+nweak_ins*ind)))],hl);
	    }
      if(IS_MASTER_THREAD) nlept_contr_tot+=4;
      
      //save decomposition Vittorio style
      for(int ig_proj=0;ig_proj<nvitt_g_proj;ig_proj++)
	NISSA_PARALLEL_LOOP(loc_t,0,loc_size[0])
	  {
	    int glb_t=(loc_t+glb_coord_of_loclx[0][0]-source_coord[0]+glb_size[0])%glb_size[0];
	    int ilnp=(glb_t>=glb_size[0]/2); //select the lepton/neutrino projector

	    spinspin td;
	    unsafe_spinspin_prod_spinspin(td,hl_loc_corr[loc_t],pronu[ilnp]);
	    spinspin dtd;
	    unsafe_spinspin_prod_spinspin(dtd,promu[ilnp],td);
	    complex hl;
	    trace_spinspin_with_dirac(hl,dtd,vitt_proj_gamma+ig_proj);
	    
	    //summ the average
	    int i=glb_t+glb_size[0]*(ig_proj+nvitt_g_proj*(ins+nweak_ins*ind));
	    complex_summassign(glb_weak_vitt_corr[i],hl);
	  }
      if(IS_MASTER_THREAD) nlept_contr_tot+=nvitt_g_proj;
      THREAD_BARRIER();
    }
  complex c[glb_size[0]];
  trace_test_lep_prop_source(c,lept,le,ilepton);
  for(int t=0;t<glb_size[0];t++)
    master_printf("%d %lg %lg\n",t,c[t][0],c[t][1]);
  //crash("");
}
THREADABLE_FUNCTION_END

//compute the total hadroleptonic correlation functions
void compute_hadroleptonic_correlation()
{
  master_printf("Computing leptonic correlation functions\n");
  lept_contr_time-=take_time();
  
  vector_reset(glb_weak_vitt_corr);
  
  int ind=0;
  for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int qins=0;qins<2;qins++)
      for(int irev=0;irev<2;irev++)
	for(int phi_eta=0;phi_eta<2;phi_eta++)
	  for(int r2=0;r2<nr;r2++)
	    {
	      //takes the index of the quarks
	      int iq1=lep_corr_iq1[ilepton];
	      int iq2=lep_corr_iq2[ilepton];
	      if(irev==1) std::swap(iq1,iq2); //select the propagator to revert
	      
	      //takes the propagators
	      qprop_t PROP1_TYPE,PROP2_TYPE;
#ifndef NOINSQ
	      if(qins==0)
		{
		  PROP1_TYPE=PROP_PHI_ETA[phi_eta];
#endif
		  PROP2_TYPE=PROP_0;
#ifndef NOINSQ
		}
	      else
		{
#endif
		  PROP1_TYPE=PROP_0;
#ifndef NOINSQ
		  PROP2_TYPE=PROP_PHI_ETA[phi_eta];
		}
#endif
	      int ip1=iqprop(iq1,PROP1_TYPE,r2); //q1 will be reverted
	      int ip2=iqprop(iq2,PROP2_TYPE,r2);
	      //compute the hadronic part
	      hadronic_part_leptonic_correlation(hadr,Q[ip1],Q[ip2]);
	      
	      for(int orie=0;orie<2;orie++)
		for(int rl=0;rl<nr;rl++)
		  {
		    //contract with lepton
		    int iprop=ilprop(ilepton,orie,!phi_eta,rl); //notice inversion of phi/eta w.r.t hadron side
		    attach_leptonic_correlation(hadr,iprop,ilepton,orie,rl,ind);
		    ind++;
		  }
	    }

  //summ to the stack
  glb_nodes_reduce_complex_vect(glb_weak_vitt_corr,glb_size[0]*nweak_ins*nvitt_g_proj*nind);
  for(int iel=0;iel<glb_size[0]*nweak_ins*nvitt_g_proj*nind;iel++)
    {
      complex_summassign(glb_weak_vitt_ave[iel],glb_weak_vitt_corr[iel]);
      complex temp={sqr(glb_weak_vitt_corr[iel][RE]),sqr(glb_weak_vitt_corr[iel][IM])};
      complex_summassign(glb_weak_vitt_err[iel],temp);
    }
  lept_contr_time+=take_time();
}

//same but projected to spin
void print_weak_proj_correlations(FILE *fout,int ind)
{
  for(int ins=0;ins<nweak_ins;ins++)
    for(int smu=0;smu<2;smu++)
      for(int snu=0;snu<2;snu++)
	{
	  master_fprintf(fout," # qins=%d lins=%d smu=%d snu=%d\n\n",list_weak_insq[ins],list_weak_insl[ins],smu,snu);
	  for(int t=0;t<glb_size[0];t++)
	    {
	      int i=t+glb_size[0]*(smu+2*(snu+2*(ins+nweak_ins*ind)));
	      master_fprintf(fout,"%+016.16lg %+016.16lg\n",glb_weak_proj_corr[i][RE]/nsources,glb_weak_proj_corr[i][IM]/nsources);
	    }
	  master_fprintf(fout,"\n");
	}
}

//same but projected vittorio way
void print_weak_vitt_correlations(FILE *fout_corr,FILE *fout_err,int ind)
{
  for(int ins=0;ins<nweak_ins;ins++)
    for(int ig_proj=0;ig_proj<nvitt_g_proj;ig_proj++)
      {
	FILE *fout[2]={fout_corr,fout_err};
	for(int i=0;i<2;i++) master_fprintf(fout[i]," # qins=%d lins=%d proj=%s\n\n",list_weak_insq[ins],list_weak_insl[ins],gtag[vitt_g_projs[ig_proj]]);
	for(int t=0;t<glb_size[0];t++)
	  {
	    int i=t+glb_size[0]*(ig_proj+nvitt_g_proj*(ins+nweak_ins*ind));
	    double A=glb_weak_vitt_ave[i][RE]/nsources;
	    double A2=glb_weak_vitt_err[i][RE]/nsources;
	    double B=glb_weak_vitt_ave[i][IM]/nsources;
	    double B2=glb_weak_vitt_err[i][IM]/nsources;
	    A2-=A*A;
	    B2-=B*B;
	    if(nsources>1)
	      {
		A2/=nsources-1;
		B2/=nsources-1;
	      }
	    A2=sqrt(fabs(A2));
	    B2=sqrt(fabs(B2));
	    master_fprintf(fout_corr,"%+016.16lg %+016.16lg\n",A,B);
	    master_fprintf(fout_err,"%+016.16lg %+016.16lg\n",A2,B2);
	  }
	for(int i=0;i<2;i++) master_fprintf(fout[i],"\n");
      }
}

//print out everything
void print_outputs()
{
  print_time-=take_time();
  
  const int nwrite=3;
  FILE *fout[3]={
    open_file(combine("%s/corr_hl_vitt",outfolder).c_str(),"w"),
    open_file(combine("%s/err_hl_vitt",outfolder).c_str(),"w"),
    open_file(combine("%s/corr_hl_proj",outfolder).c_str(),"w")};

  //write down
  int ind=0;
  for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int qins=0;qins<2;qins++)
      for(int irev=0;irev<2;irev++)
	for(int phi_eta=0;phi_eta<2;phi_eta++)
	  for(int r2=0;r2<nr;r2++)
	    {
	      //takes the index of the quarks
	      int iq1=lep_corr_iq1[ilepton];
	      int iq2=lep_corr_iq2[ilepton];
	      if(irev==1) std::swap(iq1,iq2);
	      for(int orie=0;orie<2;orie++)
		for(int rl=0;rl<nr;rl++)
		  {
		    for(int i=0;i<nwrite;i++)
		      master_fprintf(fout[i]," # mq1=%lg mq2=%lg qins=%d qrev=%d ins=%s rq1=%d rq2=%d lep_orie=%+d rl=%d\n\n",
				     qmass[iq1],qmass[iq2],qins+1,irev+1,(phi_eta==0)?"phi":"eta",!r2,r2,lepton_mom_sign[orie],rl);
		    print_weak_vitt_correlations(fout[0],fout[1],ind);
		    print_weak_proj_correlations(fout[2],ind);
		    ind++;
		  }
	    }
  for(int i=0;i<nwrite;i++) close_file(fout[i]);

  print_hadronic_correlations();
  print_time+=take_time();
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

  nissa_free(photon_eta);
  nissa_free(photon_phi);
  nissa_free(source);
  nissa_free(original_source);
  for(int iprop=0;iprop<nqprop;iprop++) nissa_free(Q[iprop]);
  nissa_free(Q);
  for(int iprop=0;iprop<nlprop;iprop++) nissa_free(L[iprop]);
  nissa_free(L);
  nissa_free(temp_lep);
  nissa_free(conf);
  nissa_free(hadr_ave);
  nissa_free(hadr_err);
  nissa_free(glb_corr);
  nissa_free(loc_corr);
  nissa_free(hadr);
  nissa_free(glb_weak_proj_corr);
  nissa_free(glb_weak_vitt_corr);
  nissa_free(glb_weak_vitt_ave);
  nissa_free(glb_weak_vitt_err);
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
      setup_conf();

      for(int isource=0;isource<nsources;isource++)
	{
	  generate_photon_stochastic_propagator();
	  generate_original_source();
	  
	  generate_lepton_propagators();
	  generate_quark_propagators();
	  
	  compute_hadroleptonic_correlation();
	  compute_hadronic_correlation();
	}
      
      //print out all the correlators
      print_outputs();
      
      //pass to the next conf if there is enough time
      char fin_file[1024];
      sprintf(fin_file,"%s/finished",outfolder);
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
