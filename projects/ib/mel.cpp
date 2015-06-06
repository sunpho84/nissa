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
complex *le_prop_ave;
complex *le_prop_err;
complex *glb_corr,*loc_corr;

//sign of the muon momentum
const int sign_orie[2]={-1,+1};

//list the 8 matrices to insert for the weak current
const int nweak_ins=16;
const int nvitt_g_proj=2,vitt_g_projs[nvitt_g_proj]={4,9};
int list_weak_insq[nweak_ins]={1,2,3,4, 6,7,8,9,  1,2,3,4, 6,7,8,9};
int list_weak_insl[nweak_ins]={1,2,3,4, 6,7,8,9,  6,7,8,9, 1,2,3,4};
int nind;
//complex *glb_weak_corr;
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

//define types of quark propagator used
const int nins_kind=6;
enum insertion_t{                    ORIGINAL,  SCALAR,  PSEUDO,  STOCH_PHI,  STOCH_ETA,  TADPOLE};
const char ins_name[nins_kind][20]={"ORIGINAL","SCALAR","PSEUDO","STOCH_PHI","STOCH_ETA","TADPOLE"};
const int nqprop_kind=7;
enum qprop_t{                           PROP_0,  PROP_S,  PROP_P,  PROP_PHI,  PROP_ETA,  PROP_PHIETA,  PROP_T};
const char prop_name[nqprop_kind][20]={"PROP_0","PROP_S","PROP_P","PROP_PHI","PROP_ETA","PROP_PHIETA","PROP_T"};
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

//#define NOINSERT
//#define INSERT_ONLY_LEPT
//#define NOPHOTON
//#define ONLYTIME

#ifdef ONLYTIME
 int dirs[4]={1,0,0,0};
#else
 int dirs[4]={1,1,1,1};
#endif

//return appropriate propagator
int nqprop,nlprop;
int iqprop(int imass,qprop_t ip,int r)
{return r+nr*(imass+nqmass*ip);}
int ilprop(int ilepton,int orie,int phi_eta,int r)
{return r+nr*(phi_eta+2*(orie+2*ilepton));}

//return appropriately modified indo
tm_quark_info get_lepton_info(int ilepton,int orie,int r)
{
  tm_quark_info le=leps[ilepton];
  le.r=r;
  for(int i=1;i<4;i++) le.bc[i]*=sign_orie[orie];
  
  return le;
}

//generate a wall-source for stochastic QCD propagator
void generate_original_source()
{
  //Source coord
  //coords M={glb_size[0]/2,glb_size[1],glb_size[2],glb_size[3]};
  for(int mu=0;mu<4;mu++) source_coord[mu]=0;//(int)(rnd_get_unif(&glb_rnd_gen,0,1)*M[mu]);
  
#ifdef POINT_SOURCE_VERSION
  master_printf("Source position: t=%d x=%d y=%d z=%d\n",source_coord[0],source_coord[1],source_coord[2],source_coord[3]);
  generate_delta_source(original_source,source_coord);
#else
  master_printf("Source position: t=%d\n",source_coord[0]);
  generate_spindiluted_source(original_source,rnd_type_map[noise_type],source_coord[0]);
#endif
}

//generate a QED stochastic propagator
void generate_photon_stochastic_propagator()
{
  photon_prop_time-=take_time();
  generate_stochastic_tlSym_gauge_propagator(photon_phi,photon_eta,photon);
  // for(int ivol=0;ivol<loc_vol;ivol++)
  //   for(int mu=0;mu<4;mu++)
  //     {
  // 	complex_put_to_real(photon_phi[ivol][mu],1);
  // 	complex_put_to_real(photon_eta[ivol][mu],1);
  //     }
  // set_borders_invalid(photon_phi);
  // set_borders_invalid(photon_eta);
  photon_prop_time+=take_time();
  nphoton_prop_tot++;
}
//generate a sequential source
void generate_source(insertion_t inser,int r,PROP_TYPE *ori)
{
  source_time-=take_time();

  switch(inser)
    {
    case ORIGINAL:prop_multiply_with_gamma(source,0,original_source);break;
    case SCALAR:prop_multiply_with_gamma(source,0,ori);break;
    case PSEUDO:prop_multiply_with_gamma(source,5,ori);break;
#if !defined NOPHOTON &&! defined INSERT_ONLY_LEPT
    case STOCH_PHI:insert_external_source(source,conf,photon_phi,ori,r);break;
    case STOCH_ETA:insert_external_source(source,conf,photon_eta,ori,r);break;
#else
    case STOCH_PHI:insert_conserved_current(source,conf,ori,r,dirs);break;
    case STOCH_ETA:insert_conserved_current(source,conf,ori,r,dirs);break;
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
	
	//rotate the source index - please note that the propagator rotate AS the sign of mass term
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

void init_simulation(char *path)
{
  //////////////////////////// read the input /////////////////////////
  
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
	  master_printf(" supposed free quark energy[%d]: %lg\n",iq,free_quark_ener[iq]);
	}
      double free_mes_ener=free_quark_ener[0]+free_quark_ener[1];
      master_printf(" supposed free meson energy: %lg\n",free_mes_ener);
      
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
  
  //Discretization for photon propagator
  char photon_discrete_str[100];
  read_str_str("PhotonDiscretization",photon_discrete_str,100);
  if(strncasecmp(photon_discrete_str,"WILSON",100)==0) photon.c1=WILSON_C1;
  else
    if(strncasecmp(photon_discrete_str,"TLSYM",100)==0) photon.c1=TLSYM_C1;
    else crash("Unkwnown photon discretization: %s",photon_discrete_str);
  
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  //Read the noise type
  read_str_int("NoiseType",&noise_type);
  
  //Read if simulating in the free theory
  read_str_int("FreeTheory",&free_theory);
  
  //Read if simulating in the free theory
  read_str_int("RandomGaugeTransform",&rnd_gauge_transform);
  
  //Read the number of sources
  read_str_int("NSources",&nsources);
  
  //Number of configurations
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
  le_prop_ave=nissa_malloc("le_prop_ave",nleptons*2*2*nr*glb_size[0],complex);
  le_prop_err=nissa_malloc("le_prop_err",nleptons*2*2*nr*glb_size[0],complex);

  glb_corr=nissa_malloc("glb_corr",glb_size[0],complex);
  loc_corr=nissa_malloc("loc_corr",glb_size[0],complex);
  nind=nleptons*nweak_ins*2*2*nr;
  //glb_weak_corr=nissa_malloc("glb_weak_corr",glb_size[0]*nweak_ins*16*nind,complex);
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

//compute phase exponent: vec{p}*\vec{x}
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
  out[RE]=cos(arg)*ext;
  out[IM]=sin(arg)*ext;
}

//compute the phase for antineutrino - the orientation is that of the muon (as above)
void get_antineutrino_source_phase_factor(complex out,int ivol,int ilepton,momentum_t bc)
{
  //compute space and time factor
  double arg=get_space_arg(ivol,bc);
  int t=(glb_coord_of_loclx[ivol][0]-source_coord[0]+glb_size[0])%glb_size[0];
  if(t>=glb_size[0]/2) t=glb_size[0]-t;
  double ext=exp(t*neu_energy[ilepton]);
  
  //compute full exponential (notice the factor +1)
  out[RE]=cos(-arg)*ext;
  out[IM]=sin(-arg)*ext;
}

//set everything to a phase factor
void set_to_lepton_sink_phase_factor(spinspin *prop,int ilepton,tm_quark_info &le)
{
  GET_THREAD_ID();
  
  vector_reset(prop);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    if(glb_coord_of_loclx[ivol][0]==((glb_size[0]/2+0+source_coord[0])%glb_size[0]))
      {
	complex ph;
	get_lepton_sink_phase_factor(ph,ivol,ilepton,le);
	spinspin_put_to_diag(prop[ivol],ph);
      }
  set_borders_invalid(prop);
}

void trace_test_lep_prop(complex *c,spinspin *prop,tm_quark_info le,coords co)
{
  GET_THREAD_ID();

  if(IS_MASTER_THREAD)
    {
      //trace
      int imom=0;
      int par_apar=0;
      memset(c,0,sizeof(complex)*glb_size[0]);
      
      //compute steps
      momentum_t steps;
      for(int mu=1;mu<4;mu++)
	steps[mu]=-le.bc[mu]*M_PI/glb_size[mu];
      
      spinspin pr;
      twisted_particle_anti_particle_projector_of_imom(pr,le,imom,par_apar);
      NISSA_LOC_VOL_LOOP(ivol)
	{
	  //compute phase exponent
	  double arg=0;
	  for(int mu=1;mu<4;mu++)
	    arg+=steps[mu]*(glb_coord_of_loclx[ivol][mu]-co[mu]);
	  
	  //compute the phase
	  complex ph={cos(arg),sin(arg)};
	  
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
      double E=lep_energy[ilepton];

      //compute steps
      momentum_t steps;
      for(int mu=1;mu<4;mu++)
	steps[mu]=-le.bc[mu]*M_PI/glb_size[mu];
      
      spinspin pr;
      twisted_particle_anti_particle_projector_of_imom(pr,le,imom,par_apar);
      NISSA_LOC_VOL_LOOP(ivol)
	{
	  //compute phase exponent
	  double arg=0;
	  for(int mu=1;mu<4;mu++)
	    arg+=steps[mu]*glb_coord_of_loclx[ivol][mu];
	  
	  //compute pure time factor
	  double ext=exp(-(glb_size[0]-glb_coord_of_loclx[ivol][0])*E);
	  
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

void print_test_lep_prop(FILE *fout,complex *c,tm_quark_info le,const char *tag="")
{
  master_fprintf(fout," # mu=%lg kappa=%lg r=%d bc=%lg %s\n\n",le.mass,le.kappa,le.r,le.bc[1],tag);
  for(int t=0;t<glb_size[0];t++) master_fprintf(fout," %+016.16lg %+016.16lg\n",c[t][0],c[t][1]);
  master_fprintf(fout,"\n");
}

void test_lep_prop(FILE *fout,tm_quark_info le)
{
  GET_THREAD_ID();

  master_printf("energy: %+016.016lg\n",tm_quark_energy(le,0));
  
  //prepare the source
  spinspin *source=nissa_malloc("source",loc_vol,spinspin);
  vector_reset(source);
  //coords co={4,5,3,2};
  coords co={0,0,0,0};
  if(rank==0) spinspin_put_to_id(source[glblx_of_coord(co)]);
  
  //compute lepton prop in various ways
  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  complex c[glb_size[0]];
  
  //propagator by multiplication via fft
  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,source,le);
  trace_test_lep_prop(c,prop,le,co);
  print_test_lep_prop(fout,c,le,"multiply_from_right_by_x_space_twisted_propagator_by_fft");
  
  //propagator via fft
  //compute_x_space_twisted_propagator_by_fft(prop,le);
  //trace_test_lep_prop(c,prop,le);
  //print_test_lep_prop(fout,c,le,"compute_x_space_twisted_propagator_by_fft");
  
  //propagator via fft
  //compute_mom_space_twisted_propagator(prop,le);
  //for(int t=0;t<loc_size[0];t++)
  //{
  //complex r;
  //trace_dirac_prod_spinspin(r,base_gamma+0,prop[t*loc_spat_vol]);
  //int glb_t=t+glb_coord_of_loclx[0][0];
  //c[2*glb_t+0]=r[RE];
  //c[2*glb_t+1]=r[IM];
  //}
  //MPI_Allreduce(MPI_IN_PLACE,c,2*glb_size[0],MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //print_test_lep_prop(fout,c,le,"mom_space_repp");
  
  //propagator with scalar insertion
  spinspin *prop2=nissa_malloc("prop2",loc_vol,spinspin);

  //insert it explicitly
  //compute_x_space_twisted_propagator_by_fft(prop,le);
  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop2,prop,le);
  trace_test_lep_prop(c,prop2,le,co);
  print_test_lep_prop(fout,c,le,"explicit insertion");

  //take the derivative
  compute_x_space_twisted_propagator_by_fft(prop,le);
  tm_quark_info te_le=le;
  double eps=0.00001*le.mass;
  te_le.mass+=eps;
  compute_x_space_twisted_propagator_by_fft(prop2,te_le);
  complex c2[glb_size[0]];
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
       spinspin_subtassign(prop2[ivol],prop[ivol]);
       spinspin_prodassign_double(prop2[ivol],1/eps);
    }
  trace_test_lep_prop(c2,prop2,te_le,co);
  print_test_lep_prop(fout,c2,le,"numerical derivative");
  
  //insert gamma0 on the left
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      unsafe_dirac_prod_spinspin(prop2[ivol],base_gamma+4,prop[ivol]);
      for(int i=1;i<4;i++)
	{
	  spinspin t;
	  unsafe_dirac_prod_spinspin(t,base_gamma+i,prop[ivol]);
	  spinspin_summassign(prop2[ivol],t);
	}
    }
  set_borders_invalid(prop2);
  //multiply_from_left_by_x_space_twisted_propagator_by_fft(prop2,prop2,le);
  trace_test_lep_prop(c2,prop2,le,co);
  //for(int t=0;t<glb_size[0];t++) for(int ri=0;ri<2;ri++) c2[2*t+ri]/=c[2*t+ri]*std::min(t,glb_size[0]-t);
  print_test_lep_prop(fout,c2,le,"g0 insertion on the left");
  
  //insert gamma0 on the right
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    unsafe_spinspin_prod_dirac(prop2[ivol],prop[ivol],base_gamma+map_mu[0]);
  set_borders_invalid(prop2);
  //multiply_from_right_by_x_space_twisted_propagator_by_fft(prop2,prop2,le);
  trace_test_lep_prop(c2,prop2,le,co);
  //for(int t=0;t<glb_size[0];t++) for(int ri=0;ri<2;ri++) c2[2*t+ri]/=c[2*t+ri]*std::min(t,glb_size[0]-t);
  print_test_lep_prop(fout,c2,le,"g0 insertion on the right");
  
  {
    compute_mom_space_twisted_propagator(prop,le);

    vector_reset(prop2);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//takes the momenta part and M
	momentum_t p;
	for(int mu=0;mu<NDIM;mu++) p[mu]=M_PI*(2*glb_coord_of_loclx[ivol][mu]+le.bc[mu])/glb_size[mu];
	
	//fill the pieces
	spinspin op;
	spinspin_put_to_zero(op);
	for(int mu=0;mu<1;mu++)
	  {
	    spinspin_dirac_summ_the_prod_double(op,base_gamma+map_mu[mu],cos(p[mu]));
	    spinspin_dirac_summ_the_prod_double(op,&base_gamma[5],-sin(p[mu])*tau3[le.r]);
	  }
	
	unsafe_spinspin_prod_spinspin(prop2[ivol],op,prop[ivol]);
      }
    set_borders_invalid(prop2);

    pass_spinspin_from_mom_to_x_space(prop2,prop2,le.bc);

    spinspin_print(prop2[1000]);
    master_printf("prop2[1000] %d %d %d %d\n",glb_coord_of_loclx[1000][0],glb_coord_of_loclx[1000][1],glb_coord_of_loclx[1000][2],glb_coord_of_loclx[1000][3]);
    
    trace_test_lep_prop(c2,prop2,le,co);
    print_test_lep_prop(fout,c2,le,"conserved g0 multiplication on the left, fourier space");
  }

  //phases
  complex phases[4];
  for(int mu=0;mu<NDIM;mu++)
    {
      phases[mu][0]=cos(le.bc[mu]*M_PI);
      phases[mu][1]=sin(le.bc[mu]*M_PI);
    }
  
  compute_x_space_twisted_propagator_by_fft(prop,le);
  master_printf("bc: %lg %lg %lg %lg\n",le.bc[0],le.bc[1],le.bc[2],le.bc[3]);
  vector_reset(prop2);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<1;mu++)
      {
	/*find neighbors*/
	int ifw=loclx_neighup[ivol][mu];
	int ibw=loclx_neighdw[ivol][mu];

	/*transport down and up*/
	spinspin fw,bw;
	if(glb_coord_of_loclx[ivol][mu]==glb_size[mu]-1) unsafe_spinspin_prod_complex(fw,prop[ifw],phases[mu]);
	else spinspin_copy(fw,prop[ifw]);
	if(glb_coord_of_loclx[ivol][mu]==0) unsafe_spinspin_prod_complex_conj2(bw,prop[ibw],phases[mu]);
	else spinspin_copy(bw,prop[ibw]);

	/*fix fw and bw*/
	spinspin_prodassign_double(fw,-0.5);
	spinspin_prodassign_double(bw,+0.5);
	
	/*summ and subtract the two*/
	spinspin bw_M_fw,bw_P_fw;
	spinspin_subt(bw_M_fw,bw,fw);
	spinspin_summ(bw_P_fw,bw,fw);

	/*put -i g5 t3 on the summ*/
	spinspin g5_bw_P_fw;
	unsafe_dirac_prod_spinspin(g5_bw_P_fw,base_gamma+5,bw_P_fw);
	spinspin_summ_the_prod_idouble(prop2[ivol],g5_bw_P_fw,-tau3[le.r]);
	
	/*put gmu on the summ*/
	spinspin gmu_bw_M_fw;
	unsafe_dirac_prod_spinspin(gmu_bw_M_fw,base_gamma+map_mu[mu],bw_M_fw);
	spinspin_summassign(prop2[ivol],gmu_bw_M_fw);
      }
  set_borders_invalid(prop2);

  //multiply_from_left_by_x_space_twisted_propagator_by_fft(prop2,prop2,le);
  trace_test_lep_prop(c2,prop2,le,co);
  spinspin_print(prop2[1000]);
  master_printf("prop2[1000]\n");
  //for(int t=0;t<glb_size[0];t++) c2[t]/=c[t]*std::min(t,glb_size[0]-t);
  print_test_lep_prop(fout,c2,le,"conserved g0 multiplication on the left");

  compute_x_space_twisted_propagator_by_fft(prop,le);
  trace_test_lep_prop(c,prop,le,co);
  multiply_from_left_by_x_space_twisted_propagator_by_fft(prop2,prop2,le);
  trace_test_lep_prop(c2,prop2,le,co);
  for(int t=0;t<glb_size[0];t++)
    for(int ri=0;ri<2;ri++) c2[t][ri]/=c[t][0]*t;
  print_test_lep_prop(fout,c2,le,"conserved g0 insertion on the left");


  //prepare each propagator for a single lepton
  //by computing i(phi(x-mu)A_mu(x-mu)(-i t3 g5-gmu)/2-phi(x+mu)A_mu(x)(-i t3 g5+gmu)/2)=
  //(ph0 A_mu(x-mu)g[r][0][mu]-ph0 A_mu(x)g[r][1][mu])=
  vector_reset(prop2);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<1;mu++)
      {
	//find neighbors
	int ifw=loclx_neighup[ivol][mu];
	int ibw=loclx_neighdw[ivol][mu];
	
	//compute phase factor
	spinspin ph_bw,ph_fw;
	//transport down and up
	if(glb_coord_of_loclx[ivol][mu]==glb_size[mu]-1) unsafe_spinspin_prod_complex(ph_fw,prop[ifw],phases[mu]);
	else spinspin_copy(ph_fw,prop[ifw]);
	if(glb_coord_of_loclx[ivol][mu]==0) unsafe_spinspin_prod_complex_conj2(ph_bw,prop[ibw],phases[mu]);
	else spinspin_copy(ph_bw,prop[ibw]);

	//fix coefficients - i is inserted here!
	spinspin_prodassign_double(ph_fw,-0.5);
	spinspin_prodassign_double(ph_bw,+0.5);
	
	//summ and subtract the two
	spinspin bw_M_fw,bw_P_fw;
	spinspin_subt(bw_M_fw,ph_bw,ph_fw);
	spinspin_summ(bw_P_fw,ph_bw,ph_fw);

	//put -i g5 t3 on the summ
	spinspin_prodassign_idouble(bw_P_fw,-tau3[le.r]);
	spinspin temp;
	unsafe_spinspin_prod_dirac(temp,bw_P_fw,base_gamma+5);
	spinspin_summassign(prop2[ivol],temp);
	
	//put gmu on the diff
	unsafe_spinspin_prod_dirac(temp,bw_M_fw,base_gamma+map_mu[mu]);
	spinspin_summassign(prop2[ivol],temp);
      }
  set_borders_invalid(prop2);
  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop2,prop2,le);
  trace_test_lep_prop(c2,prop2,le,co);
  //for(int t=0;t<glb_size[0];t++)
  //for(int ri=0;ri<2;ri++) c2[t][ri]/=c[t][0];
  print_test_lep_prop(fout,c2,le,"conserved g0 insertion on the right");
  
  nissa_free(source);
  nissa_free(prop);
  nissa_free(prop2);
}

//insert the photon on the source
void insert_photon_on_the_source(spinspin *prop,int ilepton,int phi_eta,coords dirs,tm_quark_info &le)
{ 
  GET_THREAD_ID();
  
  //select A
#ifndef NOPHOTON
  spin1field *A=(phi_eta==0)?photon_phi:photon_eta;
  communicate_lx_spin1field_borders(A);
#endif
	    
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
	if(glb_coord_of_loclx[ivol][mu]==glb_size[mu]-1) unsafe_spinspin_prod_complex(ph_fw,temp_lep[ifw],phases[mu]);
	else spinspin_copy(ph_fw,temp_lep[ifw]);
	if(glb_coord_of_loclx[ivol][mu]==0) unsafe_spinspin_prod_complex_conj2(ph_bw,temp_lep[ibw],phases[mu]);
	else spinspin_copy(ph_bw,temp_lep[ibw]);
	
	//fix coefficients - i is inserted here!
	//also dir selection is made here
#ifndef NOPHOTON
	spinspin_prodassign_idouble(ph_fw,-0.5*dirs[mu]);
	spinspin_prodassign_idouble(ph_bw,+0.5*dirs[mu]);
#else
	spinspin_prodassign_double(ph_fw,-0.5*dirs[mu]);
	spinspin_prodassign_double(ph_bw,+0.5*dirs[mu]);
#endif
	
	//fix insertion of the current
#ifndef NOPHOTON
	safe_spinspin_prod_complex(ph_fw,ph_fw,A[ivol][mu]);
	safe_spinspin_prod_complex(ph_bw,ph_bw,A[ibw][mu]);
#endif
	
	//summ and subtract the two
	spinspin bw_M_fw,bw_P_fw;
	spinspin_subt(bw_M_fw,ph_bw,ph_fw);
	spinspin_summ(bw_P_fw,ph_bw,ph_fw);
		  
	//put -i g5 t3 on the summ
	spinspin_prodassign_idouble(bw_P_fw,-tau3[le.r]);
	spinspin temp;
	unsafe_spinspin_prod_dirac(temp,bw_P_fw,base_gamma+5);
	//spinspin_summassign(prop[ivol],temp);
	
	//put gmu on the diff
	unsafe_spinspin_prod_dirac(temp,bw_M_fw,base_gamma+map_mu[mu]);
	spinspin_summassign(prop[ivol],temp);
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

  int i=0;
  for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int orie=0;orie<2;orie++)
      for(int phi_eta=0;phi_eta<2;phi_eta++)
	for(int r=0;r<nr;r++)
	  {
	    //set the properties of the meson
	    //time boundaries are anti-periodic, space are as for external line
	    tm_quark_info le=get_lepton_info(ilepton,orie,r);
	    
	    //test_lep_prop(fout,le); 
	    //crash("");
	    
	    //select the propagator
	    int iprop=ilprop(ilepton,orie,phi_eta,r);
	    spinspin *prop=L[iprop];

	    //put it to a phase
	    set_to_lepton_sink_phase_factor(prop,ilepton,le);

	    //insert the current
#if !defined NOINSERT
	    multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le);
	    insert_photon_on_the_source(prop,ilepton,phi_eta,dirs,le);
#endif
	    //multiply by the lepton propagator (from the right)
	    multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le);
#if !defined NOINSERT && defined INSERT_ONLY_LEPT
	    insert_photon_on_the_source(prop,ilepton,!phi_eta,dirs,le);
#endif
	    
	    complex test[glb_size[0]];
	    trace_test_lep_prop_source(test,prop,le,ilepton);
	    for(int t=0;t<glb_size[0];t++)
	      {
		complex te;
		complex_prod_double(te,test[t],1.0/glb_vol*glb_size[0]);
		complex_summassign(le_prop_ave[i],te);
		complex te2={sqr(te[0]),sqr(te[1])};
		complex_summassign(le_prop_err[i],te2);
		i++;
	      }
	  }
  
  if(IS_MASTER_THREAD) lepton_prop_time+=take_time();
}
THREADABLE_FUNCTION_END

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
  
  if(rnd_gauge_transform) perform_random_gauge_transform(conf,conf);
  
  //put anti-periodic boundary condition for the fermionic propagator
  old_theta[0]=old_theta[1]=old_theta[2]=old_theta[3]=0;
  put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
  adapt_theta(conf,old_theta,put_theta,0,0);

  //reset correlations
  //vector_reset(glb_weak_corr);
  vector_reset(hadr_ave);
  vector_reset(hadr_err);
  vector_reset(le_prop_ave);
  vector_reset(le_prop_err);
  vector_reset(glb_weak_vitt_ave);
  vector_reset(glb_weak_vitt_err);
  vector_reset(glb_weak_proj_corr);
}

//compute the hadronic part of the lepton correlation function
//as usual, FIRST propagator is reverted
THREADABLE_FUNCTION_3ARG(hadronic_part_leptonic_correlation, spinspin*,hadr, PROP_TYPE*,S1, PROP_TYPE*,S2)
{
  GET_THREAD_ID();
  
  vector_reset(hadr);

  //it's just the matter of inserting gamma5*gamma5=identity between S1^dag and S2, g5 on the sink index of S1 and trace on them
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int ic=0;ic<3;ic++)
#ifdef POINT_SOURCE_VERSION
      for(int jc=0;jc<3;jc++)
#endif
	for(int id=0;id<4;id++)
	  for(int jd=0;jd<4;jd++)
	    for(int kd=0;kd<4;kd++)
	      ((id<2)?complex_summ_the_conj1_prod:complex_subt_the_conj1_prod)
		(hadr[ivol][id][jd],
#ifdef POINT_SOURCE_VERSION
					    S1[ivol][ic][jc][jd][kd],S2[ivol][ic][jc][id][kd]);
#else
		                            S1[ivol][ic][jd][kd],S2[ivol][ic][id][kd]);
#endif
  
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

//compute the polarisation vector
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
THREADABLE_FUNCTION_6ARG(compute_leptonic_correlation, spinspin*,hadr, int,iprop, int,ilepton, int,orie, int,rl, int,ind)
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
  
  for(int smu=0;smu<2;smu++)
    for(int snu=0;snu<2;snu++)
      {
	spin te;
	unsafe_dirac_prod_spin(te,base_gamma+4,ubar[smu]);
	complex hl={0,0};
	for(int id_si=0;id_si<4;id_si++) complex_summ_the_prod(hl,v[snu][id_si],te[id_si]);
	//master_printf("smu=%d snu=%d proj[%s]={%+016.016lg,%+016.16lg}\n",smu,snu,gtag[4],hl[0],hl[1]);
      }
  
  //gets the projectors for Vittorio
  spinspin promu[2],pronu[2];
  twisted_on_shell_operator_of_imom(promu[0],le,0,false,-1);
  twisted_on_shell_operator_of_imom(promu[1],le,0,false,+1);
  naive_massless_on_shell_operator_of_imom(pronu[0],le.bc,0,-1);
  naive_massless_on_shell_operator_of_imom(pronu[1],le.bc,0,+1);

  //compute the right part of the leptonic loop
  dirac_matr vitt_proj_gamma[nvitt_g_proj];
  for(int ig_proj=0;ig_proj<nvitt_g_proj;ig_proj++)
    {
      int ig=vitt_g_projs[ig_proj];
      dirac_matr temp_gamma;
      dirac_herm(&temp_gamma,base_gamma+ig);
      dirac_prod(vitt_proj_gamma+ig_proj,base_gamma+map_mu[0],&temp_gamma);
    }
  //insert gamma5 on the sink-hadron-gamma - check order
  dirac_matr weak_ins_hadr_gamma[nweak_ins];
  for(int ins=0;ins<nweak_ins;ins++) dirac_prod(weak_ins_hadr_gamma+ins,base_gamma+5,base_gamma+list_weak_insl[ins]);
  
  //debug
  //twisted_on_shell_operator_of_imom(promu[0],le,0,false,-1);
  //twisted_on_shell_operator_of_imom(promu[1],le,0,false,+1);
  //naive_massless_on_shell_operator_of_imom(pronu[0],le.bc,0,+1);
  //naive_massless_on_shell_operator_of_imom(pronu[1],le.bc,0,-1);

  //master_printf("-\n");
  //spinspin out[2];
  //for(int mu=0;mu<2;mu++) spin_direct_prod(out[mu],ubar[mu],u[mu]);
  //spinspin_summassign(out[0],out[1]);
  //spinspin_print(out[0]);
  //master_printf("-\n");
  //spinspin_print(promu[0]);
  //master_printf("-\n");
  //spinspin_copy(promu[0],out[0]);

  //master_printf("-\n");
  //for(int nu=0;nu<2;nu++) spin_direct_prod(out[nu],vbar[nu],v[nu]);
  //spinspin_summassign(out[0],out[1]);
  //spinspin_print(out[0]);
  //master_printf("-\n");
  //spinspin_print(pronu[0]);
  //master_printf("-\n");
  //spinspin_copy(pronu[0],out[0]);
  
  //crash("");
  
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
	  
	  spinspin_summ_the_complex_prod(hl_loc_corr[t],l,h);
	}
      glb_threads_reduce_double_vect((double*)hl_loc_corr,loc_size[0]*sizeof(spinspin)/sizeof(double));
      
      //change sign when crossing 0 for averaging corr function properly
      for(int loc_t=0;loc_t<loc_size[0];loc_t++)
	{
	  int glb_t=loc_t+glb_coord_of_loclx[0][0];
	  int sign=1-2*(glb_t<source_coord[0]);
	  //if(rank==0 && thread_id==NACTIVE_THREADS-1) printf("ANNA SIGN %d %d %lg\n",glb_t,sign,hl_loc_corr[loc_t][0][0][RE]);
	  spinspin_prodassign_double(hl_loc_corr[loc_t],sign);
	}
      
      //reduce
      // for(int ig=0;ig<16;ig++)
      // 	for(int loc_t=0;loc_t<loc_size[0];loc_t++)
      // 	  {
      // 	    int glb_t=(loc_t+glb_coord_of_loclx[0][0]-source_coord[0]+glb_size[0])%glb_size[0];
      // 	    complex hl;
      // 	    trace_spinspin_with_dirac(hl,hl_loc_corr[loc_t],base_gamma+ig);
      // 	    if(IS_MASTER_THREAD) complex_summassign(glb_weak_corr[glb_t+glb_size[0]*(ig+16*(ins+nweak_ins*ind))],hl);
      // 	  }
      // 	if(IS_MASTER_THREAD) nlept_contr_tot+=16;

	//do it using the spinor
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

	//do it again with vittorio
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
	      if(irev==1) std::swap(iq1,iq2);
	      //takes the propagators
#ifndef NOINSERT
	      int ip1=iqprop(iq1,(qins==0)?((phi_eta==0)?PROP_PHI:PROP_ETA):PROP_0,r2); //q1 will be reverted
	      int ip2=iqprop(iq2,(qins==1)?((phi_eta==0)?PROP_PHI:PROP_ETA):PROP_0,r2);
#else
	      int ip1=iqprop(iq1,PROP_0,r2); //q1 will be reverted
	      int ip2=iqprop(iq2,PROP_0,r2);
#endif
	      //compute the hadronic part
	      hadronic_part_leptonic_correlation(hadr,Q[ip1],Q[ip2]);
	      
	      for(int orie=0;orie<2;orie++)
		for(int rl=0;rl<nr;rl++)
		  {
		    //contract with lepton
		    int iprop=ilprop(ilepton,orie,!phi_eta,rl);
		    compute_leptonic_correlation(hadr,iprop,ilepton,orie,rl,ind);
		    ind++;
		  }
	    }

  //adds id to the stack
  glb_nodes_reduce_complex_vect(glb_weak_vitt_corr,glb_size[0]*nweak_ins*nvitt_g_proj*nind);
  for(int iel=0;iel<glb_size[0]*nweak_ins*nvitt_g_proj*nind;iel++)
    {
      complex_summassign(glb_weak_vitt_ave[iel],glb_weak_vitt_corr[iel]);
      complex temp={sqr(glb_weak_vitt_corr[iel][RE]),sqr(glb_weak_vitt_corr[iel][IM])};
      complex_summassign(glb_weak_vitt_err[iel],temp);
    }
  lept_contr_time+=take_time();
}

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
      complex_prodassign_double(hadr_ave[i],n);
      complex_prodassign_double(hadr_err[i],n);
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

//print a whole weak correlation function
// void print_weak_correlations(FILE *fout,int ind)
// {
//   for(int ins=0;ins<nweak_ins;ins++)
//     for(int ig=0;ig<16;ig++)
//       {
// 	master_fprintf(fout," # qins=%d lins=%d proj=%d\n\n",list_weak_insq[ins],list_weak_insl[ins],ig);
// 	for(int t=0;t<glb_size[0];t++)
// 	  {
// 	    int i=t+glb_size[0]*(ig+16*(ins+nweak_ins*ind));
// 	    master_fprintf(fout,"%+016.16lg %+016.16lg\n",glb_weak_corr[i][RE]/nsources,glb_weak_corr[i][IM]/nsources);
// 	  }
// 	master_fprintf(fout,"\n");
//       }
// }

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

void print_le_prop()
{
  FILE *fout[2]={open_file(combine("%s/le_prop_ave",outfolder).c_str(),"w"),
		 open_file(combine("%s/le_prop_err",outfolder).c_str(),"w")};
  
  int i=0;
  for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int orie=0;orie<2;orie++)
      for(int phi_eta=0;phi_eta<2;phi_eta++)
	for(int r=0;r<nr;r++)
	  {
	    for(int ifi=0;ifi<2;ifi++)
	      {
		master_fprintf(fout[ifi]," # ilepton=%d orie=%d r=%d",ilepton,orie,r);
		master_fprintf(fout[ifi],"\n");
	      }
	    for(int t=0;t<glb_size[0];t++)
	      {
		double A=le_prop_ave[i][RE]/nsources;
		double A2=le_prop_err[i][RE]/nsources;
		double B=le_prop_ave[i][IM]/nsources;
		double B2=le_prop_err[i][IM]/nsources;
		A2-=A*A;
		B2-=B*B;
		if(nsources>1)
		  {
		    A2/=nsources-1;
		    B2/=nsources-1;
		  }
		A2=sqrt(fabs(A2));
		B2=sqrt(fabs(B2));
		master_fprintf(fout[0],"%d %+016.016lg %+016.016lg\n",t,A,B);
		master_fprintf(fout[1],"%d %+016.016lg %+016.016lg\n",t,A2,B2);
		i++;
	      }
	    for(int ifi=0;ifi<2;ifi++) master_fprintf(fout[ifi],"\n");
	  }
  
  for(int ifi=0;ifi<2;ifi++) close_file(fout[ifi]);
}

//print out everything
void print_outputs()
{
  print_time-=take_time();
  
  //final reduction
  // glb_nodes_reduce_complex_vect(glb_weak_corr,glb_size[0]*nweak_ins*16*nind);
  glb_nodes_reduce_complex_vect(glb_weak_proj_corr,glb_size[0]*nweak_ins*4*nind);
  
  //open the three outputs
  // const int nwrite=1;
  // FILE *fout[1]={open_file(combine("%s/corr_hl_vitt",outfolder).c_str(),"w")};

  const int nwrite=3;
  FILE *fout[3]={
    open_file(combine("%s/corr_hl_vitt",outfolder).c_str(),"w"),
    open_file(combine("%s/err_hl_vitt",outfolder).c_str(),"w"),
    open_file(combine("%s/corr_hl_proj",outfolder).c_str(),"w")};

  // FILE *fout[3]={
  //   open_file(combine("%s/corr_hl_vitt",outfolder).c_str(),"w"),
  //   open_file(combine("%s/corr_hl_proj",outfolder).c_str(),"w"),
  //   open_file(combine("%s/corr_hl",outfolder).c_str(),"w")};
  
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
		    // print_weak_correlations(fout[2],ind);
		    ind++;
		  }
	    }
  for(int i=0;i<nwrite;i++) close_file(fout[i]);

  print_hadronic_correlations();
  print_le_prop();
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
  nissa_free(le_prop_ave);
  nissa_free(le_prop_err);
  nissa_free(glb_corr);
  nissa_free(loc_corr);
  // nissa_free(glb_weak_corr);
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
