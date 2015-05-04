#include <nissa.hpp>

#ifdef POINT_SOURCE_VERSION
 #define PROP_TYPE su3spinspin
#else
 #define PROP_TYPE colorspinspin
#endif

using namespace nissa;

/////////////////////////////////////// data //////////////////////////////

int ninv_tot=0,nhadr_contr_tot=0,nlept_contr_tot=0,nsource_tot=0,nphoton_prop_tot=0;
double inv_time=0,hadr_contr_time=0,lept_contr_time=0;
double tot_prog_time=0,source_time=0,photon_prop_time=0,lepton_prop_time=0;

int wall_time;
int free_theory;
int ngauge_conf,nanalyzed_conf=0;
char conf_path[1024],outfolder[1024];
quad_su3 *conf;

double kappa;
double put_theta[4],old_theta[4]={0,0,0,0};

coords source_coord;
PROP_TYPE *source,*original_source;
int seed,noise_type;

int nqmass,nr;
double *qmass,*residue;
PROP_TYPE **Q;

spincolor *temp_source;
spincolor *temp_solution;
  
gauge_info photon;
double tadpole[4];
spin1field *photon_phi,*photon_eta;

complex *glb_corr,*loc_corr;

//list the 8 matrices to insert for the weak current
const int nweak_ins=16;
int list_weak_insq[nweak_ins]={1,2,3,4, 6,7,8,9,  1,2,3,4, 6,7,8,9};
int list_weak_insl[nweak_ins]={1,2,3,4, 6,7,8,9,  6,7,8,9, 1,2,3,4};
complex *glb_weak_corr;

const int V_curr_mu=0;

//define types of quark propagator used
const int nins_kind=6;
enum insertion_t{                    ORIGINAL,  SCALAR,  PSEUDO,  STOCH_PHI,  STOCH_ETA,  TADPOLE};
const char ins_name[nins_kind][20]={"ORIGINAL","SCALAR","PSEUDO","STOCH_PHI","STOCH_ETA","TADPOLE"};
const int nqprop_kind=7;
enum qprop_t{                           PROP_0,  PROP_S,  PROP_P,  PROP_A,  PROP_B,  PROP_AB,  PROP_T};
const char prop_name[nqprop_kind][20]={"PROP_0","PROP_S","PROP_P","PROP_A","PROP_B","PROP_AB","PROP_T"};
//map the source, the destination and the insertion for each propagator
const qprop_t prop_map[nqprop_kind]=         {PROP_0,   PROP_S, PROP_P, PROP_A,    PROP_B,    PROP_AB,   PROP_T};
const insertion_t insertion_map[nqprop_kind]={ORIGINAL, SCALAR, PSEUDO, STOCH_PHI, STOCH_ETA, STOCH_ETA, TADPOLE};
const qprop_t source_map[nqprop_kind]=       {PROP_0,   PROP_0, PROP_0, PROP_0,    PROP_0,    PROP_A,    PROP_0};
const char prop_abbr[]=                       "0"       "S"     "P"     "A"        "B"        "X"        "T";

//parameters of the leptons
int nleptons;
int *lep_corr_iq1;
int *lep_corr_iq2;
int lepton_mom_sign[2]={-1,+1};
tm_quark_info *leps;
double *lep_energy,*lep_mom;
spinspin **L;

//return appropriate propagator
int nqprop,nlprop;
int iqprop(int imass,qprop_t ip,int r)
{return r+nr*(imass+nqmass*ip);}
int ilprop(int ilepton,int orie,int phi_eta,int r)
{return r+nr*(phi_eta+2*(orie+2*ilepton));}

//generate a wall-source for stochastic QCD propagator
void generate_original_source()
{
#ifdef POINT_SOURCE_VERSION
  generate_delta_source(original_source,source_coord);
#else
  generate_spindiluted_source(original_source,rnd_type_map[noise_type],source_coord[0]);
#endif
}

//generate a QED stochastic propagator
void generate_photon_stochastic_propagator()
{
  photon_prop_time-=take_time();
  generate_stochastic_tlSym_gauge_propagator(photon_phi,photon_eta,photon);
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
    case STOCH_PHI:insert_external_source(source,conf,photon_phi,ori,r);break;
    case STOCH_ETA:insert_external_source(source,conf,photon_eta,ori,r);break;
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
  lep_mom=nissa_malloc("lep_mom",nleptons,double);
  expect_str("Q1Q2LepmassMesmass");
  for(int il=0;il<nleptons;il++)
    {
      //read quarks identfiying the mesons, and lepton mass
      read_int(lep_corr_iq1+il);
      read_int(lep_corr_iq2+il);
      read_double(&leps[il].mass);

      //maximal twist
      leps[il].kappa=0.125;
      leps[il].r=0;
      
      //read the mass of the meson (that must have been determined outside)
      double mes_mass;
      read_double(&mes_mass);

      //compute meson momentum and bc
      double lep_mom_mod=(pow(mes_mass,2)-pow(leps[il].mass,2))/(2*mes_mass);
      lep_mom[il]=lep_mom_mod/sqrt(3);
      leps[il].bc[0]=0;
      for(int i=1;i<4;i++) leps[il].bc[i]=lep_mom[il]/M_PI*glb_size[i];
      
      //write down energy
      lep_energy[il]=tm_quark_energy(leps[il],0);
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
  glb_corr=nissa_malloc("glb_corr",glb_size[0],complex);
  loc_corr=nissa_malloc("loc_corr",glb_size[0],complex);
  glb_weak_corr=nissa_malloc("glb_weak_corr",glb_size[0]*nweak_ins*16,complex);
  original_source=nissa_malloc("source",loc_vol,PROP_TYPE);
  source=nissa_malloc("source",loc_vol,PROP_TYPE);
  photon_eta=nissa_malloc("photon_eta",loc_vol+bord_vol,spin1field);
  photon_phi=nissa_malloc("photon_phi",loc_vol+bord_vol,spin1field);
  Q=nissa_malloc("Q*",nqprop,PROP_TYPE*);
  for(int iprop=0;iprop<nqprop;iprop++) Q[iprop]=nissa_malloc("Q",loc_vol+bord_vol,PROP_TYPE);
  L=nissa_malloc("L*",nlprop,spinspin*);
  for(int iprop=0;iprop<nlprop;iprop++) L[iprop]=nissa_malloc("L",loc_vol+bord_vol,spinspin);
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
      
      //Source coord
      read_int(&(source_coord[0]));
#ifdef POINT_SOURCE_VERSION
      for(int mu=1;mu<4;mu++) read_int(&(source_coord[mu]));
#endif

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

//compute the phase for lepton
void get_lepton_phase_factor(complex out,int ivol,int ilepton,int orie,tm_quark_info le)
{
  //fetch the energy and momentum ot the muon
  double E=lep_energy[ilepton];
  double p=lep_mom[ilepton];
  
  //compute pure time factor and summ the coordinate
  double ext=exp(glb_coord_of_loclx[ivol][0]*E);
  int summ_coord=0;
  for(int i=1;i<NDIM;i++) summ_coord+=glb_coord_of_loclx[ivol][i];
  summ_coord*=lepton_mom_sign[orie];
    
  //compute full exponential and multiply each component of the photon field
  out[RE]=cos(-summ_coord*p)*ext;
  out[IM]=sin(-summ_coord*p)*ext;
}

void trace_test_lep_prop(double *c,spinspin *prop,tm_quark_info le,int reim=RE)
{
  GET_THREAD_ID();
  
  //trace
  memset(c,0,sizeof(double)*glb_size[0]);

  if(IS_MASTER_THREAD)
    {
      NISSA_LOC_VOL_LOOP(ivol)
	{
	  //compute the phase
	  double ph=0;
	  for(int mu=1;mu<4;mu++) ph+=-M_PI*glb_coord_of_loclx[ivol][mu]*le.bc[mu]/glb_size[mu];

	  complex t;
	  trace_dirac_prod_spinspin(t,base_gamma+0,prop[ivol]);
	  //c[glb_coord_of_loclx[ivol][0]]+=sqr(cos(ph))*spinspin_norm2(prop[ivol]);
	  c[glb_coord_of_loclx[ivol][0]]+=t[reim];
	}
      MPI_Allreduce(MPI_IN_PLACE,c,glb_size[0],MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
  THREAD_BARRIER();
}

void print_test_lep_prop(FILE *fout,double *c,tm_quark_info le,const char *tag="")
{
  master_fprintf(fout," # mu=%lg kappa=%lg r=%d bc=%lg %s\n\n",le.mass,le.kappa,le.r,le.bc[1],tag);
  for(int t=0;t<glb_size[0];t++) master_fprintf(fout," %+016.16lg\n",c[t]);
  master_fprintf(fout,"\n");
}

void test_lep_prop(FILE *fout,tm_quark_info le)
{
  GET_THREAD_ID();
  
  //prepare the source
  spinspin *source=nissa_malloc("source",loc_vol,spinspin);
  vector_reset(source);
  if(rank==0) spinspin_put_to_id(source[0]);
  
  //compute lepton prop in various ways
  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  double c[glb_size[0]];
  
  //propagator by multiplication via fft
  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,source,le);
  trace_test_lep_prop(c,prop,le);
  print_test_lep_prop(fout,c,le,"multiply_from_right_by_x_space_twisted_propagator_by_fft");

  //propagator via fft
  compute_x_space_twisted_propagator_by_fft(prop,le);
  trace_test_lep_prop(c,prop,le);
  print_test_lep_prop(fout,c,le,"compute_x_space_twisted_propagator_by_fft");

  //compute_x_space_twisted_propagator_by_inv(prop,le);
  //trace_test_lep_prop(c,prop,le);
  //print_test_lep_prop(fout,c,"multiply_from_right_by_x_space_twisted_propagator_by_fft");

  //propagator with scalar insertion
  spinspin *prop2=nissa_malloc("prop2",loc_vol,spinspin);

  //insert it explicitly
  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop2,prop,le);
  trace_test_lep_prop(c,prop2,le);
  print_test_lep_prop(fout,c,le,"explicit insertion");

  //take the derivative
  tm_quark_info te_le=le;
  double eps=0.001*le.mass;
  te_le.mass+=eps;
  compute_x_space_twisted_propagator_by_fft(prop2,te_le);
  trace_test_lep_prop(c,prop,le);
  double c2[glb_size[0]];
  trace_test_lep_prop(c2,prop2,te_le);
  for(int t=0;t<glb_size[0];t++) c2[t]=(c[t]-c2[t])/eps;
  print_test_lep_prop(fout,c2,le,"numerical derivative");
  
  //insert gamma0 on the left
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    unsafe_dirac_prod_spinspin(prop2[ivol],base_gamma+4,prop[ivol]);
  set_borders_invalid(prop2);
  multiply_from_left_by_x_space_twisted_propagator_by_fft(prop2,prop2,le);
  trace_test_lep_prop(c2,prop2,le);
  for(int t=0;t<glb_size[0];t++) c2[t]/=c[t]*t;
  print_test_lep_prop(fout,c2,le,"g0 insertion on the left");
  
  //insert gamma0 on the right
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    unsafe_spinspin_prod_dirac(prop2[ivol],prop[ivol],base_gamma+4);
  set_borders_invalid(prop2);
  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop2,prop2,le);
  trace_test_lep_prop(c2,prop2,le);
  for(int t=0;t<glb_size[0];t++) c2[t]/=c[t]*t;
  print_test_lep_prop(fout,c2,le,"g0 insertion on the right");

  //conserved current on the left
  spinspin g[2][2][4];
  double mps[2]={-0.5,+0.5}; //signs inverted?
  for(int r=0;r<2;r++)
    for(int mp=0;mp<2;mp++)
      for(int mu=0;mu<4;mu++)
	{
	  spinspin_dirac_prod_double(g[r][mp][mu],base_gamma+5,0.5*tau3[r]);
	  spinspin_dirac_summ_the_prod_idouble(g[r][mp][mu],base_gamma+mu,mps[mp]);
	}
  vector_reset(prop2);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      {
	int ifw=loclx_neighup[ivol][mu];
	int ibw=loclx_neighdw[ivol][mu];
	spinspin_summ_the_spinspin_prod(prop2[ivol],g[le.r][0][mu],prop[ibw]);
	spinspin_summ_the_spinspin_prod(prop2[ivol],g[le.r][1][mu],prop[ifw]);
      }
  set_borders_invalid(prop2);
  multiply_from_left_by_x_space_twisted_propagator_by_fft(prop2,prop2,le);
  trace_test_lep_prop(c2,prop2,le,1);
  for(int t=0;t<glb_size[0];t++) c2[t]/=c[t]*t;
  print_test_lep_prop(fout,c2,le,"conserved g0 insertion on the left");
  
  nissa_free(source);
  nissa_free(prop);
  nissa_free(prop2);
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

  FILE *fout=open_file(combine("%s/le_prop",outfolder).c_str(),"w");
  
  //communicate borders of photon
  communicate_lx_spin1field_borders(photon_eta);
  communicate_lx_spin1field_borders(photon_phi);
  
  //for the time being the boundaries for internal propagator are periodic
  tm_quark_info le;
  le.kappa=0.125;
  for(int mu=0;mu<4;mu++) le.bc[mu]=0;
  
  //compute the 4 pairs of gammas: i(-i t3 g5-+gmu)/2=(t3 g5-+i gmu)/2
  spinspin g[2][2][4]; //[t3][-+][mu]
  double mps[2]={-0.5,+0.5};
  for(int r=0;r<2;r++)
    for(int mp=0;mp<2;mp++)
      for(int mu=0;mu<4;mu++)
	{
	  spinspin_dirac_prod_double(g[r][mp][mu],base_gamma+5,0.5*tau3[r]);
	  spinspin_dirac_summ_the_prod_idouble(g[r][mp][mu],base_gamma+mu,mps[mp]);
	}
  
  for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int orie=0;orie<2;orie++)
      for(int phi_eta=0;phi_eta<2;phi_eta++)
	for(int r=0;r<nr;r++)
	  {
	    //select the propagator and reset it
	    int iprop=ilprop(ilepton,orie,phi_eta,r);
	    vector_reset(L[iprop]);
	    
	    //select A and fix lepton energy
	    spin1field *A=(phi_eta==0)?photon_phi:photon_eta;
	    le.mass=leps[ilepton].mass;
	    le.r=r;
	    
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
		  complex ph0,ph1;
		  get_lepton_phase_factor(ph0,ibw,ilepton,orie,le);
		  get_lepton_phase_factor(ph1,ifw,ilepton,orie,le);
		  
		  //multiply by A
		  complex phA0,phA1;
		  unsafe_complex_prod(phA0,ph0,A[ibw][mu]);
		  unsafe_complex_prod(phA1,ph0,A[ifw][mu]);
		  
		  //summ the two pieces
		  spinspin_summ_the_complex_prod(L[iprop][ivol],g[r][0][mu],phA0);
		  spinspin_subt_the_complex_prod(L[iprop][ivol],g[r][1][mu],phA1);
		}
	    set_borders_invalid(L[iprop]);

	    test_lep_prop(fout,le);
	    
	    //multiply by the lepton propagator (from the right)
	    multiply_from_right_by_x_space_twisted_propagator_by_fft(L[iprop],L[iprop],le);
	  }

  close_file(fout);
  
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
      master_printf("plaq: %.16g\n",global_plaquette_lx_conf(conf));
    }
  else  generate_cold_lx_conf(conf);

  //put anti-periodic boundary condition for the fermionic propagator
  old_theta[0]=old_theta[1]=old_theta[2]=old_theta[3]=0;
  put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
  adapt_theta(conf,old_theta,put_theta,0,0);
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

  nissa_free(photon_eta);
  nissa_free(photon_phi);
  nissa_free(source);
  nissa_free(original_source);
  for(int iprop=0;iprop<nqprop;iprop++) nissa_free(Q[iprop]);
  nissa_free(Q);
  for(int iprop=0;iprop<nlprop;iprop++) nissa_free(L[iprop]);
  nissa_free(L);
  nissa_free(conf);
  nissa_free(glb_corr);
  nissa_free(loc_corr);
  nissa_free(glb_weak_corr);
  nissa_free(temp_source);
  nissa_free(temp_solution);
  nissa_free(lep_corr_iq1);
  nissa_free(lep_corr_iq2);
  nissa_free(leps);
  nissa_free(lep_energy);
  nissa_free(lep_mom);
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
					    S1[ivol][ic][jc][kd][id],S2[ivol][ic][jc][kd][jd]);
#else
		                            S1[ivol][ic][kd][id],S2[ivol][ic][kd][jd]);
#endif
  
  THREAD_BARRIER();
}
THREADABLE_FUNCTION_END

//compute the full leptonic correlation function
THREADABLE_FUNCTION_2ARG(compute_leptonic_correlation, spinspin*,hadr, spinspin*,lept)
{
  GET_THREAD_ID();
  
  vector_reset(loc_corr);

    for(int ins=0;ins<nweak_ins;ins++)
      {
	//define a local storage
	spinspin thread_corr[glb_size[0]];
	memset(thread_corr,0,sizeof(spinspin)*glb_size[0]);
	
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  {
	    //multiply lepton side
	    spinspin l;
	    unsafe_spinspin_prod_dirac(l,lept[ivol],base_gamma+list_weak_insl[ins]);
	    
	    //trace hadron side
	    complex h;
	    trace_spinspin_with_dirac(h,hadr[ivol],base_gamma+list_weak_insq[ins]);
	    
	    //put in the local stack
	    int t=glb_coord_of_loclx[ivol][0];
	    spinspin_summ_the_complex_prod(thread_corr[t],l,h);
	  }
	THREAD_BARRIER();

	//reduce
	for(int ig=0;ig<16;ig++)
	  for(int t=0;t<glb_size[0];t++)
	    {
	      complex lh;
	      trace_spinspin_with_dirac(lh,thread_corr[t],base_gamma+ig);
	      glb_reduce_complex(glb_weak_corr[t+glb_size[0]*(ig+16*ins)],lh);
	    }
	if(IS_MASTER_THREAD) nlept_contr_tot+=16;
      }
}
THREADABLE_FUNCTION_END

//print a whole weak correlation function
void print_weak_correlations(FILE *fout)
{
  for(int ins=0;ins<nweak_ins;ins++)
    for(int ig=0;ig<16;ig++)
      {
	master_fprintf(fout," # qins=%d lins=%d proj=%d\n\n",list_weak_insq[ins],list_weak_insl[ins],ig);
	for(int t=source_coord[0];t<source_coord[0]+glb_size[0];t++)
	  {
	    int i=(t%glb_size[0])+glb_size[0]*(ig+16*ins);
	    master_fprintf(fout,"%+016.16lg %+016.16lg\n",glb_weak_corr[i][RE],glb_weak_corr[i][IM]);
	  }
	master_fprintf(fout,"\n");
      }
}

//compute the total hadroleptonic correlation functions
void compute_leptonic_correlation()
{
  master_printf("Computing leptonic correlation functions\n");
  lept_contr_time-=take_time();
  
  FILE *fout=open_file(combine("%s/hl_corr",outfolder).c_str(),"w");
  spinspin *hadr=nissa_malloc("hadr",loc_vol,spinspin);

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
	      int ip1=iqprop(iq1,(qins==0)?((phi_eta==0)?PROP_A:PROP_B):PROP_0,r2); //q1 will be reverted
	      int ip2=iqprop(iq2,(qins==1)?((phi_eta==0)?PROP_A:PROP_B):PROP_0,r2);
	      
	      //compute the hadronic part
	      hadronic_part_leptonic_correlation(hadr,Q[ip1],Q[ip2]);
	      
	      for(int orie=0;orie<2;orie++)
		for(int rl=0;rl<nr;rl++)
		  {
		    //contract with lepton
		    int il=ilprop(ilepton,orie,!phi_eta,rl);
		    compute_leptonic_correlation(hadr,L[il]);

		    //write down
		    master_fprintf(fout," # mq1=%lg mq2=%lg qins=%d qrev=%d ins=%s rq1=%d rq2=%d lep_orie=%+d\n\n",
			    qmass[iq1],qmass[iq2],qins+1,irev+1,(phi_eta==0)?"phi":"eta",!r2,r2,lepton_mom_sign[orie]);
		    print_weak_correlations(fout);
		  }
	    }
  
  close_file(fout);
  lept_contr_time+=take_time();
  nissa_free(hadr);
}

//wrapper to compute a single correlation function
void compute_hadronic_correlation(complex *corr,int *ig_so,PROP_TYPE *s1,int *ig_si,PROP_TYPE *s2,int ncontr)
{
  hadr_contr_time-=take_time();
  meson_two_points_Wilson_prop(glb_corr,loc_corr,ig_so,s1,ig_si,s2,ncontr);
  hadr_contr_time+=take_time();
  nhadr_contr_tot++;
}

//compute all the hadronic correlations
void compute_hadronic_correlation()
{
  master_printf("Computing hadronic correlation functions\n");
  
  const int ncombo_corr=9;
  qprop_t prop1_map[ncombo_corr]={PROP_0,PROP_0,PROP_0,PROP_0, PROP_0,PROP_A};
  qprop_t prop2_map[ncombo_corr]={PROP_0,PROP_S,PROP_P,PROP_AB,PROP_T,PROP_B};
  
  for(int icombo=0;icombo<ncombo_corr;icombo++)
    {
      FILE *fout=open_file(combine("%s/corr_%c%c",outfolder,prop_abbr[prop1_map[icombo]],prop_abbr[prop2_map[icombo]]).c_str(),"w");
      
      for(int imass=0;imass<nqmass;imass++)
	for(int jmass=0;jmass<nqmass;jmass++)
	  for(int r=0;r<nr;r++)
	    {
	      //compute the correlation function
	      int ig_so=5,ig_si=5;
	      compute_hadronic_correlation(glb_corr,&ig_so,Q[iqprop(imass,prop1_map[icombo],r)],&ig_si,Q[iqprop(jmass,prop2_map[icombo],r)],1);
	      
	      //print out
	      master_fprintf(fout," # m1(rev)=%lg m2(ins)=%lg r=%d\n",qmass[imass],qmass[jmass],r);
	      print_contractions_to_file(fout,1,&ig_so,&ig_si,glb_corr,source_coord[0],"",1.0);
	      master_fprintf(fout,"\n");		  
	    }
      
      //close the file
      close_file(fout);
    }
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
      generate_original_source();
      generate_photon_stochastic_propagator();
      
      generate_lepton_propagators();
      generate_quark_propagators();
      
      compute_leptonic_correlation();
      compute_hadronic_correlation();
            
      //pass to the next conf if there is enough time
      char fin_file[1024],run_file[1024];
      sprintf(fin_file,"%s/finished",outfolder);
      sprintf(run_file,"%s/running",outfolder);
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
