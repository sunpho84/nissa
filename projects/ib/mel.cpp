#include <nissa.hpp>

#ifdef POINT_SOURCE_VERSION
 #define PROP_TYPE su3spinspin
#else
 #define PROP_TYPE colorspinspin
#endif

using namespace nissa;

/////////////////////////////////////// data //////////////////////////////

int ninv_tot=0,ncontr_tot=0,nsource_tot=0,nphoton_prop_tot=0;
double inv_time=0,contr_time=0,tot_prog_time=0,source_time=0,photon_prop_time=0;

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

complex *corr,*loc_corr;

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
double *lep_mass;
double *lep_energy;
double *lep_mom;
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
void get_prop(PROP_TYPE *out,PROP_TYPE *in,int imass,bool r,int rotate=true)
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

//compute the energy of a particle according to lattice dispertion relation
double lep_energy_fun(double m,double p,double k=0.125)
{
  /*
  double sinp=sin(p);
  double sinph=sin(p/2);
  
  double M=1/(2*k)-4;

  double Mp=0,sin2_mom=0;
  for(int mu=0;mu<NDIM;mu++)
    {
      Mp+=sinph*sinph;
      sin2_mom+=sinp*sinp;
    }
  Mp=2*Mp+M;
    
  double den=sin2_mom+Mp*Mp+m*m;

  return den;
  */

  return sqrt(m*m+3*p*p);
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
  lep_mass=nissa_malloc("lep_mass",nleptons,double);
  lep_energy=nissa_malloc("lep_energy",nleptons,double);
  lep_mom=nissa_malloc("lep_mom",nleptons,double);
  expect_str("Q1Q2LMesMassTheta");
  for(int il=0;il<nleptons;il++)
    {
      read_int(lep_corr_iq1+il);
      read_int(lep_corr_iq2+il);
      read_double(lep_mass+il);
      double mes_mass;
      read_double(&mes_mass);

      //compute meson momentum and energy
      double lep_mom_mod=(pow(mes_mass,2)-pow(lep_mass[il],2))/(2*mes_mass);
      lep_mom[il]=lep_mom_mod/sqrt(3);
      lep_energy[il]=lep_energy_fun(lep_mass[il],lep_mom[il]);
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
  nqprop=nr*nqmass*nqprop_kind;
  nlprop=2*2*nr*nleptons;
  
  //allocate temporary vectors
  temp_source=nissa_malloc("temp_source",loc_vol,spincolor);
  temp_solution=nissa_malloc("temp_solution",loc_vol,spincolor);
  corr=nissa_malloc("corr",glb_size[0],complex);
  loc_corr=nissa_malloc("loc_corr",glb_size[0],complex);
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
	    get_prop(Q[iqprop(imass,prop_map[ip],r)],source,imass,r);
	  }
    }
}

//compute the phase for lepton
void get_lepton_phase_factor(complex out,int ivol,int ilepton,quark_info le)
{
  //fetch the energy and momentum ot the muon
  double E=lep_energy[ilepton];
  double p=lep_mom[ilepton];
  
  //compute pure time factor and summ the coordinate
  double ext=exp(glb_coord_of_loclx[ivol][0]*E);
  int summ_coord=0;
  for(int i=1;i<NDIM;i++) summ_coord+=glb_coord_of_loclx[ivol][i];
  
  //compute full exponential and multiply each component of the photon field
  out[RE]=cos(-summ_coord*p)*ext;
  out[IM]=sin(-summ_coord*p)*ext;
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

  //for the time being the boundaries are periodic
  quark_info le;
  le.kappa=0.125;
  memset(le.bc,0,sizeof(momentum_t));
  
  //compute the 4 pairs of gammas: i(-i t3 g5-+gmu)/2=(t3 g5-+i gmu)/2=
  spinspin g[2][2][4]; //[t3][-+][mu]
  double mps[2]={-0.5,+0.5};
  memset(g,0,2*2*4*sizeof(spinspin));
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
	    le.mass=lep_energy[ilepton];
	    
	    //prepare each propagator for a single lepton, by computing
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
		  get_lepton_phase_factor(ph0,ibw,ilepton,le);
		  get_lepton_phase_factor(ph1,ifw,ilepton,le);
		  
		  //multiply by A
		  complex phA0,phA1;
		  unsafe_complex_prod(phA0,ph0,A[ibw][mu]);
		  unsafe_complex_prod(phA1,ph0,A[ivol][mu]);
		  
		  //summ the two pieces
		  spinspin_summ_the_complex_prod(L[iprop][ivol],g[r][0][mu],phA0);
		  spinspin_subt_the_complex_prod(L[iprop][ivol],g[r][1][mu],phA1);
		}
	    set_borders_invalid(L[iprop]);
	    
	    //multiply by the lepton propagator (from the right)
	    multiply_from_right_by_x_space_twisted_propagator_by_fft(L[iprop],L[iprop],le);
	  }
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
  master_printf("Total time: %g",tot_prog_time);
  master_printf(", of which:\n");
  master_printf(" - %02.2f%s to prepare %d photon stochastic propagators (%2.2gs avg)\n",photon_prop_time/tot_prog_time*100,"%",nphoton_prop_tot,photon_prop_time/nphoton_prop_tot);
  master_printf(" - %02.2f%s to prepare %d generalized sources (%2.2gs avg)\n",source_time/tot_prog_time*100,"%",nsource_tot,source_time/nsource_tot);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf("  of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d contractions (%2.2gs avg)\n",contr_time/tot_prog_time*100,"%",ncontr_tot,contr_time/ncontr_tot);

  nissa_free(photon_eta);
  nissa_free(photon_phi);
  nissa_free(source);
  nissa_free(original_source);
  for(int iprop=0;iprop<nqprop;iprop++) nissa_free(Q[iprop]);
  nissa_free(Q);
  nissa_free(conf);
  nissa_free(corr);
  nissa_free(loc_corr);
  nissa_free(temp_source);
  nissa_free(temp_solution);
  nissa_free(lep_corr_iq1);
  nissa_free(lep_corr_iq2);
  nissa_free(lep_mass);
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

//wrapper to compute a single correlation function
void compute_corr(complex *corr,int *ig_so,PROP_TYPE *s1,int *ig_si,PROP_TYPE *s2,int ncontr)
{
  contr_time-=take_time();
  meson_two_points_Wilson_prop(corr,loc_corr,ig_so,s1,ig_si,s2,ncontr);
  contr_time+=take_time();
  ncontr_tot++;
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
      
      generate_quark_propagators();
      generate_lepton_propagators();
      
      //compute all correlations
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
		  const int ncontr=3;
		  int ig_so[ncontr]={5,4,5},ig_si[ncontr]={5,5,4};
		  compute_corr(corr,ig_so,Q[iqprop(imass,prop1_map[icombo],r)],ig_si,Q[iqprop(jmass,prop2_map[icombo],r)],ncontr);
		  
		  //print out
		  master_fprintf(fout," # m1(rev)=%lg m2(ins)=%lg r=%d\n",qmass[imass],qmass[jmass],r);
		  print_contractions_to_file(fout,ncontr,ig_so,ig_si,corr,source_coord[0],"",1.0);
		  master_fprintf(fout,"\n");		  
		}
	  
	  //close the file
	  close_file(fout);
	}
      
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
