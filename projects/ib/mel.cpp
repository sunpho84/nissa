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

int nmass,nr;
double *mass,*residue;
PROP_TYPE **S;

spincolor *temp_source;
spincolor *temp_solution;
  
gauge_info photon;
double tadpole[4];
spin1field *photon_phi,*photon_eta;

complex *corr,*loc_corr;

const int V_curr_mu=0;
const double tau3[2]={-1,+1};


//#define MEL_DEBUG

#ifdef MEL_DEBUG
//define types of propagator used
const int nprop_kind=3;
enum insertion_t                     {ORIGINAL,  CONS_VEC0,  WALL5 };
const char ins_name[nprop_kind][20]={"ORIGINAL","CONS_VEC0","WALL5"};
enum prop_t{                           PROP_0,  PROP_CV,  PROP_W};
const char prop_name[nprop_kind][20]={"PROP_0","PROP_CV","PROP_W"};
//map the source, the destination and the insertion for each propagator
const prop_t prop_map[nprop_kind]=          {PROP_0,   PROP_CV,   PROP_W };
const insertion_t insertion_map[nprop_kind]={ORIGINAL, CONS_VEC0, WALL5 };
const prop_t source_map[nprop_kind]=        {PROP_0,   PROP_0,    PROP_0};
const char prop_abbr[]=                      "0"       "C"        "W";
#else
//define types of propagator used
const int nins_kind=9;
enum insertion_t{                    ORIGINAL,  SCALAR,  PSEUDO,  VEC0,  CONS_VEC0,  STOCH_PHI,  STOCH_ETA,  TADPOLE,  WALL5 };
const char ins_name[nins_kind][20]={"ORIGINAL","SCALAR","PSEUDO","VEC0","CONS_VEC0","STOCH_PHI","STOCH_ETA","TADPOLE","WALL5"};
const int nprop_kind=10;
enum prop_t{                           PROP_0,  PROP_S,  PROP_P,  PROP_A,  PROP_B,  PROP_AB,  PROP_T,  PROP_CV, PROP_V,  PROP_W};
const char prop_name[nprop_kind][20]={"PROP_0","PROP_S","PROP_P","PROP_A","PROP_B","PROP_AB","PROP_T","PROP_CV","PROP_V","PROP_W"};
//map the source, the destination and the insertion for each propagator
const prop_t prop_map[nprop_kind]=          {PROP_0,   PROP_S, PROP_P, PROP_A,    PROP_B,    PROP_AB,   PROP_T,  PROP_CV,   PROP_V, PROP_W };
const insertion_t insertion_map[nprop_kind]={ORIGINAL, SCALAR, PSEUDO, STOCH_PHI, STOCH_ETA, STOCH_ETA, TADPOLE, CONS_VEC0, VEC0,   WALL5 };
const prop_t source_map[nprop_kind]=        {PROP_0,   PROP_0, PROP_0, PROP_0,    PROP_0,    PROP_A,    PROP_0,  PROP_0,    PROP_0, PROP_0};
const char prop_abbr[]=                      "0"       "S"     "P"     "A"        "B"        "X"        "T"      "C"        "V"     "W";
#endif
   
//return appropriate propagator
int nprop;
int iprop(int imass,prop_t ip,int r)
{return r+nr*(imass+nmass*ip);}

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

//compute the tadpole by taking the zero momentum ft of momentum prop
void compute_tadpole()
{
  double tad_time=-take_time();
  spin1prop *gprop=nissa_malloc("gprop",loc_vol,spin1prop);
  compute_x_space_tlSym_gauge_propagator_by_fft(gprop,photon);
  for(int mu=0;mu<4;mu++) tadpole[mu]=broadcast(gprop[0][mu][mu][RE]);
  nissa_free(gprop);
  tad_time+=take_time();
  master_printf("Tadpole: (%lg,%lg,%lg,%lg), time to compute: %lg s\n",tadpole[0],tadpole[1],tadpole[2],tadpole[3],tad_time);
}

/*
  insert the operator:  \sum_mu A_{x,mu} [
   f_fw * ( - i t3 g5 - gmu) U_{x,mu} \delta{x',x+mu} + f_bw * ( - i t3 g5 + gmu) U^+_{x-mu,mu} \delta{x',x-mu}]
*/
void insert_operator(PROP_TYPE *out,spin1field *curr,PROP_TYPE *in,complex fact_fw,complex fact_bw,int r,void(*get_curr)(complex,spin1field *,int,int))
{
  GET_THREAD_ID();
  
  //reset the output and communicate borders
  vector_reset(out);
  NAME3(communicate_lx,PROP_TYPE,borders)(in);
  communicate_lx_quad_su3_borders(conf);
  if(curr) communicate_lx_spin1field_borders(curr);

  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      {
	//find neighbors
	int ifw=loclx_neighup[ivol][mu];
	int ibw=loclx_neighdw[ivol][mu];
	
	//transport down and up
	PROP_TYPE fw,bw;
	NAME2(unsafe_su3_prod,PROP_TYPE)(fw,conf[ivol][mu],in[ifw]);
	NAME2(unsafe_su3_dag_prod,PROP_TYPE)(bw,conf[ibw][mu],in[ibw]);
	
	//weight the two pieces
	NAME2(PROP_TYPE,prodassign_complex)(fw,fact_fw);
	NAME2(PROP_TYPE,prodassign_complex)(bw,fact_bw);
	
	//insert the current
	complex fw_curr,bw_curr;
	get_curr(fw_curr,curr,ivol,mu);
	get_curr(bw_curr,curr,ibw,mu);
	NAME2(PROP_TYPE,prodassign_complex)(fw,fw_curr);
	NAME2(PROP_TYPE,prodassign_complex)(bw,bw_curr);
	
	//summ and subtract the two
	PROP_TYPE bw_M_fw,bw_P_fw;
	NAME2(PROP_TYPE,subt)(bw_M_fw,bw,fw);
	NAME2(PROP_TYPE,summ)(bw_P_fw,bw,fw);
	
	//put -i g5 t3 on the summ
	PROP_TYPE g5_bw_P_fw;
	NAME2(unsafe_dirac_prod,PROP_TYPE)(g5_bw_P_fw,base_gamma+5,bw_P_fw);
	NAME2(PROP_TYPE,summ_the_prod_idouble)(out[ivol],g5_bw_P_fw,-tau3[r]);
	
	//put gmu on the difference
	PROP_TYPE gmu_bw_M_fw;
	NAME2(unsafe_dirac_prod,PROP_TYPE)(gmu_bw_M_fw,base_gamma+map_mu[mu],bw_M_fw);
	NAME2(PROP_TYPE,summassign)(out[ivol],gmu_bw_M_fw);
    }

  set_borders_invalid(out);
}

/*
  insert the operator:  \sum_mu A_{x,mu} [
   f_fw * ( 1 - gmu) U_{x,mu} \delta{x',x+mu} + f_bw * ( 1 + gmu) U^+_{x-mu,mu} \delta{x',x-mu}]
*/
void insert_operator_twisted(PROP_TYPE *out,spin1field *curr,PROP_TYPE *in,complex fact_fw,complex fact_bw,void(*get_curr)(complex,spin1field *,int,int))
{
  GET_THREAD_ID();
  
  //reset the output and communicate borders
  vector_reset(out);
  NAME3(communicate_lx,PROP_TYPE,borders)(in);
  communicate_lx_quad_su3_borders(conf);
  if(curr) communicate_lx_spin1field_borders(curr);

  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      {
	//find neighbors
	int ifw=loclx_neighup[ivol][mu];
	int ibw=loclx_neighdw[ivol][mu];
	
	//transport down and up
	PROP_TYPE fw,bw;
	NAME2(unsafe_su3_prod,PROP_TYPE)(fw,conf[ivol][mu],in[ifw]);
	NAME2(unsafe_su3_dag_prod,PROP_TYPE)(bw,conf[ibw][mu],in[ibw]);
	
	//weight the two pieces
	NAME2(PROP_TYPE,prodassign_complex)(fw,fact_fw);
	NAME2(PROP_TYPE,prodassign_complex)(bw,fact_bw);
	
	//insert the current
	complex fw_curr,bw_curr;
	get_curr(fw_curr,curr,ivol,mu);
	get_curr(bw_curr,curr,ibw,mu);
	NAME2(PROP_TYPE,prodassign_complex)(fw,fw_curr);
	NAME2(PROP_TYPE,prodassign_complex)(bw,bw_curr);
	
	//summ and subtract the two
	PROP_TYPE bw_M_fw,bw_P_fw;
	NAME2(PROP_TYPE,subt)(bw_M_fw,bw,fw);
	NAME2(PROP_TYPE,summ)(bw_P_fw,bw,fw);
	
	//add the summ
	NAME2(PROP_TYPE,summassign)(out[ivol],bw_P_fw);

	//put gmu on the difference
	PROP_TYPE gmu_bw_M_fw;
	NAME2(unsafe_dirac_prod,PROP_TYPE)(gmu_bw_M_fw,base_gamma+map_mu[mu],bw_M_fw);
	NAME2(PROP_TYPE,summassign)(out[ivol],gmu_bw_M_fw);
    }

  set_borders_invalid(out);
}

//insert the conserved vector current
void insert_conserved_vector_current_handle(complex out,spin1field *aux,int ivol,int mu)
{out[RE]=(mu==V_curr_mu);out[IM]=0;}
THREADABLE_FUNCTION_3ARG(insert_conserved_vector_current, PROP_TYPE*,out, PROP_TYPE*,in, int,r)
{
  //call with no source insertion, minus between fw and bw, and a global -0.25
  complex fw_factor={-0.25,0},bw_factor={+0.25,0};
  insert_operator(out,NULL,in,fw_factor,bw_factor,r,insert_conserved_vector_current_handle);
}
THREADABLE_FUNCTION_END

//insert the tadpole
void insert_tadpole_handle(complex out,spin1field *aux,int ivol,int mu){out[RE]=tadpole[mu];out[IM]=0;}
THREADABLE_FUNCTION_3ARG(insert_tadpole, PROP_TYPE*,out, PROP_TYPE*,in, int,r)
{
  //call with no source insertion, plus between fw and bw, and a global -0.25
  complex fw_factor={-0.25,0},bw_factor={-0.25,0};
  insert_operator(out,NULL,in,fw_factor,bw_factor,r,insert_tadpole_handle);
}
THREADABLE_FUNCTION_END

//insert the external source, that is one of the two extrema of the stoch prop
void insert_external_source_handle(complex out,spin1field *aux,int ivol,int mu){complex_copy(out,aux[ivol][mu]);}
THREADABLE_FUNCTION_4ARG(insert_external_source, PROP_TYPE*,out, spin1field*,curr, PROP_TYPE*,in, int,r)
{
  //call with no source insertion, minus between fw and bw, and a global -i*0.5
  complex fw_factor={0,-0.5},bw_factor={0,+0.5};
  insert_operator(out,curr,in,fw_factor,bw_factor,r,insert_external_source_handle);
}
THREADABLE_FUNCTION_END

//multiply with gamma
THREADABLE_FUNCTION_4ARG(prop_multiply_with_gamma, PROP_TYPE*,out, int,ig, PROP_TYPE*,in, int,twall)
{
  GET_THREAD_ID();
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      NAME2(safe_dirac_prod,PROP_TYPE)(out[ivol],base_gamma+ig,in[ivol]);
      NAME2(PROP_TYPE,prodassign_double)(out[ivol],(twall==-1||glb_coord_of_loclx[ivol][0]==twall));
    }
  set_borders_invalid(out);
}
THREADABLE_FUNCTION_END
void prop_multiply_with_gamma(PROP_TYPE *out,int ig,PROP_TYPE *in)
{prop_multiply_with_gamma(out,ig,in,-1);}

//multiply with an imaginary factor
THREADABLE_FUNCTION_2ARG(prop_multiply_with_idouble, PROP_TYPE*,out, double,f)
{
  GET_THREAD_ID();
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    NAME2(PROP_TYPE,prodassign_idouble)(out[ivol],f);
  set_borders_invalid(out);
}
THREADABLE_FUNCTION_END

//generate a sequential source
void generate_source(insertion_t inser,int r,PROP_TYPE *ori)
{
  source_time-=take_time();
  
  switch(inser)
    {
    case ORIGINAL:prop_multiply_with_gamma(source,0,original_source);break;
#ifndef MEL_DEBUG
    case SCALAR:prop_multiply_with_gamma(source,0,ori);break;
    case PSEUDO:prop_multiply_with_gamma(source,5,ori);break;
#endif
    case CONS_VEC0:insert_conserved_vector_current(source,ori,r);break;
#ifndef MEL_DEBUG
    case VEC0:prop_multiply_with_gamma(source,4,ori);break;
    case STOCH_PHI:insert_external_source(source,photon_phi,ori,r);break;
    case STOCH_ETA:insert_external_source(source,photon_eta,ori,r);break;
    case TADPOLE:insert_tadpole(source,ori,r);break;
#endif
    case WALL5:prop_multiply_with_gamma(source,5,ori,(glb_size[0]/2+source_coord[0])%glb_size[0]);break;
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
	inv_tmD_cg_eoprec_eos(temp_solution,NULL,conf,kappa,tau3[r]*mass[imass],100000,residue[imass],temp_source);
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
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L);
  //Wall time
  read_str_int("WallTime",&wall_time);
  //Kappa
  read_str_double("Kappa",&kappa);
  //One or two r
  read_str_int("NR",&nr);
  //Masses and residue
  read_list_of_double_pairs("MassResidues",&nmass,&mass,&residue);
  
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
  compute_tadpole();
  
  //Allocate
  nprop=2*nmass*nprop_kind;
  //allocate temporary vectors
  temp_source=nissa_malloc("temp_source",loc_vol,spincolor);
  temp_solution=nissa_malloc("temp_solution",loc_vol,spincolor);
  corr=nissa_malloc("corr",glb_size[0],complex);
  loc_corr=nissa_malloc("loc_corr",glb_size[0],complex);
  original_source=nissa_malloc("source",loc_vol,PROP_TYPE);
  source=nissa_malloc("source",loc_vol,PROP_TYPE);
  photon_eta=nissa_malloc("photon_eta",loc_vol+bord_vol,spin1field);
  photon_phi=nissa_malloc("photon_phi",loc_vol+bord_vol,spin1field);
  S=nissa_malloc("S*",nprop,PROP_TYPE*);
  for(int iprop=0;iprop<nprop;iprop++)
    S[iprop]=nissa_malloc("S",loc_vol+bord_vol,PROP_TYPE);
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
  for(int iprop=0;iprop<nprop;iprop++) nissa_free(S[iprop]);
  nissa_free(S);
  nissa_free(conf);
  nissa_free(corr);
  nissa_free(loc_corr);
  nissa_free(temp_source);
  nissa_free(temp_solution);
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
void compute_corr(complex *corr,PROP_TYPE *s1,PROP_TYPE *s2,int ig_so,int ig_si)
{
  contr_time-=take_time();
  meson_two_points_Wilson_prop(corr,loc_corr,&ig_so,s1,&ig_si,s2,1);
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
      //smear the conf and generate the source
      setup_conf();
      generate_original_source();
      generate_photon_stochastic_propagator();
      
      //generate all the propagators
      for(int ip=0;ip<nprop_kind;ip++)
	{
	  master_printf("Generating propagtor of type %s inserting %s on source %s\n",prop_name[prop_map[ip]],ins_name[insertion_map[ip]],prop_name[source_map[ip]]);
	  for(int imass=0;imass<nmass;imass++)
	    for(int r=0;r<nr;r++)
	      {
		master_printf(" mass[%d]=%lg, r=%d\n",imass,mass[imass],r);
		generate_source(insertion_map[ip],r,S[iprop((prop_map[ip]==PROP_W)?0:imass,source_map[ip],r)]);
		get_prop(S[iprop(imass,prop_map[ip],r)],source,imass,(prop_map[ip]==PROP_W)?(!r):r);
	      }
	}
      
      //compute all correlations
#ifdef MEL_DEBUG
      const int ncombo_corr=3;
      prop_t prop1_map[ncombo_corr]={PROP_0,PROP_0, PROP_0};
      prop_t prop2_map[ncombo_corr]={PROP_0,PROP_CV,PROP_W};
#else
      const int ncombo_corr=9;
      prop_t prop1_map[ncombo_corr]={PROP_0,PROP_0,PROP_0,PROP_0, PROP_0,PROP_A,PROP_0,PROP_0, PROP_0};
      prop_t prop2_map[ncombo_corr]={PROP_0,PROP_S,PROP_P,PROP_AB,PROP_T,PROP_B,PROP_V,PROP_CV,PROP_W};
#endif
      for(int icombo=0;icombo<ncombo_corr;icombo++)
	{
	  FILE *fout=open_file(combine("%s/corr_%c%c",outfolder,prop_abbr[prop1_map[icombo]],prop_abbr[prop2_map[icombo]]).c_str(),"w");
	  
	  for(int imass=0;imass<nmass;imass++)
	    for(int jmass=0;jmass<nmass;jmass++)
	      for(int r=0;r<nr;r++)
		{
		  //compute the correlation function
		  int gso=5,gsi=(prop2_map[icombo]==PROP_W)?map_mu[V_curr_mu]:5;
		  compute_corr(corr,S[iprop(imass,prop1_map[icombo],r)],S[iprop(jmass,prop2_map[icombo],r)],gso,gsi);
		  
		  //print out
		  master_fprintf(fout," # m1(rev)=%lg m2(ins)=%lg r=%d\n",mass[imass],mass[jmass],r);
		  print_contractions_to_file(fout,1,&gso,&gsi,corr,source_coord[0],"",1.0);
		  master_fprintf(fout,"\n");
		  
		  //consider also the conserved current
		  if(prop2_map[icombo]==PROP_W)
		    {
		      generate_source(CONS_VEC0,!r,S[iprop(imass,PROP_W,r)]);
		      
		      //compute the correlation function
		      int gso=5,gsi=0;
		      compute_corr(corr,S[iprop(imass,PROP_0,r)],source,gso,gsi);
		      
		      //print out
		      print_contractions_to_file(fout,1,&gso,&gsi,corr,source_coord[0],"C_",1.0);
		      master_fprintf(fout,"\n");
		    }
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
