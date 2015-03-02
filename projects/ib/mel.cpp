#include <nissa.hpp>

#ifdef POINT_SOURCE_VERSION
 #define PROP_TYPE su3spinspin
#else
 #define PROP_TYPE colorspinspin
#endif

using namespace nissa;

/////////////////////////////////////// data //////////////////////////////

int ninv_tot=0;
double inv_time=0,tot_prog_time=0;

int wall_time;
int ngauge_conf,nanalyzed_conf=0;
char conf_path[1024],outfolder[1024];
quad_su3 *conf;

double kappa;
double put_theta[4],old_theta[4]={0,0,0,0};

coords source_coord;
PROP_TYPE *source,*original_source;
int seed,noise_type;

int nmass;
double *mass,*residue;
PROP_TYPE **S;

gauge_info photon;
double tadpole;
spin1field *photon_phi,*photon_eta;

complex *corr,*loc_corr;

//return appropriate propagator
enum insertion_t{ORIGINAL,SCALAR,PSEUDO,TIME_VECTOR_CURRENT,CONSERVED_CURRENT,STOCH_A,STOCH_B,TADPOLE};
enum prop_t{PROP_SIMPLE,PROP_SCALAR,PROP_PSEUDO,PROP_A,PROP_B,PROP_AB,PROP_T};
int nprop;
int iprop(int imass,prop_t ip,int r)
{return r+2*(imass+nmass*ip);}

//generate a wall-source for stochastic QCD propagator
void generate_original_source()
{
#ifdef POINT_SOURCE_VERSION
  generate_delta_source(original_source,source_coord);
#else
  generate_spindiluted_source(original_source,rnd_type_map[noise_type],source_coord[0]);
#endif
}

//compute the tadpole by summing all momenta
void compute_tadpole()
{
  double loc_tadpole=0;
  NISSA_LOC_VOL_LOOP(imom)
  {
    spin1prop prop;
    mom_space_tlSym_gauge_propagator_of_imom(prop,photon,imom);
    loc_tadpole+=prop[0][0][RE];
  }
  MPI_Allreduce(&loc_tadpole,&tadpole,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

//generate a stochastic QED propagator

//insert the conserved current
void insert_operator(PROP_TYPE *out,spin1field *curr,PROP_TYPE *in,complex fact_fw,complex fact_bw)
{
  //reset the output and communicate borders
  vector_reset(out);
  NAME3(communicate_lx,PROP_TYPE,borders)(in);
  communicate_lx_quad_su3_borders(conf);
  communicate_lx_spin1field_borders(curr);

  NISSA_LOC_VOL_LOOP(ivol)
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
	NAME2(PROP_TYPE,prodassign_complex)(fw,curr[ifw][mu]);
	NAME2(PROP_TYPE,prodassign_complex_conj)(bw,curr[ibw][mu]);
	
	//summ the two from the output
	NAME2(PROP_TYPE,summassign)(out[ivol],fw);
	NAME2(PROP_TYPE,summassign)(out[ivol],bw);
	
	//fw comes with -g, bw with +g
	PROP_TYPE temp1,temp2;
	NAME2(PROP_TYPE,subt)(temp1,bw,fw);
	NAME2(unsafe_dirac_prod,PROP_TYPE)(temp2,base_gamma+map_mu[mu],temp1);
	NAME2(PROP_TYPE,summassign)(out[ivol],temp2);
    }

  set_borders_invalid(out);
}

//insert the conserved vector current
void insert_conserved_vector_current(PROP_TYPE *out,PROP_TYPE *in)
{
  //prepare current insertion
  spin1field *aux=nissa_malloc("aux",loc_vol+bord_vol,spin1field);
  NISSA_LOC_VOL_LOOP(ivol) for(int mu=0;mu<4;mu++) complex_put_to_real(aux[ivol][mu],1);
  set_borders_invalid(aux);
  
  //call with no source insertion, minus between fw and bw, and a global -i*0.5
  complex fw_factor={0,-0.5},bw_factor={0,+0.5};
  insert_operator(out,aux,in,fw_factor,bw_factor);
  
  nissa_free(aux);
}

//insert the tadpole
void insert_tadpole(PROP_TYPE *out,PROP_TYPE *in)
{
  //prepare current insertion
  spin1field *aux=nissa_malloc("aux",loc_vol+bord_vol,spin1field);
  NISSA_LOC_VOL_LOOP(ivol) for(int mu=0;mu<4;mu++) complex_put_to_real(aux[ivol][mu],tadpole);
  set_borders_invalid(aux);
  
  //call with no source insertion, plus between fw and bw, and a global -0.25
  complex fw_factor={0,-0.25},bw_factor={0,-0.25};
  insert_operator(out,aux,in,fw_factor,bw_factor);
  
  nissa_free(aux);
}

//insert the external source, that is one of the two extrema of the stoch prop
void insert_external_source(PROP_TYPE *out,spin1field *curr,bool dag,PROP_TYPE *in)
{
  //prepare current insertion
  spin1field *aux=nissa_malloc("aux",loc_vol+bord_vol,spin1field);
  void (*fix_curr)(complex,complex)=(dag?complex_conj:complex_copy);
   NISSA_LOC_VOL_LOOP(ivol) for(int mu=0;mu<4;mu++) fix_curr(aux[ivol][mu],curr[ivol][mu]);
  set_borders_invalid(aux);
  
  //call with no source insertion, minus between fw and bw, and a global -i*0.5
  complex fw_factor={0,-0.5},bw_factor={0,+0.5};
  insert_operator(out,aux,in,fw_factor,bw_factor);
  
  nissa_free(aux);
}

//multiply with gamma
void prop_multiply_with_gamma(PROP_TYPE *out,int ig,PROP_TYPE *in)
{
  NISSA_LOC_VOL_LOOP(ivol)
    NAME2(safe_dirac_prod,PROP_TYPE)(out[ivol],base_gamma+ig,in[ivol]);
  set_borders_invalid(out);
}

//multiply with an imaginary factor
void prop_multiply_with_idouble(PROP_TYPE *out,double f)
{
  NISSA_LOC_VOL_LOOP(ivol)
    NAME2(PROP_TYPE,prodassign_idouble)(out[ivol],f);
  set_borders_invalid(out);
}

//generate a sequential source
void generate_source(int imass,int r,insertion_t inser,prop_t ORI=PROP_SIMPLE)
{
  PROP_TYPE *ori=S[iprop(imass,ORI,r)];
  switch(inser)
    {
    case ORIGINAL:prop_multiply_with_gamma(source,0,original_source);break;
    case SCALAR:prop_multiply_with_gamma(source,0,ori);break;
    case PSEUDO:prop_multiply_with_gamma(source,5,ori);break;
    case CONSERVED_CURRENT:insert_conserved_vector_current(source,ori);break;
    case TIME_VECTOR_CURRENT:prop_multiply_with_gamma(source,4,ori);break;
    case STOCH_A:insert_external_source(source,photon_eta,true,ori);break;
    case STOCH_B:insert_external_source(source,photon_phi,false,ori);break;
    case TADPOLE:insert_tadpole(source,ori);break;
    }
}

//invert on top of a source, putting all needed for the appropriate quark
void get_prop(PROP_TYPE *out,PROP_TYPE *in,int imass,bool r)
{
  //sign and rotators
  double sign_r[2]={-1,+1};
  dirac_matr *rot[2]={&Pminus,&Pplus};
  
  //allocate temporary vectors
  spincolor *temp_source=nissa_malloc("temp_source",loc_vol,spincolor);
  spincolor *temp_solution=nissa_malloc("temp_solution",loc_vol,spincolor);
  
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
	
	//rotate the source index
	safe_dirac_prod_spincolor(temp_source,rot[r],temp_source);
	
	//invert
	inv_time-=take_time();
	inv_tmD_cg_eoprec_eos(temp_solution,NULL,conf,kappa,sign_r[r]*mass[imass],100000,residue[imass],temp_source);
	ninv_tot++;inv_time+=take_time();
	
	//rotate the sink index
	safe_dirac_prod_spincolor(temp_solution,rot[r],temp_solution);      
	
	//put the output on place
#ifdef POINT_SOURCE_VERSION
	master_printf("  finished the inversion dirac index %d, color %d\n",id,ic);
	put_spincolor_into_su3spinspin(out,temp_solution,id,ic);
#else
	master_printf("  finished the inversion dirac index %d\n",id);	
	put_spincolor_into_colorspinspin(out,temp_solution,id);
#endif
      }
  
  //free temporary vectors
  nissa_free(temp_source);
  nissa_free(temp_solution);
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
  
  //Number of configurations
  read_str_int("NGaugeConf",&ngauge_conf);
  
  ///////////////////// finihed reading apart from conf list ///////////////
  
  //compute the tadpole summing all momentum
  compute_tadpole();
  
  //Allocate
  nprop=2*nmass*7;
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
  read_ildg_gauge_conf(conf,conf_path);
  master_printf("plaq: %.16g\n",global_plaquette_lx_conf(conf));
  
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
#ifdef BENCH
  master_printf(", of which:\n");
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf("  of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
#endif

  nissa_free(photon_eta);
  nissa_free(photon_phi);
  nissa_free(source);
  nissa_free(original_source);
  for(int iprop=0;iprop<nprop;iprop++) nissa_free(S[iprop]);
  nissa_free(S);
  nissa_free(conf);
  nissa_free(corr);
  nissa_free(loc_corr);
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
{meson_two_points_Wilson_prop(corr,loc_corr,&ig_so,s1,&ig_si,s2,1);}

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
      generate_stochastic_tlSym_gauge_propagator(photon_phi,photon_eta,photon);
      
      //map the source, the destination and the insertion for each propagator
      const insertion_t insertion_map[7]={ORIGINAL,SCALAR,PSEUDO,STOCH_A,STOCH_B,STOCH_B,TADPOLE};
      const prop_t prop_map[7]={PROP_SIMPLE,PROP_SCALAR,PROP_PSEUDO,PROP_A,PROP_B,PROP_AB,PROP_T};
      const prop_t source_map[7]={PROP_SIMPLE,PROP_SIMPLE,PROP_SIMPLE,PROP_SIMPLE,PROP_SIMPLE,PROP_A,PROP_SIMPLE};
      const char prop_name[7][20]={"SIMPLE","SCALAR","PSEUDO","STOCH_A","STOCH_B","STOCH_AB","TADPOLE"};
      const char prop_abbr[]="0SPABXT";
      
      //generate all the propagators	
      for(int ip=0;ip<7;ip++)
	{
	  master_printf("Generating propagtor of type %s using as a source %s\n",prop_name[ip],prop_name[source_map[ip]]);
	  for(int imass=0;imass<nmass;imass++)
	    for(int r=0;r<2;r++)
	      {
		master_printf(" mass[%d]=%lg, r=%d\n",imass,mass[imass],r);
		generate_source(imass,r,insertion_map[ip],source_map[ip]);
		get_prop(S[iprop(imass,prop_map[ip],r)],source,imass,r);
	      }
	}
      
      //compute all correlations
      const int ncombo_corr=6;
      prop_t prop1_map[ncombo_corr]={PROP_SIMPLE,PROP_SIMPLE,PROP_SIMPLE,PROP_SIMPLE,PROP_SIMPLE,PROP_A};
      prop_t prop2_map[ncombo_corr]={PROP_SIMPLE,PROP_SCALAR,PROP_PSEUDO,PROP_AB,PROP_T,PROP_B};
      for(int icombo=0;icombo<ncombo_corr;icombo++)
	{
	  FILE *fout=open_file(combine("%s/corr_%c%c",outfolder,prop_abbr[prop1_map[icombo]],prop_abbr[prop2_map[icombo]]).c_str(),"w");
	  
	  for(int imass=0;imass<nmass;imass++)
	    for(int jmass=0;jmass<nmass;jmass++)
	      for(int r=0;r<2;r++)
		{
		  //compute the correlation function
		  compute_corr(corr,S[iprop(imass,prop1_map[icombo],r)],S[iprop(imass,prop2_map[icombo],r)],5,5);
		  
		  //print out
		  master_fprintf(fout," # m1=%lg m2=%lg r=%d\n\n",mass[imass],mass[jmass],r);
		  for(int t=source_coord[0];t<glb_size[0]+source_coord[0];t++)
		    master_fprintf(fout,"%d %+015.15lg\t%+015.15lg\n",t-source_coord[0],corr[t%glb_size[0]][RE],corr[t%glb_size[0]][IM]);
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
