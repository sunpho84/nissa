/*
  This program can be used to generate gauge configurations
  according to the rooted-staggered action in presence of 
  electromagnetic fields and/or imaginary chemical potentials.
  
  The molecular dynamic routines are in the file:
   ../../src/hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.cpp
*/

#include <math.h>

#include "nissa.hpp"

using namespace nissa;

//observables
char gauge_obs_path[1024];

//input and output path for confs
char conf_path[1024];
char store_conf_path[1024];

//conf and staples
quad_su3 *conf;
squared_staples_t *squared_staples;    
rectangular_staples_t *rectangular_staples;

//evol pars
double beta;
gauge_action_name_t gauge_action_name;
pure_gauge_evol_pars_t evol_pars;

//traj
int cube_size,prod_ntraj;
int itraj,max_ntraj;
int store_conf_each;
int store_running_temp_conf;
int conf_created;
int seed;

const int HOT=0,COLD=1;

//write a conf adding info
void write_conf()
{
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  char text[1024];

  //traj id
  sprintf(text,"%d",itraj);
  ILDG_string_message_append_to_last(&mess,"MD_traj",text);
  
  //skip 10 random numbers
  for(int iskip=0;iskip<10;iskip++) rnd_get_unif(&glb_rnd_gen,0,1);

  //glb_rnd_gen status
  convert_rnd_gen_to_text(text,&glb_rnd_gen);
  ILDG_string_message_append_to_last(&mess,"RND_gen_status",text);
  
  //write the conf
  write_ildg_gauge_conf(conf_path,conf,64,&mess);
  
  //free messages
  ILDG_message_free_all(&mess);
}

//read conf
void read_conf()
{
  //init messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);

  //read the conf
  read_ildg_gauge_conf(conf,conf_path,&mess);
  
  //scan messages
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    {  
      if(strcasecmp(cur_mess->name,"MD_traj")==0) sscanf(cur_mess->data,"%d",&itraj);
      if(strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
    }
  
  //if message with string not found start from input seed
  if(loc_rnd_gen_inited==0)
    {
      master_printf("RND_gen status not found inside conf, starting from input passed seed\n");
      start_loc_rnd_gen(seed);
    }
  
  ILDG_message_free_all(&mess);
}

//initialize the simulation
void init_simulation(char *path)
{
  //////////////////////////// read the input /////////////////////////
  
  //open input file
  open_input(path);
  
  //init the grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);  
  
  read_str_str("GaugeObsPath",gauge_obs_path,1024); //gauge observables path
  read_str_int("MaxNTraj",&max_ntraj); //number of trajectory to evolve
  read_str_int("Seed",&seed); //seed

  //kind of action and evolution pars
  char gauge_action_name_str[1024];
  read_str_str("GaugeAction",gauge_action_name_str,1024);
  if(strcmp(gauge_action_name_str,"Wilson")==0) gauge_action_name=Wilson_action;
  else
    if(strcmp(gauge_action_name_str,"tlSym")==0) gauge_action_name=tlSym_action;
    else crash("unknown gauge action: %s",gauge_action_name_str);
  read_str_double("Beta",&beta);
  read_pure_gauge_evol_pars(evol_pars);
  
  //read in and out conf path
  read_str_str("ConfPath",conf_path,1024);
  read_str_str("StoreConfPath",store_conf_path,1024);
  read_str_int("StoreConfEach",&store_conf_each);
  read_str_int("StoreRunningTempConf",&store_running_temp_conf);
  
  //read if configuration must be generated cold or hot
  char start_conf_cond_str[1024];
  read_str_str("StartConfCond",start_conf_cond_str,1024);
  int start_conf_cond=-1;
  if(strcasecmp(start_conf_cond_str,"HOT")==0) start_conf_cond=HOT;
  if(strcasecmp(start_conf_cond_str,"COLD")==0) start_conf_cond=COLD;
  if(start_conf_cond==-1) crash("unknown starting condition cond %s, expected 'HOT' or 'COLD'",start_conf_cond_str);
  
  close_input();
  
  ////////////////////////// allocate stuff ////////////////////////
   
  //allocate conf and staples
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  rectangular_staples=nissa_malloc("rectangular_staples",loc_vol+bord_vol,rectangular_staples_t);
  squared_staples=nissa_malloc("squared_staples",loc_vol+bord_vol,squared_staples_t);  

  //load conf or generate it
  if(file_exists(conf_path))
    {
      master_printf("File %s found, loading\n",conf_path);
      read_conf();
      conf_created=false;
    }
  else
    {
      conf_created=true;
      
      //start the random generator using passed seed
      start_loc_rnd_gen(seed);
      
      //generate hot or cold conf
      if(start_conf_cond==HOT)
	{
	  master_printf("File %s not found, generating hot conf\n",conf_path);
	  generate_hot_lx_conf(conf);
	}
      else
	{
	  master_printf("File %s not found, generating cold conf\n",conf_path);
	  generate_cold_lx_conf(conf);
	}
      
      //reset conf id
      itraj=0;
    }  
  
  /*
  //find cube size
  int max_cube_size=std::min(loc_size[0],std::min(loc_size[1],std::min(loc_size[2],loc_size[3])));
  int cube_fit=0,cube_numb[4],cube_size=3,
  while(cube_size<max_cube_size && !cube_fit)
    {
      cube_fit=1;
      for(int mu=0;mu<4;mu++)
	{
	  cube_numb[mu]=loc_size[mu]/cube_size;
	  cube_fit&=(cube_numb[mu]*cube_size==loc_size[mu]);
	}
      if(!cube_fit) cube_size++;
    }
  if(!cube_fit) crash("no possible cube partitioning of local size");
  */
}

//finalize everything
void close_simulation()
{
  if(!store_running_temp_conf) write_conf();
  nissa_free(conf);
  nissa_free(squared_staples);
  nissa_free(rectangular_staples);
}

//updated all sites
void sweep_this_rank(int ho)
{
/*
  int ho=(isweep<evol_pars.nhb_sweeps);
  
  //receive staple from neighbours
  
  for(int mu=0;mu<4;mu++)
    for(int cube_par=0;cube_par<cube_size;cube_par++)
      {
	NISSA_PARALLEL_LOOP(ired,0,loc_vol/cube_size)
	  {
	    int ivol=ired*cube_size+cube_par;
	    
	    //compute staples
	    su3 staples;
	    
	    //compute heatbath or overrelax link and change it
	    su3 link;
	    if(ho==0) su3_find_heatbath(link,conf[ivol][mu],staples,beta,evol_pars.nhb_hits,loc_rnd_gen+ivol);
	    else      su3_find_overrelaxed(link,conf[ivol][mu],staples,evol_pars.nov_hits);
	    su3_copy(conf[ivol][mu],link);
	  }
	THREAD_BARRIER();
      }
*/
}

//compute for other ranks
void compute_for_other_ranks()
{
}

//heatbath or overrelax algorithm for the quenched simulation case, Wilson action
void generate_new_conf()
{
  /*
  GET_THREAD_ID();
  
  if(nranks>1 && nranks%2!=0) crash("expected 1 ranks or nranks multiple of 2");
  
  for(int isweep=0;isweep<evol_pars.nhb_sweeps+evol_pars.nov_sweeps;isweep++)
    //only half of the ranks can be simultaneously updated
    for(int rank_par=0;rank_par<2;rank_par++)
      {
	if(rank_par==(rank%2)) sweep_this_rank();
	else compute_for_other_ranks();

	//sync ranks
	MPI_BARRIER(MPI_COMM_WORLD);
      }
  */
  set_borders_invalid(conf); 
}

//measure plaquette and polyakov loop
void measure_gauge_obs()
{
  //open creating or appending
  FILE *file=open_file(gauge_obs_path,conf_created?"w":"a");

  //paths
  double paths[2];
  memset(conf+loc_vol,0,sizeof(quad_su3)*bord_vol);
  set_borders_valid(conf);
  
  all_to_all_gathering_list_t gl;
  for(int ibord=0;ibord<bord_vol;ibord++)
    for(int mu=0;mu<4;mu++)
      {
	int idest=mu+4*(ibord+loc_vol);
	int iglb=glblx_of_bordlx[ibord];
	coords c;
	glb_coord_of_glblx(c,iglb);
	
	int isrc,irank;
	get_loclx_and_rank_of_coord(&isrc,&irank,c);
	isrc=isrc*4+mu;
	
	gl[std::make_pair(rank,isrc)]=idest;
      }
  all_to_all_comm_t comm(gl);
  comm.communicate(conf,conf,sizeof(su3));
  global_plaquette_and_rectangles_lx_conf(paths,conf);
  
  //polyakov loop
  complex pol;
  average_polyakov_loop_lx_conf(pol,conf,0);
  
  master_fprintf(file,"%d\t%016.16lg\t%016.16lg\t%+016.16lg\t%+016.16lg\n",itraj,paths[0],paths[1],pol[0],pol[1]);
  
  if(rank==0) fclose(file);
}

//store conf when appropriate
void store_conf_if_necessary()
{
  if(store_conf_each!=0 && itraj%store_conf_each==0)
    {
      char path[1024];
      sprintf(path,"%s.%05d",store_conf_path,itraj);
      write_conf();
    }
}

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  ///////////////////////////////////////
  
  //evolve for the required number of traj
  prod_ntraj=0;
  master_printf("\n");
  do
    {
      // 1) produce new conf
      if(max_ntraj!=0)
	{
	  generate_new_conf();
	  prod_ntraj++;
	  itraj++;
	}
      
      // 2) measure
      measure_gauge_obs();
      
      // 3) increment id and write conf
      if(store_running_temp_conf) write_conf();
      
      // 4) if conf is multiple of store_conf_each copy it
      store_conf_if_necessary();
      
      // 5) spacing between output
      master_printf("\n");
      
      //surely now we have created conf
      conf_created=0;
    }
  while(prod_ntraj<max_ntraj && !file_exists("stop") && !file_exists("restart"));
  
  /////////////////////////////////////// timings /////////////////////////////////
  
  close_simulation();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
