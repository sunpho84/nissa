#include "nissa.hpp"

struct source_t
{
  char name[1024];
  int type;
  char get_source_from_propgroup_of_name[1024];
  double get_source_from_propgroup_theta;
  double get_source_from_propgroup_mass;
  int select_timeslice;
  int gamma_to_be_inserted;
  int noise_type;
  int save_source;
  colorspinspin *data;
};

struct propgroup_t
{
  char name[1024];
  int nmass;
  double *mass;
  int ntheta;
  double *theta;
  int both_r;
  int load_props;
  int save_props;
  colorspinspin **data;
};

int wall_time,seed,ape_niter;
double kappa,ape_alpha;

int nsource;
source_t *sources;

//read a string and return its position inside a list, crash if not found
int read_str_enum(int n,char **tags,int stl)
{
  char read[stl];
  read_str(read,stl);
  
  int found=-1;
  for(int i=0;i<n;i++)
    if(strcasecmp(read,tags[i])==0)
      found=i;
  
  if(found==-1)
    {
      master_fprintf(stderr,"unkwnown tag %s, expected:\n",read);
      for(int i=0;i<n;i++)
	master_fprintf(stderr," %s\n",tags[i]);
      CRASH("check input");
    }
  
  return found;
}

//initialize the program
void init_semileptonic(int narg,char **arg)
{
  init_nissa();
  
  if(narg<2) CRASH("ue: %s input",arg[0]);
  
  open_input(arg[1]);

  // 1) Read information about the lattice
  
  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Smearing parameters of the conf
  read_str_double("ApeAlpha",&ape_alpha);
  read_str_int("ApeNiter",&ape_niter);
  
  // 2) Information about the sources to be generated
  
  //Number of different sources to be generated
  read_str_int("NSource",&nsource);
  sources=nissa_malloc("Sources_info",nsource,source_t);
  for(int isource=0;isource<nsource;isource++)
    {
      //Source name
      read_str_str("SourceName",sources[isource].name,1024);
      //Type: generate(0), load(1) or sequential(2)
      char type_tag[3][20]={"generate","load","sequential"};
      int type=read_str_enum(3,(char**)type_tag,20);
      //Noise type
      read_str_int("NoiseType",&sources[isource].noise_type);
      //Save the source?
      read_str_int("SaveSource",&sources[isource].save_source);
    }
  
  // 3) Information about propagators to be generated
  
  //Number of different propagators group
  
  
  // 3) Read information about the run
  
  //Wall time
  read_str_int("WallTime",&wall_time);
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  //Kappa
  read_str_double("Kappa",&kappa);
  
}

//close the program
void close_semileptonic()
{
  close_nissa();
}

int main(int narg,char **arg)
{
  //init the program
  init_semileptonic(narg,arg);
  
  //close the program
  close_semileptonic();
  
  return 0;
}
