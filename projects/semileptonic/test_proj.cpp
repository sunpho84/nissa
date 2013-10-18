#include <nissa.hpp>

#define ASCII 1
#define NGAMMA 16

typedef double proj_prop_t[NGAMMA];

//pars
int T,L,seed,tsource,npassed_conf,ntot_conf,rprop,binary_ascii;
double kappa,mass,residue,init_time,wall_time;
char conf_path[1024],out_path[1024];

//fields
quad_su3 *conf;
su3spinspin *source,**solution[2],*prop;

//scan input file
void read_input_header(char *input_path)
{
  open_input(input_path);
  
  read_str_int("L",&L);
  read_str_int("T",&T);
  read_str_double("Kappa",&kappa);
  read_str_double("Mass",&mass);
  read_str_int("R",&rprop);
  read_str_double("Residue",&residue);
  read_str_str("OutPath",out_path,1024);
  read_str_double("Walltime",&wall_time);
  read_str_int("Seed",&seed);
  read_str_int("BinaryAscii",&binary_ascii);
  read_str_int("NTotConf",&ntot_conf);
}

//compute the number of compute configurations
int count_npassed_conf()
{
  //if out file does not exist return 0
  if(!file_exists(out_path)) return 0;
  
  //find conf number according to saving mode
  int nconf;
  if(binary_ascii==ASCII)
    {
      //nconf is given by the ratio between file length and nlines
      int nlines=count_file_lines(out_path);
      nconf=nlines/T;
      
      //check that number of lines is a multipe of T
      if(nconf*T!=nlines) crash("outfile %s contains %d lines which is not multiple of T=%d",out_path,nlines,T);
    }
  else
    {
      //nconf is given by the ratio between file_size and conf_size
      int conf_size=T*sizeof(proj_prop_t),file_size=get_file_size(out_path);
      nconf=file_size/conf_size;
      
      //check that file size is a multiple of conf_size
      if(nconf*conf_size!=file_size)
	crash("outfile %s has size %d which is not multiple of conf size %d",out_path,file_size,conf_size);
    }
  
  return nconf;
}

//initalize the program
void init(int narg,char **arg)
{
  //take initial time needed to check against wall_time
  init_time=take_time();
  
  //check arguments and read input header
  if(narg<2) crash("use %s input",arg[0]);
  read_input_header(arg[1]);
  
  //init cartesian grid
  init_grid(T,L);

  //count passed conf
  npassed_conf=count_npassed_conf();  
  master_printf("Already computed: %d conf\n",npassed_conf);
  
  //start local random generator
  start_loc_rnd_gen(seed);

  //allocate conf, source and solution
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  source=nissa_malloc("source",loc_vol+bord_vol,su3spinspin);
  for(int r=0;r<2;r++)
    {
      solution[r]=nissa_malloc("solution_r",1,su3spinspin*);
      solution[r][0]=nissa_malloc("solution",loc_vol+bord_vol,su3spinspin);  
    }
  
  //bind prop to one of the two r
  prop=solution[rprop][0];
}

//read conf path checking conf existence
void read_conf_path()
{
  read_str(conf_path,1024);
  if(!file_exists(conf_path)) crash("conf %s does not exist");
}

//generate a point source
void generate_source()
{
  //generate source position
  tsource=(int)rnd_get_unif(&glb_rnd_gen,0,T);
  coords O={tsource,0,0,0};
  for(int i=1;i<4;i++) O[i]=(int)rnd_get_unif(&glb_rnd_gen,0,L);
  
  //generate source
  generate_delta_source(source,O);
}

//load conf, compute plaquette and put anti-periodic boundary conditions in time
void prepare_conf()
{
  read_ildg_gauge_conf(conf,conf_path);
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  double put_theta[4]={1,0,0,0},old_theta[4]={0,0,0,0};
  adapt_theta(conf,old_theta,put_theta,1,1);
}

//invert the prop
void get_prop()
{
  compute_su3spinspin_tm_propagators_multi_mass((su3spinspin***)solution,conf,kappa,&mass,1,10000,&residue,source);
  rotate_vol_su3spinspin_to_physical_basis(prop,!rprop,!rprop);
}

//project a single su3spinspin over all gamma matrices
void point_gamma_proj_summ(proj_prop_t out,su3spinspin in)
{
  //loop over all the gamma and tracing over color the square of the trace over dirac
  for(int igamma=0;igamma<16;igamma++)
    for(int ic1=0;ic1<3;ic1++)
      for(int ic2=0;ic2<3;ic2++)
	  {
	    //dirac trace
	    complex t={0,0};
	    for(int id1=0;id1<4;id1++)
	      {
		int id2=base_gamma[igamma].pos[id1];
		complex_summ_the_prod(t,base_gamma[igamma].entr[id1],in[ic1][ic2][id2][id1]);
	      }
	    
	    //square
	    out[igamma]+=squared_complex_norm(t)/4;
	  }
}

//projects over all gamma matrices
void gamma_proj()
{
  //buffers
  proj_prop_t loc_out[T],glb_out[T];
  memset(loc_out,0,sizeof(proj_prop_t)*T);
  
  //summ all lattice
  NISSA_LOC_VOL_LOOP(ivol)
    point_gamma_proj_summ(loc_out[(glb_coord_of_loclx[ivol][0]-tsource+T)%T],prop[ivol]);  
  MPI_Reduce(loc_out,glb_out,T*NGAMMA,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  //print or write
  FILE *fout=open_file(out_path,"aw");
  if(binary_ascii==ASCII)  
    for(int t=0;t<T;t++)
      {
	for(int igamma=0;igamma<NGAMMA;igamma++) master_fprintf(fout,"%16.16lg\t",glb_out[t][igamma]);
	master_fprintf(fout,"\n");
      }
  else
    {
      if(little_endian) doubles_to_doubles_changing_endianess((double*)glb_out,(double*)glb_out,NGAMMA*T);
      if(fwrite(glb_out,sizeof(proj_prop_t),T,fout)!=(size_t)T) crash("while writing data");
    }  
  close_file(fout);
}

//close the program
void finish()
{
  close_input();

  //free
  nissa_free(source);
  for(int r=0;r<2;r++)
    {
      nissa_free(solution[r][0]);
      nissa_free(solution[r]);
    }
  nissa_free(conf);
}

void in_main(int narg,char **arg)
{
  init(narg,arg);
  
  //start from conf 0 for coherence
  int iconf=0,go_on=true;
  
  //work on all conf up to reaching ntot_conf or walltime
  while(go_on)
    {
      go_on=(take_time()-init_time<wall_time) && iconf<ntot_conf;
      
      //read the conf path and generate the source
      read_conf_path();
      generate_source();
      
      if(go_on)
	{
	  //works only on non passed conf
	  if(iconf>=npassed_conf)
	    {
	      prepare_conf();
	      
	      //get prop and trace over all gamma
	      get_prop();
	      gamma_proj();
	    }
	  else master_printf("Skipping conf %s because already passed\n",conf_path);
	}
      
      //increment conf id
      iconf++;
    }

  //check if we are 
  if(iconf==ntot_conf) master_printf("Finished all the %d confs\n",ntot_conf);
  else master_printf("Passed walltime\n");

  finish();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);  
  close_nissa();
  
  return 0;
}
