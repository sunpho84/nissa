#include <nissa.h>

/*
int T=48,L=24;
int tsource=0;
double kappa=0.160856,mass=0.0100,residue=1e-14;
char conf_path[]="test/data/L4T8conf";
*/

typedef double compo[16];

//pars
int T,L,seed,tsource,npassed_conf,ntot_conf;
int rprop=1;
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
  read_str_double("Residue",&residue);
  read_str_str("OutPath",out_path,1024);
  read_str_double("Walltime",&wall_time);
  read_str_int("NTotConf",&ntot_conf);
}

//compute the number of compute configurations
int count_npassed_conf()
{
  //if out file does not exist return 0
  if(!file_exists(out_path)) return 0;
  
  //check that number of lines is a multipe of T
  int n=count_file_lines(out_path);
  if(n%T) crash("outfile %s contains %d lines which is not multiple of T=%d",out_path,n,T);
  
  return n/T;
}

//initalize the program
void init(int narg,char **arg)
{
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
  prop=solution[rprop][0];
}

//read a conf
void prepare_conf()
{
  //load conf, compute plaquette and put anti-periodic boundary conditions in time
  read_ildg_gauge_conf(conf,conf_path);
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  double put_theta[4]={1,0,0,0},old_theta[4]={0,0,0,0};
  adapt_theta(conf,old_theta,put_theta,1,1);
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

//invert the prop
void get_prop()
{
  compute_su3spinspin_tm_propagators_multi_mass((su3spinspin***)solution,conf,kappa,&mass,1,10000,&residue,source);
  rotate_vol_su3spinspin_to_physical_basis(prop,!rprop,!rprop);
}

//project a single su3spinspin over all gamma matrices
void point_gamma_proj_summ(compo out,su3spinspin in)
{
  for(int igamma=0;igamma<16;igamma++)
    for(int ic1=0;ic1<3;ic1++)
      for(int ic2=0;ic2<3;ic2++)
	  {
	    complex t={0,0};
	    for(int id1=0;id1<4;id1++)
	      {
		int id2=base_gamma[igamma].pos[id1];
		complex_summ_the_prod(t,base_gamma[igamma].entr[id1],in[ic1][ic2][id2][id1]);
	      }
	    out[igamma]+=squared_complex_norm(t);
	  }
}

//projects over all gamma matrices
void gamma_proj()
{
  FILE *fout=open_file(out_path,"aw");
  
  //buffers
  compo loc_out[T],glb_out[T];
  memset(loc_out,0,sizeof(compo)*T);
  
  //summ all lattice
  nissa_loc_vol_loop(ivol)
    point_gamma_proj_summ(loc_out[(loc_coord_of_loclx[ivol][0]-tsource+T)%T],prop[ivol]);  
  MPI_Reduce(loc_out,glb_out,T*16,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  //print
  for(int t=0;t<T;t++)
    {
      for(int igamma=0;igamma<16;igamma++) master_fprintf(fout,"%16.16lg\t",glb_out[t][igamma]/4);
      master_fprintf(fout,"\n");
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
      read_str(conf_path,1024);
      generate_source();
      
      //works only on non passed conf
      if(go_on && iconf>=npassed_conf)
	{
	  prepare_conf();
	  
	  //get prop and trace over all gamma
	  get_prop();
	  gamma_proj();
	}
      
      iconf++;
    }
  
  finish();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  
  close_nissa();
  
  return 0;
}
