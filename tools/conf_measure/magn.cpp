#include "nissa.hpp"

using namespace nissa;

//fixed to Z magnetization
const int mu=1,nu=2;

int L,T;
int nconfs,nflavs,nintervals,nhits;
double b_min,b_max;

char out_path[1024];
quark_content_t *quark_content;
stout_pars_t stout_pars;
quad_su3 *base_conf_eo[2],*stout_conf_eo[2];
color *chi[2],*rnd[2],*guess;
coords *phases;
complex *point_magn;
quad_u1 *u1b[2];

FILE *outfile;

//compute the magnetization
THREADABLE_FUNCTION_4ARG(compute_mag, quark_content_t*,qc, int,iconf, int,ihit, int,iflav)
{
  //reset initial guess
  vector_reset(guess);
  
  for(int ip=0;ip<=nintervals;ip++)
    {
      //define the value of b
      double b=b_min+(b_max-b_min)*ip/nintervals;
      master_printf("  B=%d/%d (%lg)\n",ip,nintervals,b);
      
      //initialize background field to id, then add magnetic field
      init_backfield_to_id(u1b);
      add_em_field_to_backfield(u1b,qc,b,mu,nu);
      
      //set the configuration background
      add_backfield_to_conf(stout_conf_eo,u1b);
      
      //invert and copy the guess for next iteration
      inv_stD_cg(chi,guess,stout_conf_eo,qc->mass,1000000,1e-12,rnd);
      vector_copy(guess,chi[EVN]);
      
      //compute the magnetization
      complex magn;
      magnetization(&magn,stout_conf_eo,qc,rnd,chi,point_magn,phases,mu,nu);
      
      //remove the background field
      rem_backfield_from_conf(stout_conf_eo,u1b);
      
      master_fprintf(outfile,"%d %d %d %lg %16.16lg\n",iconf,ihit,iflav,b,magn[0]);
    }
}
THREADABLE_FUNCTION_END

//analyze a single conf
void analyze_conf(int iconf)
{
  //read the configuration
  char conf_path[1024];
  read_str(conf_path,1024);
  read_ildg_gauge_conf_and_split_into_eo_parts(base_conf_eo,conf_path);
  
  //smear the configuration
  stout_smear(stout_conf_eo,base_conf_eo,&stout_pars);
  
  //loop over the number of hits
  for(int ihit=0;ihit<nhits;ihit++)
    {
      master_printf("Hit %d/%d\n",ihit+1,nhits);
      
      //generate the source
      generate_fully_undiluted_eo_source(rnd,RND_GAUSS,-1);
      
      //loop over all the flavors
      for(int iflav=0;iflav<nflavs;iflav++)
	{
	  master_printf(" Flav %d/%d\n",iflav+1,nflavs);
	  compute_mag(quark_content+iflav,iconf,ihit,iflav);
	}
    }
}

//initialize
void init(const char *input_path)
{
  GET_THREAD_ID();
  
  //open input file
  open_input(input_path);
  
  //init the grid 
  read_str_int("L",&L);
  read_str_int("T",&T);

  //read seed
  int seed;
  read_str_int("Seed",&seed);
  
  //read flavor parameters
  read_str_int("NDiffFlavs",&nflavs);
  quark_content=nissa_malloc("quark_content",nflavs,quark_content_t);
  for(int iflav=0;iflav<nflavs;iflav++) read_quark_content(quark_content[iflav]);

  read_stout_pars(stout_pars);  //stout parameters
  read_str_double("BMin",&b_min); //extrema1
  read_str_double("BMax",&b_max); //extrema2
  read_str_int("NIntervals",&nintervals); //number of points
  read_str_str("OutPath",out_path,1024); //output path
  read_str_int("NHits",&nhits);   //read the number of hits
  read_str_int("NConfs",&nconfs); //read the number of configurations
  
  //start grid and loc rnd gen
  init_grid(T,L);  
  start_loc_rnd_gen(seed);
  
  //allocate
  point_magn=nissa_malloc("app",loc_vol,complex);
  phases=nissa_malloc("arg",loc_vol+bord_vol,coords);
  guess=nissa_malloc("guess",loc_volh+bord_volh,color);
  for(int eo=0;eo<2;eo++)
    {
      u1b[eo]=nissa_malloc("u1b",loc_vol+bord_vol,quad_u1);  
      chi[eo]=nissa_malloc("chi",loc_volh+bord_volh,color);
      rnd[eo]=nissa_malloc("rnd",loc_volh+bord_volh,color);
      base_conf_eo[eo]=nissa_malloc("base_conf_eo",loc_volh+bord_volh+edge_volh,quad_su3);
      stout_conf_eo[eo]=nissa_malloc("stout_conf_eo",loc_volh+bord_volh+edge_volh,quad_su3);
    }
  
  //open output
  outfile=open_file(out_path,"w");

  //store phases
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol+bord_vol)
    get_args_of_one_over_L2_quantization(phases[ivol],ivol,mu,nu);    
}

//close
void stop()
{
  close_input();
  nissa_free(quark_content);
  nissa_free(point_magn);
  nissa_free(phases);
  nissa_free(guess);
  for(int eo=0;eo<2;eo++)
    {
      nissa_free(u1b[eo]);
      nissa_free(chi[eo]);
      nissa_free(rnd[eo]);
      nissa_free(base_conf_eo[eo]);
      nissa_free(stout_conf_eo[eo]);
    }
}

void in_main(int narg,char **arg)
{
  //initialize
  if(narg<2) crash("use: %s input",arg[0]);
  init(arg[1]);

  //analyze the configurations one by one
  for(int iconf=0;iconf<nconfs;iconf++) analyze_conf(iconf);
  
  //close everything
  stop();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
