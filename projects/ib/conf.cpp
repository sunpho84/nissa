#include <nissa.hpp>

namespace nissa
{
  int nanalyzed_conf=0;
  double tot_prog_time=0,wall_time;
  
  //init the MPI grid
  void read_init_grid()
  {
    int L,T;
    read_str_int("L",&L);
    read_str_int("T",&T);
    
    init_grid(T,L);
  }
  
  //read the conf and setup it
  void setup_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta,const char *conf_path,int rnd_gauge_transform,int free_theory)
  {
    //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
    if(!free_theory)
      {
	read_ildg_gauge_conf(conf,conf_path);
	master_printf("plaq: %+016.016g\n",global_plaquette_lx_conf(conf));
      }
    else generate_cold_lx_conf(conf);
    
    //if asked, randomly transform the configurations
    if(rnd_gauge_transform) perform_random_gauge_transform(conf,conf);
    
    //put anti-periodic boundary condition for the fermionic propagator
    old_theta[0]=old_theta[1]=old_theta[2]=old_theta[3]=0;
    put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
    adapt_theta(conf,old_theta,put_theta,0,0);
  }
  
  //used to shift the configuration
  void index_shift(int &irank_out,int &ivol_out,int ivol_in,void *pars)
  {
    int *source_coord=(int*)pars;
    coords co;
    for(int nu=0;nu<NDIM;nu++) co[nu]=(glb_coord_of_loclx[ivol_in][nu]+source_coord[nu])%glb_size[nu];
    get_loclx_and_rank_of_coord(&ivol_out,&irank_out,co);
  }
  
  //perform a random shift
  void random_shift_gauge_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta)
  {
    //remove phase
    put_theta[0]=0;put_theta[1]=put_theta[2]=put_theta[3]=0;
    adapt_theta(conf,old_theta,put_theta,0,0);
    
    //source coord
    coords shift_coord;
    generate_random_coord(shift_coord);
    
    //shift the configuration
    double shift_time=-take_time();
    vector_remap_t shifter(loc_vol,index_shift,(void*)shift_coord);
    shifter.remap(conf,conf,sizeof(quad_su3));
    shift_time+=take_time();
    master_printf("Shifted of %d %d %d %d in %lg sec, plaquette after shift: %+016.016lg\n",shift_coord[0],shift_coord[1],shift_coord[2],shift_coord[3],shift_time,global_plaquette_lx_conf(conf));
    
    //put back the phase
    put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
    adapt_theta(conf,old_theta,put_theta,0,0);
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
}
