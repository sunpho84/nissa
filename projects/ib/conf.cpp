#include <nissa.hpp>

namespace nissa
{
  
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
  
  //generate a random postion
  void generate_random_coord(coords c)
  {for(int mu=0;mu<NDIM;mu++) c[mu]=(int)(rnd_get_unif(&glb_rnd_gen,0,glb_size[mu]));}
  
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
}
