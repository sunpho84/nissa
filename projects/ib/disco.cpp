#include "nissa.hpp"

#include "conf.hpp"
#include "pars.hpp"
#include "prop.hpp"

using namespace nissa;

qprop_t q;

//compute the matrix element of the conserved current between two propagators
THREADABLE_FUNCTION_4ARG(calc_cur, spin1field*,cur, spincolor*,source, quad_su3*,conf, spincolor*,prop)
{
  GET_THREAD_ID();
  
  vector_reset(cur);
  
  communicate_lx_spincolor_borders(source);
  communicate_lx_spincolor_borders(prop);
  communicate_lx_quad_su3_borders(conf);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<NDIM;mu++)
      {
	// int ifw=loclx_neighup[ivol][mu];
	// spincolor f;
	
	// //piece psi_ivol U_ivol psi_fw
	// unsafe_su3_prod_spincolor(f,conf[ivol][mu],prop[ifw]);
	// spincolor_scalar_prod(cur[ivol][mu],source[ivol],f);
	
	int ifw=loclx_neighup[ivol][mu];
	spincolor f,Gf;
	complex c;
	
	//piece psi_ivol U_ivol psi_fw
	unsafe_su3_prod_spincolor(f,conf[ivol][mu],prop[ifw]);
	spincolor_copy(Gf,f);
	dirac_subt_the_prod_spincolor(Gf,base_gamma+igamma_of_mu[mu],f);
	spincolor_scalar_prod(c,source[ivol],Gf);
	complex_summ_the_prod_idouble(cur[ivol][mu],c,-0.5);
	
	//piece psi_fw U_ivol^dag psi_ivol
	unsafe_su3_dag_prod_spincolor(f,conf[ivol][mu],prop[ivol]);
	spincolor_copy(Gf,f);
	dirac_summ_the_prod_spincolor(Gf,base_gamma+igamma_of_mu[mu],f);
	spincolor_scalar_prod(c,source[ifw],Gf);
	complex_summ_the_prod_idouble(cur[ivol][mu],c,+0.5);
      }
  set_borders_invalid(cur);
}
THREADABLE_FUNCTION_END

void in_main_Marcus(int narg,char **arg)
{
  init_grid(16,8);
  
  spincolor *source=nissa_malloc("source",loc_vol,spincolor);
  spincolor *prop=nissa_malloc("prop",loc_vol,spincolor);
  quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);
  
  read_ildg_gauge_conf(conf,"conf");
  momentum_t th={1,0,0,0};
  put_boundaries_conditions(conf,th,0,0);
  
  read_real_vector(source,"source.0000.00000","scidac-binary-data");
  read_real_vector(prop,"source.0000.00000.inverted","scidac-binary-data");
  
  spin1field *cur=nissa_malloc("cur",loc_vol,spin1field);
  calc_cur(cur,source,conf,prop);
  
  FILE *fout=open_file("cur.txt","w");
  master_fprintf(fout,"# [loops_em]  norm square of source     = %.16lg\n",double_vector_glb_norm2(source,loc_vol));
  master_fprintf(fout,"# [loops_em]  norm square of propagator = %.16lg\n",double_vector_glb_norm2(prop,loc_vol));
  NISSA_LOC_VOL_LOOP(ivol)
  {
    master_fprintf(fout,"# [loops_em] x\t");
    for(int mu=0;mu<NDIM;mu++) master_fprintf(fout,"%d ",loc_coord_of_loclx[ivol][mu]);
    master_fprintf(fout,"\n");
    
    for(int mu=0;mu<NDIM;mu++) master_fprintf(fout,"%d" "\t" "%.16lg" "\t" "%.16lg\n",mu,cur[ivol][mu][RE],cur[ivol][mu][IM]);
  }
  close_file(fout);
  
  nissa_free(cur);
  nissa_free(conf);
  nissa_free(prop);
  nissa_free(source);
}

void init_simulation(int narg,char **arg)
{
  std::string input_path="";
  
  //parse opts
  int c;
  while((c=getopt(narg,arg,"i:"))!= -1)
    switch (c)
      {
      case 'i': input_path=optarg; break;
      default: crash("Unknown option -%c",optopt);
      }
  
  if(input_path=="") crash("Please specify -i");
  open_input(input_path);
  
    //init the grid
  read_init_grid();
  
  //Wall time
  read_str_double("WallTime",&wall_time);
  
  //Sources
  read_seed_start_random();
  read_stoch_source();
  set_diluted_spin(1);
  set_diluted_color(1);
  set_diluted_space(1);
  read_nhits();

  read_twisted_run();
  read_clover_run();
  
  /////////////////////////////////////////////////////////////////

  int tins=-1;
  double kappa=0.125,mass=0.0,charge=0,theta[NDIM],residue=1e-16;
  theta[0]=temporal_bc;
  for(int mu=1;mu<NDIM;mu++) theta[mu]=0;
  int r=0,store_prop=0;
      
  read_str_double("Kappa",&kappa);
  if(twisted_run)
    {
      read_str_double("Mass",&mass);
      master_printf("Read variable 'Mass' with value: %lg\n",mass);
      read_int(&r);
      master_printf("Read variable 'R' with value: %d\n",r);
      
      //include tau in the mass
      mass*=tau3[r];
    }
  
  read_str_int("StoreProp",&store_prop);
  master_printf("Read variable 'Store' with value: %d\n",store_prop);
  q.init_as_propagator(ins_from_tag("-"),{{"source",{1.0,0.0}}},tins,residue,kappa,mass,r,charge,theta,store_prop);
  read_double(&residue);
  master_printf("Read variable 'Residue' with value: %lg\n",residue);
  
  read_photon_pars();
  
  read_free_theory_flag();
  read_random_gauge_transform();
  read_Landau_gauge_fix();
  read_store_conf();
  
  read_loc_hadr_curr();
  
  read_ngauge_conf();
}

void close()
{
  close_input();
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //init simulation according to input file
  init_simulation(narg,arg);
  
  //loop over the configs
  int iconf=0;
  while(read_conf_parameters(iconf,finish_file_present))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      for(int ihit=0;ihit<nhits;ihit++)
	{
	  start_hit(ihit);
	  generate_propagators(ihit);
	  //compute_contractions();
	  propagators_fft(ihit);
	}
      //print_contractions();
      
      mark_finished();
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
