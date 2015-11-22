#include "nissa.hpp"

#define PROP_TYPE colorspinspin

#include "conf.hpp"
#include "pars.hpp"
#include "prop.hpp"

using namespace nissa;

spincolor *eta,*phi,*temp;
spin1field **bubble,*temp_bubble;
spin1field **tot_bubble;
std::vector<double> tadpole_val;

int nB;
int iB(int imass,int r)
{return r+2*imass;}

//initialize the simulation
void init_simulation(const char *path)
{
  open_input(path);
  
  read_input_preamble();
  read_photon_pars();
  read_seed_start_random();
  read_noise_type();
  read_free_theory_flag();
  read_random_gauge_transform();
  read_loc_pion_curr();
  read_nsources();
  read_ngauge_conf();
  
  //allocate
  nB=iB(nqmass-1,nr-1)+1;
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  eta=nissa_malloc("eta",loc_vol+bord_vol,spincolor);
  phi=nissa_malloc("phi",loc_vol+bord_vol,spincolor);
  temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  temp_bubble=nissa_malloc("temp_bubble",loc_vol+bord_vol,spin1field);
  bubble=nissa_malloc("bubble*",nB,spin1field*);
  tot_bubble=nissa_malloc("tot_bubble*",nB,spin1field*);
  for(int iB=0;iB<nB;iB++)
    {
      bubble[iB]=nissa_malloc("bubble",loc_vol+bord_vol,spin1field);
      tot_bubble[iB]=nissa_malloc("tot_bubble",loc_vol+bord_vol,spin1field);
    }
}

//init a new conf
void start_new_conf()
{
  setup_conf(conf,old_theta,put_theta,conf_path,rnd_gauge_transform,free_theory);
}

//what to do to skip a configuration
void skip_conf()
{
}

//close deallocating everything
void close()
{
  //free
  nissa_free(eta);
  nissa_free(phi);
  nissa_free(temp);
  nissa_free(conf);
  for(int iB=0;iB<nB;iB++)
    {
      nissa_free(tot_bubble[iB]);
      nissa_free(bubble[iB]);
    }
  nissa_free(temp_bubble);
  nissa_free(tot_bubble);
  nissa_free(bubble);
  
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g, of which:\n",tot_prog_time);
}

//tadpole
THREADABLE_FUNCTION_3ARG(compute_tadpole, spincolor*,phi, spincolor*,eta, int,r)
{
  GET_THREAD_ID();
  
  if(!pure_wilson) insert_tm_tadpole(temp,conf,phi,r,tadpole,-1);
  else             insert_wilson_tadpole(temp,conf,phi,tadpole,-1);
  
  //trace
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<NCOL;ic++)
	safe_complex_conj1_prod(temp[ivol][id][ic],eta[ivol][id][ic],temp[ivol][id][ic]);
  THREAD_BARRIER();
  
  //reduce and store
  complex tadpole_res;
  complex_vector_glb_collapse(tadpole_res,(complex*)temp,sizeof(spincolor)/sizeof(complex)*loc_vol);
  tadpole_val.push_back(tadpole_res[RE]);
  tadpole_val.push_back(tadpole_res[IM]);
  
  master_printf("%+16.16lg %+16.16lg\n",tadpole_res[0],tadpole_res[1]);
}
THREADABLE_FUNCTION_END

//matrix element of the vector current as a function of the photon position
THREADABLE_FUNCTION_5ARG(vector_matrix_element, spin1field*,out, spincolor*,sink, quad_su3*,conf, int,r, spincolor*,source)
{
  GET_THREAD_ID();
  
  //set parameters
  dirac_matr GAMMA;
  if(!pure_wilson) dirac_prod_idouble(&GAMMA,base_gamma+5,-tau3[r]);
  else             GAMMA=base_gamma[0];
  
  //reset the output and communicate borders
  vector_reset(out);
  communicate_lx_spincolor_borders(sink);
  communicate_lx_spincolor_borders(source);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<NDIM;mu++)
      {
	//find neighbors
	int ifw=loclx_neighup[ivol][mu];
	
	//transport down and up
	spincolor fw,bw;
	unsafe_su3_prod_spincolor(fw,conf[ivol][mu],source[ifw]);
	unsafe_su3_dag_prod_spincolor(bw,conf[ivol][mu],source[ivol]);
	
	//put GAMMA-gmu on the forward piece
	spincolor GA_fw,gmu_fw,GA_M_gmu_fw;
	unsafe_dirac_prod_spincolor(GA_fw,&GAMMA,fw);
	unsafe_dirac_prod_spincolor(gmu_fw,base_gamma+map_mu[mu],fw);
	spincolor_subt(GA_M_gmu_fw,GA_fw,gmu_fw);
	//put GAMMA+gmu on the backward piece
	spincolor GA_bw,gmu_bw,GA_P_gmu_bw;
	unsafe_dirac_prod_spincolor(GA_bw,&GAMMA,bw);
	unsafe_dirac_prod_spincolor(gmu_bw,base_gamma+map_mu[mu],bw);
	spincolor_summ(GA_P_gmu_bw,GA_bw,gmu_bw);
	
	//take the scalar product
	complex c_fw,c_bw;
	spincolor_scalar_prod(c_fw,sink[ivol],GA_M_gmu_fw);
	spincolor_scalar_prod(c_bw,sink[ifw],GA_P_gmu_bw);
	
	//summ the result
	complex fact_fw={0,-0.5},fact_bw={0,+0.5};
	complex_summ_the_prod(out[ivol][mu],c_fw,fact_fw);
	complex_summ_the_prod(out[ivol][mu],c_bw,fact_bw);
      }
  
  set_borders_invalid(out);
}
THREADABLE_FUNCTION_END

//compute everything with a single mass
void single_mass(int imass,int r)
{
  master_printf(" imass %d/%d, r %d/%d\n",imass+1,nqmass,r,nr);
  
  //solve for phi for each quark
  get_qprop(phi,eta,imass,r);
  
  //compute the tadpole
  compute_tadpole(phi,eta,r);
  
  //compute the single bubble
  vector_matrix_element(bubble[iB(imass,r)],eta,conf,r,phi);
  
  //debug
  complex tot={0,0};
  for(int mu=0;mu<4;mu++)
    for(int ivol=0;ivol<loc_vol;ivol++)
      complex_summassign(tot,bubble[iB(imass,r)][ivol][mu]);
  master_printf("AAAAA %lg %lg\n",tot[RE],tot[IM]);
  
  insert_tm_external_source(temp,conf,NULL,phi,r,-1);
  complex tot2={0,0};
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<NCOL;ic++)
	complex_summ_the_conj1_prod(tot2,eta[ivol][id][ic],temp[ivol][id][ic]);
  master_printf("BBBB %lg %lg\n",tot2[RE],tot2[IM]);
}

//compute everything with a single source
void single_source(int isource)
{
  master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
  
  //generate the source
  generate_undiluted_source(eta,rnd_type_map[noise_type],-1);
  
  for(int imass=0;imass<nqmass;imass++)
    for(int r=0;r<nr;r++)
      single_mass(imass,r);
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
  while(iconf<ngauge_conf && enough_time && !file_exists("stop") && read_conf_parameters(iconf,skip_conf))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      //loop over sources
      for(int isource=0;isource<nsources;isource++)
	single_source(isource);
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
