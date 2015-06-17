#include "nissa.hpp"

using namespace nissa;

tm_quark_info le;
double lep_energy,neu_energy;

//compute phase exponent for space part: vec{p}*\vec{x}
double get_space_arg(int ivol,momentum_t bc)
{
  double arg=0;
  for(int mu=1;mu<4;mu++)
    {
      double step=bc[mu]*M_PI/glb_size[mu];
      arg+=step*glb_coord_of_loclx[ivol][mu];
    }
  return arg;
}

//compute the phase for lepton on its sink
void get_lepton_sink_phase_factor(complex out,int ivol)
{
  //compute space and time factor
  double arg=get_space_arg(ivol,le.bc);
  int t=glb_coord_of_loclx[ivol][0];
  double ext=exp(t*lep_energy);
  
  //compute full exponential (notice the factor -1)
  out[RE]=cos(-arg)*ext;
  out[IM]=sin(-arg)*ext;
}

//set everything to a phase factor
void set_to_lepton_sink_phase_factor(spinspin *prop,int twall)
{
  GET_THREAD_ID();
  
  vector_reset(prop);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    if(twall==-1||glb_coord_of_loclx[ivol][0]==twall)
      //if(glb_coord_of_loclx[ivol][0]==0)
      {
	complex ph;
	get_lepton_sink_phase_factor(ph,ivol);
	spinspin_put_to_diag(prop[ivol],ph);
      }
  set_borders_invalid(prop);
}

//compute the phase for antineutrino - the orientation is that of the muon (as above)
void get_antineutrino_source_phase_factor(complex out,int ivol,momentum_t bc)
{
  //compute space and time factor
  double arg=get_space_arg(ivol,bc);
  int t=glb_coord_of_loclx[ivol][0];
  if(t>=glb_size[0]/2) t=glb_size[0]-t;
  double ext=exp(t*neu_energy);
  
  //compute full exponential (notice the factor +1)
  out[RE]=cos(+arg)*ext;
  out[IM]=sin(+arg)*ext;
}

THREADABLE_FUNCTION_0ARG(compute_lepton_free_loop)
{
  GET_THREAD_ID();
  
  FILE *fout=open_file("corr_l_free","w");
  
  spinspin *prop=nissa_malloc("prop",loc_vol+bord_vol,spinspin);
  complex *corr=nissa_malloc("corr",glb_size[0],complex);
  
  //put it to a phase
  int twall=glb_size[0]/2;
  set_to_lepton_sink_phase_factor(prop,twall);
  
  //multiply with the prop
  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le);
  
  double summu=lep_energy+neu_energy;
  double A02=2*sqr(lep_energy)*(1-sqr(lep_energy/summu));
  master_printf("%+016.016lg\n",A02);
  
  //get the projectors
  spinspin promu[2],pronu[2];
  twisted_on_shell_operator_of_imom(promu[0],le,0,false,-1);
  twisted_on_shell_operator_of_imom(promu[1],le,0,false,+1);
  naive_massless_on_shell_operator_of_imom(pronu[0],le.bc,0,-1);
  naive_massless_on_shell_operator_of_imom(pronu[1],le.bc,0,+1);
  
  //compute the right part of the leptonic loop: G0 G^dag
  const int nhadrolept_proj=2,hadrolept_projs[nhadrolept_proj]={9,4};
  dirac_matr hadrolept_proj_gamma[nhadrolept_proj];
  for(int ig_proj=0;ig_proj<nhadrolept_proj;ig_proj++)
    {
      int ig=hadrolept_projs[ig_proj];
      dirac_matr temp_gamma;
      dirac_herm(&temp_gamma,base_gamma+ig);
      dirac_prod(hadrolept_proj_gamma+ig_proj,base_gamma+map_mu[0],&temp_gamma);
    }
  
  const int nweak_ins=2;
  int list_weak_insl[nweak_ins]={4,9};
  for(int ins=0;ins<nweak_ins;ins++)
    {
      //define a local storage
      spinspin l_loc_corr[loc_size[0]];
      for(int i=0;i<loc_size[0];i++) spinspin_put_to_zero(l_loc_corr[i]);
      
      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	{
	  int t=loc_coord_of_loclx[ivol][0];
	  
	  //multiply lepton side on the right (source) side
	  spinspin l;
	  unsafe_spinspin_prod_dirac(l,prop[ivol],base_gamma+list_weak_insl[ins]);
	  
	  //add the neutrino phase
	  complex ph;
	  get_antineutrino_source_phase_factor(ph,ivol,le.bc);
	  spinspin_summ_the_complex_prod(l_loc_corr[t],l,ph);
	}
      glb_threads_reduce_double_vect((double*)l_loc_corr,loc_size[0]*sizeof(spinspin)/sizeof(double));
      
      //save projection on LO
      for(int ig_proj=0;ig_proj<nhadrolept_proj;ig_proj++)
	{
	  vector_reset(corr);
	  NISSA_PARALLEL_LOOP(loc_t,0,loc_size[0])
	    {
	      int glb_t=loc_t+glb_coord_of_loclx[0][0];
	      int ilnp=(glb_t>=glb_size[0]/2); //select the lepton/neutrino projector
	      
	      spinspin td;
	      unsafe_spinspin_prod_spinspin(td,l_loc_corr[loc_t],pronu[ilnp]);
	      //unsafe_dirac_prod_spinspin(td,base_gamma+list_weak_insl[ins],pronu[ilnp]);
	      spinspin dtd;
	      unsafe_spinspin_prod_spinspin(dtd,promu[ilnp],td);
	      trace_spinspin_with_dirac(corr[glb_t],dtd,hadrolept_proj_gamma+ig_proj);
	    }
	  THREAD_BARRIER();
	  
	  if(IS_MASTER_THREAD)
	    {
	      glb_nodes_reduce_complex_vect(corr,glb_size[0]);
	      master_fprintf(fout," # ins=%s, ig_proj=%s\n\n",gtag[list_weak_insl[ins]],gtag[hadrolept_projs[ig_proj]]);
	      for(int t=0;t<glb_size[0];t++)
		{
		  int mt=(t<glb_size[0]/2?t:glb_size[0]-t);
		  double n=exp(-mt*(lep_energy+neu_energy))/glb_vol*glb_size[0];
		  //n=1;
		  master_fprintf(fout,"%+016.016lg %+016.016lg\n",corr[t][0]*n/A02,corr[t][1]*n/A02);
		}
	      master_fprintf(fout,"\n");
	    }
	}
    }
  
  nissa_free(prop);
  nissa_free(corr);
  close_file(fout);
}
THREADABLE_FUNCTION_END

void in_main(int narg,char **arg)
{
  if(narg<2) crash("use %s nx",arg[0]);
  int X=atoi(arg[1]);
  int T=X,L=X;
  init_grid(T,L);
  
  //prepare the quark info
  le.r=0;
  le.mass=0.1;
  le.kappa=0.125;
  le.bc[0]=1;
  for(int mu=1;mu<4;mu++) le.bc[mu]=0.3;
  //compute energies
  lep_energy=tm_quark_energy(le,0);
  neu_energy=naive_massless_quark_energy(le.bc,0);
  master_printf("mcrit: %lg, Emu: %+016.016lg, Enu: %+016.016lg\n",m0_of_kappa(le.kappa),lep_energy,neu_energy);
    
  compute_lepton_free_loop();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
