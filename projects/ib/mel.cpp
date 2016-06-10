#include <nissa.hpp>

#include "conf.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

/* the loop is normalised such that the physical rate at leading order
   is obtained multiplying the loop by Gf^2 fpi^2 * phi2 (space phase
   factor) which is (1-rl^2)/(16 pi mpi) where rl=ml/mpi, whereas the
   interference is obtained by the full mesolepton contrelation
   multiplied by 4*mpi*fpi*Gf^2*phi2 */

using namespace nissa;

/////////////////////////////////////// data //////////////////////////////

//list the 8 matrices to insert for the weak current
const int nweak_ins=17;
const int nweak_ind=9;
//const int nmeslep_proj=4,meslep_projs[nmeslep_proj]={9,4,5,0};
const int nmeslep_proj=1,meslep_projs[nmeslep_proj]={4};
const int list_weak_insq[nweak_ins]=     {1,2,3,4, 6,7,8,9,  1,2,3,4, 6,7,8,9, 5};
const int list_weak_insl[nweak_ins]=     {1,2,3,4, 6,7,8,9,  6,7,8,9, 1,2,3,4, 5};
const int list_weak_ind_contr[nweak_ins]={0,0,0,1, 2,2,2,3,  4,4,4,5, 6,6,6,7, 8};
const char list_weak_ind_nameq[nweak_ind][3]={"VK","V0","AK","A0","VK","V0","AK","A0","P5"};
const char list_weak_ind_namel[nweak_ind][3]={"VK","V0","AK","A0","AK","A0","VK","V0","V0"};
int nind;
spinspin *meslep_hadr_part;
complex *meslep_contr;

//parameters of the leptons
int *lep_contr_iq1;
int *lep_contr_iq2;

///////////////////////////////// initialise the library, read input file, allocate /////////////////////////////////////

//set the list of gamma
void set_mes_gamma_contr_list()
{
  for(int ig=0;ig<16;ig++) mes_gamma_list.push_back(idirac_pair_t(5,ig));    //P5GI
  mes_gamma_list.push_back(idirac_pair_t(4,4));                              //V0V0
  for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu,mu));   //VKVK
  for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu,mu+9)); //VKTK
  for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu+9,mu)); //TKVK
  for(int ig=10;ig<=15;ig++) mes_gamma_list.push_back(idirac_pair_t(ig,ig)); //TKTK,BKBK
}

void init_simulation(char *path)
{
  //open input file
  open_input(path);
  read_input_preamble();
  read_dilutions();
  read_stoch_source();
  
  //Leptons
  read_str_int("LeptonicContr",&nleptons);
  lep_contr_iq1=nissa_malloc("lep_contr_iq1",nleptons,int);
  lep_contr_iq2=nissa_malloc("lep_contr_iq2",nleptons,int);
  leps=nissa_malloc("leps",nleptons,tm_quark_info);
  lep_energy=nissa_malloc("lep_energy",nleptons,double);
  neu_energy=nissa_malloc("neu_energy",nleptons,double);
  if(!pure_wilson) expect_str("Q1Q2LepmassMesmass");
  else             expect_str("Q1Q2LepkappaMesmass");
  for(int il=0;il<nleptons;il++)
    {
      //read quarks identfiying the mesons
      read_int(lep_contr_iq1+il);
      read_int(lep_contr_iq2+il);
      
      //if not pure wilson read mass
      if(pure_wilson) leps[il].mass=0;
      else            read_double(&leps[il].mass);
      
      //antiperiodic or periodic
      leps[il].bc[0]=QUARK_BOUND_COND;
      
      //maximal twist (if tm), otherwise read kappa
      if(pure_wilson) read_double(&leps[il].kappa);
      else            leps[il].kappa=0.125;
      leps[il].r=0;
      
      //read the mass of the meson (that must have been determined outside)
      double mes_mass;
      read_double(&mes_mass);
      
      //set initial value of bc and check kinematic
      for(int i=1;i<NDIM;i++) leps[il].bc[i]=0;
      if(tm_quark_energy(leps[il],0)>=mes_mass) crash("initial state is lighter (%lg) than final state at rest (%lg)!",mes_mass,tm_quark_energy(leps[il],0));
      
      //compute meson momentum and bc
      double err;
      master_printf("Resolving kinematical condition for combination of quarks %d/%d\n",il+1,nleptons);
      do
      	{
      	  //compute the error
	  double lep_energy=tm_quark_energy(leps[il],0);
	  double neu_energy=naive_massless_quark_energy(leps[il].bc,0);
      	  err=lep_energy+neu_energy-mes_mass;
      	  //compute the derivative
      	  double eps=1e-8;
      	  for(int i=1;i<NDIM;i++) leps[il].bc[i]+=eps;
      	  double der=(tm_quark_energy(leps[il],0)+naive_massless_quark_energy(leps[il].bc,0)-mes_mass-err)/eps;
      	  for(int i=1;i<NDIM;i++) leps[il].bc[i]-=eps+err/der;
	  
      	  master_printf("  lep_e: %+010.10lg, neu_e: %+010.10lg, mes_mass: %lg, error: %lg, der: %lg\n",lep_energy,neu_energy,mes_mass,err,der);
      	}
      while(fabs(err)>1e-14);
      
      //write down energy
      lep_energy[il]=tm_quark_energy(leps[il],0);
      neu_energy[il]=naive_massless_quark_energy(leps[il].bc,0);
      master_printf(" ilepton %d, lepton energy: %lg, neutrino energy: %lg\n",il,lep_energy[il],neu_energy[il]);
      master_printf(" lep+neut energy: %lg\n",lep_energy[il]+neu_energy[il]);
      master_printf(" bc: %+016.016lg\n\n",leps[il].bc[1]);
    }
  
  read_photon_pars();
  read_use_photon_field();
  read_seed_start_random();
  read_noise_type();
  read_free_theory_flag();
  read_gospel_convention();
  read_random_gauge_transform();
  read_loc_hadr_curr();
  read_loc_muon_curr();
  read_nsources();
  read_ngauge_conf();
  
  set_inversions();
  set_mes_gamma_contr_list();
  set_mes_contr_list();
  
  ///////////////////// finished reading apart from conf list ///////////////
  
  nind=nleptons*nweak_ind*norie*nr*nins;
  meslep_hadr_part=nissa_malloc("hadr",loc_vol,spinspin);
  meslep_contr=nissa_malloc("meslep_contr",glb_size[0]*nweak_ind*nmeslep_proj*nind,complex);
  allocate_source();
  allocate_photon_fields();
  allocate_mes_contr();
  allocate_Q_prop();
  allocate_L_prop();
  temp_lep=nissa_malloc("temp_lep",loc_vol+bord_vol,spinspin);
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
}

//handle to discard the source
void skip_conf()
{
  for(int isource=0;isource<nsources;isource++)
    {
      coords coord;
      generate_random_coord(coord);
      generate_stochastic_tlSym_gauge_propagator_source(photon_eta);
      generate_original_source();
    }
}

//init a new conf
void start_new_conf()
{
  setup_conf(conf,old_theta,put_theta,conf_path,rnd_gauge_transform,free_theory);
  
  //reset contractions
  vector_reset(mes_contr);
  vector_reset(meslep_contr);
}

/////////////////////////////////////////////// mesoleptonic contractions //////////////////////////////////////////

//compute the meson part of the lepton contraction function
//as usual, FIRST propagator is reverted
THREADABLE_FUNCTION_8ARG(meson_part_leptonic_contr, spinspin*,hadr, int,iq1, int,prop1_type, int,r1, int,iq2, int,prop2_type, int,r2, int,irev)
{
  GET_THREAD_ID();
  
  vector_reset(hadr);
  
  //it's just the matter of inserting gamma5*gamma5=identity between S1^dag and S2
  //on the sink gamma5 must still be inserted!
    for(int id_so=0;id_so<NDIRAC;id_so++)
      for(int ic_so=0;ic_so<nso_col;ic_so++)
	{
	  int ip1=iqprop(iq1,prop1_type,r1,id_so,ic_so);
	  int ip2=iqprop(iq2,prop2_type,r2,id_so,ic_so);
	  if(irev==1) std::swap(ip1,ip2); //select the propagator to revert
	  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	    for(int ic_si=0;ic_si<nso_spi;ic_si++)
	      for(int id_si1=0;id_si1<NDIRAC;id_si1++)
		for(int id_si2=0;id_si2<NDIRAC;id_si2++)
		  //this way when taking the trace with dirac matrix, that is acting on S2, as it should
		  complex_summ_the_conj1_prod(hadr[ivol][id_si2][id_si1],Q[ip1][ivol][id_si1][ic_si],Q[ip2][ivol][id_si2][ic_si]);
	}
  THREAD_BARRIER();
}
THREADABLE_FUNCTION_END

//compute the leptonic part of the contraction function
THREADABLE_FUNCTION_6ARG(attach_leptonic_contr, spinspin*,hadr, int,iprop, int,ilepton, int,orie, int,rl, int,ext_ind)
{
  GET_THREAD_ID();
  
  complex *loc_contr=nissa_malloc("loc_contr",glb_size[0],complex);
  vector_reset(loc_contr);
  
  //get the lepton info and prop
  tm_quark_info le=get_lepton_info(ilepton,orie,rl);
  spinspin *lept=L[iprop];
  
  //get the projectors
  spinspin promu[2],pronu[2];
  twisted_on_shell_operator_of_imom(promu[0],le,0,false,-1,base);
  if(follow_chris_or_nazario==follow_nazario) twisted_on_shell_operator_of_imom(promu[1],le,0,false,+1,base);
  else twisted_on_shell_operator_of_imom(promu[1],le,0,false,-1,base);
  naive_massless_on_shell_operator_of_imom(pronu[0],le.bc,0,-1);
  if(follow_chris_or_nazario==follow_nazario) naive_massless_on_shell_operator_of_imom(pronu[1],le.bc,0,+1);
  else naive_massless_on_shell_operator_of_imom(pronu[1],le.bc,0,-1);
  if(follow_chris_or_nazario==follow_chris)
    for(int i=0;i<2;i++)
      safe_spinspin_prod_dirac(promu[i],promu[i],base_gamma+map_mu[0]);
  
  //compute the right part of the leptonic loop: G0 G^dag
  dirac_matr meslep_proj_gamma[nmeslep_proj];
  for(int ig_proj=0;ig_proj<nmeslep_proj;ig_proj++)
    {
      int ig=meslep_projs[ig_proj];
      dirac_matr temp_gamma;
      dirac_herm(&temp_gamma,base_gamma+ig);
      dirac_prod(meslep_proj_gamma+ig_proj,base_gamma+map_mu[0],&temp_gamma);
    }
  //insert gamma5 on the sink-hadron-gamma: S1^dag G5 GW S2 (G5 G5) - no dagger, no commutator because it's on the LO leptonic part
  dirac_matr weak_ins_hadr_gamma[nweak_ins];
  for(int ins=0;ins<nweak_ins;ins++) dirac_prod(weak_ins_hadr_gamma+ins,base_gamma+5,base_gamma+list_weak_insq[ins]);
  
  //define the combined weak projectors (see below)
  dirac_matr neutr_1m_g5_proj;
  dirac_subt(&neutr_1m_g5_proj,base_gamma+0,base_gamma+5);
  
  for(int ins=0;ins<nweak_ins;ins++)
    {
      //define a local storage
      spinspin mesolep_loc_contr[loc_size[0]];
      for(int i=0;i<loc_size[0];i++) spinspin_put_to_zero(mesolep_loc_contr[i]);
      
      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	{
	  int t=loc_coord_of_loclx[ivol][0];
	  
	  //multiply lepton side on the right (source) side
	  spinspin la;
	  unsafe_spinspin_prod_dirac(la,lept[ivol],base_gamma+list_weak_insl[ins]);
	  
	  //include 4*(1-5)/2/2=(1-5) coming from the two neturino projector+(1-g5) weak lepton structure
	  //the second /2 comes from sqr(1/sqrt(2)) of 1502.00257
	  spinspin l;
	  unsafe_spinspin_prod_dirac(l,la,&neutr_1m_g5_proj);
	  
	  //get the neutrino phase (multiply hadron side) - notice that the sign of momentum is internally reversed
	  complex ph;
	  get_antineutrino_source_phase_factor(ph,ivol,ilepton,le.bc);
	  
	  //trace hadron side
	  complex h;
	  trace_spinspin_with_dirac(h,hadr[ivol],weak_ins_hadr_gamma+ins);
	  
	  //combine mesolep
	  complex_prodassign(h,ph);
	  spinspin_summ_the_complex_prod(mesolep_loc_contr[t],l,h);
	}
      glb_threads_reduce_double_vect((double*)mesolep_loc_contr,loc_size[0]*sizeof(spinspin)/sizeof(double));
      
      //save projection on LO
      for(int ig_proj=0;ig_proj<nmeslep_proj;ig_proj++)
	NISSA_PARALLEL_LOOP(loc_t,0,loc_size[0])
	  {
	    int glb_t=loc_t+rank_coord[0]*loc_size[0];
	    int ilnp=(glb_t>=glb_size[0]/2); //select the lepton/neutrino projector
	    
	    spinspin td;
	    unsafe_spinspin_prod_spinspin(td,mesolep_loc_contr[loc_t],pronu[ilnp]);
	    spinspin dtd;
	    unsafe_spinspin_prod_spinspin(dtd,promu[ilnp],td);
	    complex mesolep;
	    trace_spinspin_with_dirac(mesolep,dtd,meslep_proj_gamma+ig_proj);
	    
	    //summ the average
	    int i=glb_t+glb_size[0]*(ig_proj+nmeslep_proj*(list_weak_ind_contr[ins]+nweak_ind*ext_ind));
	    complex_summ_the_prod_double(meslep_contr[i],mesolep,1.0/glb_spat_vol); //here to remove the statistical average on xw
	  }
      if(IS_MASTER_THREAD) nmeslep_contr+=nmeslep_proj;
      THREAD_BARRIER();
    }
  
  nissa_free(loc_contr);
}
THREADABLE_FUNCTION_END

//return the index of the combination of r, orientation, etc
int meslep_contrpack_ind(int rl,int orie,int r2,int irev,int qins,int ilepton)
{return rl+nr*(orie+norie*(r2+nr*(irev+nrev*(qins+nins*ilepton))));}

//compute the total meslep contraction functions
void compute_meslep_contr()
{
  master_printf("Computing leptonic contraction functions\n");
  meslep_contr_time-=take_time();
  
 for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int qins=0;qins<nins;qins++)
      for(int irev=0;irev<nrev;irev++)
	for(int r2=0;r2<nr;r2++)
	  {
	    //takes the index of the quarks
	    int iq1=lep_contr_iq1[ilepton];
	    int iq2=lep_contr_iq2[ilepton];
	    
	    //takes the propagators
	    int PROP1_TYPE=PROP_0,PROP2_TYPE=PROP_0;
	    if(qins==1) PROP1_TYPE=PROP_PHOTON_A;
	    if(qins==2) PROP2_TYPE=PROP_PHOTON_A;
	    
	    //compute the meson part
	    meson_part_leptonic_contr(meslep_hadr_part, iq1,PROP1_TYPE,r2, iq2,PROP2_TYPE,r2, irev);
	    
	    for(int orie=0;orie<norie;orie++)
	      for(int rl=0;rl<nr;rl++)
		{
		  //contract with lepton
		  int ilins=(qins!=0);
		  int iprop=ilprop(ilepton,ilins,orie,rl);
		  int ind=meslep_contrpack_ind(rl,orie,r2,irev,qins,ilepton);
		  attach_leptonic_contr(meslep_hadr_part,iprop,ilepton,orie,rl,ind);
		}
	  }
  
  meslep_contr_time+=take_time();
}

//print out contractions
void print_meslep_contr()
{
  contr_print_time-=take_time();
  
  //open file and reduce
  FILE *fout=open_file(combine("%s/mesolep_contr",outfolder).c_str(),"w");
  glb_nodes_reduce_complex_vect(meslep_contr,glb_size[0]*nweak_ind*nmeslep_proj*nind);
  
  //write down
  for(int ilepton=0;ilepton<nleptons;ilepton++)
    for(int qins=0;qins<nins;qins++)
      for(int irev=0;irev<nrev;irev++)
	for(int r2=0;r2<nr;r2++)
	  for(int orie=0;orie<norie;orie++)
	    for(int rl=0;rl<nr;rl++)
	      {
		int contrpack_ind=meslep_contrpack_ind(rl,orie,r2,irev,qins,ilepton);
		
		if(!pure_wilson) master_fprintf(fout," # mlept[%d]=%lg mq1=%lg mq2=%lg qins=%d qrev=%d rq1=%d rq2=%d lep_orie=%+d rl=%d\n\n",
						ilepton,leps[ilepton].mass,qmass[lep_contr_iq1[ilepton]],qmass[lep_contr_iq2[ilepton]],qins,irev+1,!r2,r2,sign_orie[orie],rl);
		else             master_fprintf(fout," # klept[%d]=%lg kappaq1=%lg kappaq2=%lg qins=%d qrev=%d lep_orie=%+d\n\n",
						ilepton,leps[ilepton].kappa,qkappa[lep_contr_iq1[ilepton]],qkappa[lep_contr_iq2[ilepton]],qins,irev+1,sign_orie[orie]);
		for(int ind=0;ind<nweak_ind;ind++)
		  for(int ig_proj=0;ig_proj<nmeslep_proj;ig_proj++)
		    {
		      master_fprintf(fout," # qins=%s lins=%s proj=%s\n\n",list_weak_ind_nameq[ind],list_weak_ind_namel[ind],gtag[meslep_projs[ig_proj]]);
		      for(int t=0;t<glb_size[0];t++)
			{
			  int i=t+glb_size[0]*(ig_proj+nmeslep_proj*(ind+nweak_ind*contrpack_ind));
			  master_fprintf(fout,"%+016.16lg %+016.16lg\n",meslep_contr[i][RE]/nsources,meslep_contr[i][IM]/nsources);
			}
		      master_fprintf(fout,"\n");
		    }
	      }
  close_file(fout);
  
  contr_print_time+=take_time();
}

//close deallocating everything
void close()
{
  print_statistics();
  
  free_photon_fields();
  free_source();
  free_Q_prop();
  free_L_prop();
  nissa_free(conf);
  free_mes_contr();
  nissa_free(meslep_hadr_part);
  nissa_free(meslep_contr);
  nissa_free(lep_contr_iq1);
  nissa_free(lep_contr_iq2);
  nissa_free(leps);
  nissa_free(lep_energy);
  nissa_free(neu_energy);
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
  int iconf=0;
  while(read_conf_parameters(iconf,skip_conf,finish_file_present))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      for(int isource=0;isource<nsources;isource++)
	{
	  master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
	  
	  //init
	  random_shift_gauge_conf(conf,old_theta,put_theta);
	  generate_photon_stochastic_propagator();
	  generate_original_source();
	  
	  generate_lepton_propagators();
	  generate_quark_propagators(isource);
	  
	  compute_meslep_contr();
	  compute_mes_contr();
	}
      
      //print out contractions
      print_meslep_contr();
      print_mes_contr();
      
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
