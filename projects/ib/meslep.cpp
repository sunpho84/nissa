#include <nissa.hpp>

#define EXTERN_MESLEP
 #include "meslep.hpp"

#include "prop.hpp"

namespace nissa
{
  //read all the parameters to contract with leptons
  void read_meslep_contr_pars()
  {
    //Leptons
    if(twisted_run) read_str_int("NMesLepQ1Q2LepmassMesMass",&nquark_lep_combos);
    else            read_str_int("NMesLepQ1Q2LepkappaMesMass",&nquark_lep_combos);
    
    if(nquark_lep_combos)  read_gospel_convention();
    
    lep_contr_iq1=nissa_malloc("lep_contr_iq1",nquark_lep_combos,int);
    lep_contr_iq2=nissa_malloc("lep_contr_iq2",nquark_lep_combos,int);
    leps=nissa_malloc("leps",nquark_lep_combos,tm_quark_info);
    lep_energy=nissa_malloc("lep_energy",nquark_lep_combos,double);
    neu_energy=nissa_malloc("neu_energy",nquark_lep_combos,double);
    for(int il=0;il<nquark_lep_combos;il++)
      {
	//read quarks identfiying the mesons
	read_int(lep_contr_iq1+il);
	read_int(lep_contr_iq2+il);
	
	//if not pure wilson read mass
	if(not twisted_run) leps[il].mass=0;
	else                read_double(&leps[il].mass);
	
	//antiperiodic or periodic
	leps[il].bc[0]=temporal_bc;
	
	//maximal twist (if tm), otherwise read kappa
	if(not twisted_run) read_double(&leps[il].kappa);
	else                leps[il].kappa=0.125;
	leps[il].r=0;
	
	//read the mass of the meson (that must have been determined outside)
	double mes_mass;
	read_double(&mes_mass);
	
	//set initial value of bc and check kinematic
	for(int i=1;i<NDIM;i++) leps[il].bc[i]=0;
	if(tm_quark_energy(leps[il],0)>=mes_mass) crash("initial state is lighter (%lg) than final state at rest (%lg)!",mes_mass,tm_quark_energy(leps[il],0));
	
	//compute meson momentum and bc
	double err;
	master_printf("Resolving kinematical condition for combination of quarks %d/%d\n",il+1,nquark_lep_combos);
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
	    
	    master_printf("  lep_e: %+10.10lg, neu_e: %+10.10lg, mes_mass: %lg, error: %lg, der: %lg\n",lep_energy,neu_energy,mes_mass,err,der);
	  }
	while(fabs(err)>1e-14);
	
	//write down energy
	lep_energy[il]=tm_quark_energy(leps[il],0);
	neu_energy[il]=naive_massless_quark_energy(leps[il].bc,0);
	master_printf(" ilepton %d, lepton energy: %lg, neutrino energy: %lg\n",il,lep_energy[il],neu_energy[il]);
	master_printf(" lep+neut energy: %lg\n",lep_energy[il]+neu_energy[il]);
	master_printf(" bc: %+16.16lg\n\n",leps[il].bc[1]);
      }
  }
  
  //allocate all leptonic propagators
  void allocate_L_prop()
  {
    // temp_lep=nissa_malloc("temp_lep",loc_vol+bord_vol,spinspin);
    // nlprop=ilprop(nquark_lep_combos-1,nlins-1,norie-1,nr_lep-1)+1;
    // L=nissa_malloc("L*",nlprop,spinspin*);
    // for(int iprop=0;iprop<nlprop;iprop++) L[iprop]=nissa_malloc("L",loc_vol+bord_vol,spinspin);
  }
  
  //free all leptonic propagators
  void free_L_prop()
  {
    // for(int iprop=0;iprop<nlprop;iprop++) nissa_free(L[iprop]);
    // nissa_free(L);
    // nissa_free(temp_lep);
  }
  
  //return appropriately modified info
  tm_quark_info get_lepton_info(int ilepton,int orie,int r)
  {
    tm_quark_info le=leps[ilepton];
    le.r=r;
    le.bc[0]=temporal_bc;
    for(int i=1;i<NDIM;i++) le.bc[i]*=sign_orie[orie];
    
    return le;
  }
  
  //compute phase exponent for space part: vec{p}*\vec{x}
  CUDA_HOST_AND_DEVICE double get_space_arg(int ivol,const momentum_t& bc)
  {
    double arg=0;
    for(int mu=1;mu<NDIM;mu++)
      {
	double step=bc[mu]*M_PI/glbSize[mu];
	arg+=step*rel_coord_of_loclx(ivol,mu);
      }
    return arg;
  }
  
  //compute the phase for lepton on its sink
  CUDA_HOST_AND_DEVICE void get_lepton_sink_phase_factor(complex out,int ivol,int ilepton,tm_quark_info le)
  {
    //compute space and time factor
    double arg=get_space_arg(ivol,le.bc);
    int t=rel_time_of_loclx(ivol);
    if(follow_chris_or_nazario==follow_nazario and t>=glbSize[0]/2) t=glbSize[0]-t;
    double ext=exp(t*lep_energy[ilepton]);
    
    //compute full exponential (notice the factor -1)
    out[RE]=cos(-arg)*ext;
    out[IM]=sin(-arg)*ext;
  }
  
  //compute the phase for antineutrino - the orientation is that of the muon (as above)
  CUDA_HOST_AND_DEVICE void get_antineutrino_source_phase_factor(complex out,int ivol,int ilepton,const momentum_t& bc)
  {
    //compute space and time factor
    double arg=get_space_arg(ivol,bc);
    int t=rel_time_of_loclx(ivol);
    if(follow_chris_or_nazario==follow_nazario and t>=glbSize[0]/2) t=glbSize[0]-t;
    double ext=exp(t*neu_energy[ilepton]);
    
    //compute full exponential (notice the factor +1)
    out[RE]=cos(+arg)*ext;
    out[IM]=sin(+arg)*ext;
  }
  
  //set everything to a phase factor
  void set_to_lepton_sink_phase_factor(spinspin *prop,int ilepton,tm_quark_info &le)
  {
    
    vector_reset(prop);
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	complex ph;
	get_lepton_sink_phase_factor(ph,ivol,ilepton,le);
	spinspin_put_to_diag(prop[ivol],ph);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(prop);
  }
  
  //insert the photon on the source side
  void insert_photon_on_the_source(spinspin* prop,spin1field* A,int* dirs,tm_quark_info le,int twall)
  {
    
    //select A
    communicate_lx_spin1field_borders(A);
    
    //copy on the temporary and communicate borders
    vector_copy(temp_lep,prop);
    communicate_lx_spinspin_borders(temp_lep);
    vector_reset(prop);
    
    if(!loc_muon_curr)
      {
	dirac_matr GAMMA;
	if(twisted_run>0) GAMMA=dirac_prod_double(base_gamma[0],1);
	else GAMMA=dirac_prod_idouble(base_gamma[5],-tau3[le.r]);
	
	//prepare each propagator for a single lepton
	//by computing i(phi(x-mu)A_mu(x-mu)(-i t3 g5-gmu)/2-phi(x+mu)A_mu(x)(-i t3 g5+gmu)/2)=
	//(ph0 A_mu(x-mu)g[r][0][mu]-ph0 A_mu(x)g[r][1][mu])=
	for(int mu=0;mu<NDIM;mu++)
	  if(dirs[mu])
	    NISSA_PARALLEL_LOOP(ivol,0,locVol)
	      if(twall==-1 or rel_time_of_loclx(ivol)==twall)
		{
		  //phases
		  complex phase;
		  phase[0]=cos(le.bc[mu]*M_PI);
		  phase[1]=sin(le.bc[mu]*M_PI);
		  
		  //find neighbors
		  int ifw=loclxNeighup[ivol][mu];
		  int ibw=loclxNeighdw[ivol][mu];
		  
		  //compute phase factor
		  spinspin ph_bw,ph_fw;
		  
		  //transport down and up
		  if(rel_coord_of_loclx(ivol,mu)==glbSize[mu]-1) unsafe_spinspin_prod_complex_conj2(ph_fw,temp_lep[ifw],phase);
		  else spinspin_copy(ph_fw,temp_lep[ifw]);
		  if(rel_coord_of_loclx(ivol,mu)==0) unsafe_spinspin_prod_complex(ph_bw,temp_lep[ibw],phase);
		  else spinspin_copy(ph_bw,temp_lep[ibw]);
		  
		  //fix coefficients, i is inserted here!
		  //also dir selection is made here
		  spinspin_prodassign_idouble(ph_fw,+0.5*dirs[mu]);
		  spinspin_prodassign_idouble(ph_bw,-0.5*dirs[mu]);
		  
		  //fix insertion of the current
		  safe_spinspin_prod_complex(ph_fw,ph_fw,A[ivol][mu]);
		  safe_spinspin_prod_complex(ph_bw,ph_bw,A[ibw][mu]);
		  
		  //summ and subtract the two
		  spinspin fw_M_bw,fw_P_bw;
		  spinspin_subt(fw_M_bw,ph_fw,ph_bw);
		  spinspin_summ(fw_P_bw,ph_fw,ph_bw);
		  
		  //put GAMMA on the summ
		  spinspin temp_P;
		  unsafe_spinspin_prod_dirac(temp_P,fw_P_bw,GAMMA);
		  spinspin_summassign(prop[ivol],temp_P);
		  
		  //put gmu on the diff
		  spinspin temp_M;
		  unsafe_spinspin_prod_dirac(temp_M,fw_M_bw,base_gamma[igamma_of_mu[mu]]);
		  spinspin_summassign(prop[ivol],temp_M);
		}
	NISSA_PARALLEL_LOOP_END;
      }
    else
      {
	for(int mu=0;mu<NDIM;mu++)
	  if(dirs[mu])
	    NISSA_PARALLEL_LOOP(ivol,0,locVol)
	      if(twall==-1 or rel_time_of_loclx(ivol)==twall)
		{
		  spinspin temp1,temp2;
		  unsafe_spinspin_prod_dirac(temp1,temp_lep[ivol],base_gamma[igamma_of_mu[mu]]);
		  unsafe_spinspin_prod_complex(temp2,temp1,A[ivol][mu]);
		  spinspin_summ_the_prod_idouble(prop[ivol],temp2,1);
		}
	NISSA_PARALLEL_LOOP_END;
      }
    
    set_borders_invalid(prop);
  }
  
  //generate all the lepton propagators, pointing outward
  //the computations is done by:
  // 1)putting the correct phase in x space, given by exp(E_mu*t-i*vec(p)*vec(x))
  // 2)multiplying it by the conserved current inserting eta or phi
  // 3)going to momentum space
  // 4)multiplying by the lepton propagator in momentum space
  // 5)coming back to x space
  void generate_lepton_propagators()
  {
    
    // if(IS_MASTER_THREAD) lepton_prop_time-=take_time();
    
    // for(int ilepton=0;ilepton<nquark_lep_combos;ilepton++)
    //   for(int ilins=0;ilins<nlins;ilins++)
    // 	for(int ori=0;ori<norie;ori++)
    // 	  for(int r=0;r<nr_lep;r++)
    // 	    {
    // 	      //set the properties of the meson
    // 	      //time boundaries are anti-periodic, space are as for external line
    // 	      tm_quark_info le=get_lepton_info(ilepton,ori,r);
	      
    // 	      //select the propagator
    // 	      int iprop=ilprop(ilepton,ilins,ori,r);
    // 	      spinspin *prop=L[iprop];
	      
    // 	      //put it to a phase
    // 	      set_to_lepton_sink_phase_factor(prop,ilepton,le);
	      
    // 	      //if we are doing Nazario's way (with the external line) add it
    // 	      if(follow_chris_or_nazario==follow_nazario)
    // 		{
    // 		  //select only the wall
    // 		  int tmiddle=glb_size[0]/2;
    // 		  select_propagator_timeslice(prop,prop,tmiddle);
    // 		  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le,base);
    // 		}
	      
    // 	      //insert or not photon
    // 	      if(ilins)
    // 		{
    // 		  //insert photon and prolong
    // 		  insert_photon_on_the_source(prop,photon_field,all_dirs,le,-1); //all times
    // 		  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le,base);
    // 		}
    // 	    }
    
    // if(IS_MASTER_THREAD) lepton_prop_time+=take_time();
  }
  
  /*
    the loop is normalised such that the physical rate at leading order
    is obtained multiplying the loop by Gf^2 fpi^2 * phi2 (space phase
    factor) which is (1-rl^2)/(16 pi mpi) where rl=ml/mpi, whereas the
    interference is obtained by the full mesolepton contrelation
    multiplied by 4*mpi*fpi*Gf^2*phi2 */
  
  //compute the meson part of the lepton contraction function
  //as usual, FIRST propagator is reverted
  void meson_part_leptonic_contr(spinspin* hadr,int iq1,int prop1_type,int iq2,int prop2_type,int irev)
  {
    
    // vector_reset(hadr);
    
    // //it's just the matter of inserting gamma5*gamma5=identity between S1^dag and S2
    // //on the sink gamma5 must still be inserted!
    // for(int id_so=0;id_so<nso_spi;id_so++)
    //   for(int ic_so=0;ic_so<nso_col;ic_so++)
    // 	{
    // 	  int isou=so_sp_col_ind(id_so,ic_so);
    // 	  if(irev==1) std::swap(iq1,iq2); //select the propagator to revert
    // 	  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    // 	    for(int ic_si=0;ic_si<NCOL;ic_si++)
    // 	      for(int id_si1=0;id_si1<NDIRAC;id_si1++)
    // 		for(int id_si2=0;id_si2<NDIRAC;id_si2++)
    // 		  //this way when taking the trace with dirac matrix, that is acting on S2, as it should
    // 		  complex_summ_the_conj1_prod(hadr[ivol][id_si2][id_si1],Q[iq1][isou][ivol][id_si1][ic_si],Q[iq2][isou][ivol][id_si2][ic_si]);
    //NISSA_PARALLEL_LOOP_END;
    // 	}
    // THREAD_BARRIER();
  }
  
  //compute the leptonic part of the contraction function
  void attach_leptonic_contr(spinspin* hadr,int iprop,int ilepton,int orie,int rl,int ext_ind)
  {
    
    complex *loc_contr=nissa_malloc("loc_contr",glbSize[0],complex);
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
	safe_spinspin_prod_dirac(promu[i],promu[i],base_gamma[igamma_of_mu[0]]);
    
    //compute the right part of the leptonic loop: G0 G^dag
    dirac_matr meslep_proj_gamma[nmeslep_proj];
    for(int ig_proj=0;ig_proj<nmeslep_proj;ig_proj++)
      {
	int ig=meslep_projs[ig_proj];
	const dirac_matr temp_gamma=dirac_herm(base_gamma[ig]);
	meslep_proj_gamma[ig_proj]=base_gamma[igamma_of_mu[0]]*temp_gamma;
      }
    //insert gamma5 on the sink-hadron-gamma: S1^dag G5 GW S2 (G5 G5) - no dagger, no commutator because it's on the LO leptonic part
    dirac_matr weak_ins_hadr_gamma[nmeslep_weak_ins];
    for(int ins=0;ins<nmeslep_weak_ins;ins++)
      weak_ins_hadr_gamma[ins]=base_gamma[5]*base_gamma[list_weak_insq[ins]];
    
    //define the combined weak projectors (see below)
    dirac_matr neutr_1m_g5_proj=dirac_subt(base_gamma[0],base_gamma[5]);
    
    for(int ins=0;ins<nmeslep_weak_ins;ins++)
      {
	//define a local storage
	spinspin mesolep_loc_contr[locSize[0]];
	for(int i=0;i<locSize[0];i++) spinspin_put_to_zero(mesolep_loc_contr[i]);
	
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  {
	    [[maybe_unused]] int t=locCoordOfLoclx[ivol][0];
	    
	    //multiply lepton side on the right (source) side
	    spinspin la;
	    unsafe_spinspin_prod_dirac(la,lept[ivol],base_gamma[list_weak_insl[ins]]);
	    
	    //include 4*(1-5)/2/2=(1-5) coming from the two neturino projector+(1-g5) weak lepton structure
	    //the second /2 comes from sqr(1/sqrt(2)) of 1502.00257
	    spinspin l;
	    unsafe_spinspin_prod_dirac(l,la,neutr_1m_g5_proj);
	    
	    //get the neutrino phase (multiply hadron side) - notice that the sign of momentum is internally reversed
	    complex ph;
	    get_antineutrino_source_phase_factor(ph,ivol,ilepton,le.bc);
	    
	    //trace hadron side
	    complex h;
	    trace_spinspin_with_dirac(h,hadr[ivol],weak_ins_hadr_gamma[ins]);
	    
	    //combine mesolep
	    complex_prodassign(h,ph);
	    spinspin_summ_the_complex_prod(mesolep_loc_contr[t],l,h);
	  }
	NISSA_PARALLEL_LOOP_END;
	crash("#warning glb_threads_reduce_double_vect((double*)mesolep_loc_contr,loc_size[0]*sizeof(spinspin)/sizeof(double));");
	
	//save projection on LO
	for(int ig_proj=0;ig_proj<nmeslep_proj;ig_proj++)
	  NISSA_PARALLEL_LOOP(loc_t,0,locSize[0])
	    {
	      int glb_t=loc_t+rank_coord[0]*locSize[0];
	      int ilnp=(glb_t>=glbSize[0]/2); //select the lepton/neutrino projector
	      
	      spinspin td;
	      unsafe_spinspin_prod_spinspin(td,mesolep_loc_contr[loc_t],pronu[ilnp]);
	      spinspin dtd;
	      unsafe_spinspin_prod_spinspin(dtd,promu[ilnp],td);
	      complex mesolep;
	      trace_spinspin_with_dirac(mesolep,dtd,meslep_proj_gamma[ig_proj]);
	      
	      //summ the average
	      int i=glb_t+glbSize[0]*(ig_proj+nmeslep_proj*(list_weak_ind_contr[ins]+nindep_meslep_weak*ext_ind));
	      complex_summ_the_prod_double(meslep_contr[i],mesolep,1.0/glbSpatVol); //here to remove the statistical average on xw
	    }
	NISSA_PARALLEL_LOOP_END;
	if(IS_MASTER_THREAD) nmeslep_contr_made+=nmeslep_proj;
	THREAD_BARRIER();
      }
    
    nissa_free(loc_contr);
  }
    
  //compute the total meslep contraction functions
  void compute_meslep_contr()
  {
    crash("to be reviewed");
    // master_printf("Computing leptonic contraction functions\n");
    // meslep_contr_time-=take_time();
    
    // for(int ilepton=0;ilepton<nquark_lep_combos;ilepton++)
    //   for(int qins=0;qins<nins;qins++)
    // 	for(int irev=0;irev<nrev;irev++)
    // 	  {
    // 	    //takes the index of the quarks
    // 	    int iq1=lep_contr_iq1[ilepton];
    // 	    int iq2=lep_contr_iq2[ilepton];
	    
    // 	    //takes the propagators
    // 	    int prop1_type=PROP_0,prop2_type=PROP_0;
    // 	    if(qins==1) prop1_type=PROP_PHOTON_A;
    // 	    if(qins==2) prop2_type=PROP_PHOTON_A;
	    
    // 	    //compute the meson part
    // 	    meson_part_leptonic_contr(meslep_hadr_part, iq1,prop1_type, iq2,prop2_type, irev);
	    
    // 	    for(int orie=0;orie<norie;orie++)
    // 	      for(int rl=0;rl<nr_lep;rl++)
    // 		{
    // 		    //contract with lepton
    // 		    int ilins=(qins!=0);
    // 		    int iprop=ilprop(ilepton,ilins,orie,rl);
    // 		    int ind=meslep_contrpack_ind(rl,orie,irev,qins,ilepton);
    // 		    attach_leptonic_contr(meslep_hadr_part,iprop,ilepton,orie,rl,ind);
    // 		  }
    // 	    }
    
    // meslep_contr_time+=take_time();
  }
  
  //print out contractions
  void print_meslep_contr()
  {
    crash("to be reviewd");
    // contr_print_time-=take_time();
    
    // //open file and reduce
    // FILE *fout=open_file(combine("%s/mesolep_contr",outfolder).c_str(),"w");
    // glb_nodes_reduce_complex_vect(meslep_contr,glb_size[0]*nindep_meslep_weak*nmeslep_proj*nmeslep_corr);
    
    // //write down
    // for(int ilepton=0;ilepton<nquark_lep_combos;ilepton++)
    //   for(int qins=0;qins<nins;qins++)
    // 	for(int irev=0;irev<nrev;irev++)
    // 	  for(int orie=0;orie<norie;orie++)
    // 	    for(int rl=0;rl<nr_lep;rl++)
    // 	      {
    // 		int contrpack_ind=meslep_contrpack_ind(rl,orie,irev,qins,ilepton);
    // 		int iq1=lep_contr_iq1[ilepton];
    // 		int iq2=lep_contr_iq2[ilepton];
		
    // 		if(twisted_run) master_fprintf(fout," # mlept[%d]=%lg mq1=%lg mq2=%lg qins=%d qrev=%d rq1=%d rq2=%d lep_orie=%+d rl=%d\n\n",
    // 					       ilepton,leps[ilepton].mass,qmass[iq1],qmass[iq2],qins,irev+1,!qr[iq1],qr[iq2],sign_orie[orie],rl);
    // 		else            master_fprintf(fout," # klept[%d]=%lg kappaq1=%lg kappaq2=%lg qins=%d qrev=%d lep_orie=%+d\n\n",
    // 					       ilepton,leps[ilepton].kappa,qkappa[iq1],qkappa[iq2],qins,irev+1,sign_orie[orie]);
    // 		for(int ind=0;ind<nindep_meslep_weak;ind++)
    // 		  for(int ig_proj=0;ig_proj<nmeslep_proj;ig_proj++)
    // 		    {
    // 		      master_fprintf(fout," # qins=%s lins=%s proj=%s\n\n",list_weak_ind_nameq[ind],list_weak_ind_namel[ind],gtag[meslep_projs[ig_proj]]);
    // 		      for(int t=0;t<glb_size[0];t++)
    // 			{
    // 			  int i=t+glb_size[0]*(ig_proj+nmeslep_proj*(ind+nindep_meslep_weak*contrpack_ind));
    // 			  master_fprintf(fout,"%+16.016lg %+16.016lg\n",meslep_contr[i][RE]/nsources,meslep_contr[i][IM]/nsources);
    // 			}
    // 		      master_fprintf(fout,"\n");
    // 		    }
    // 	      }
    // close_file(fout);
    
    // contr_print_time+=take_time();
  }
}
