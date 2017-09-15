#include <nissa.hpp>

#define EXTERN_MESLEP
 #include "meslep.hpp"

#include "prop.hpp"

namespace nissa
{
  /////////////////////////////////////////////// mesoleptonic contractions //////////////////////////////////////////
  
  /*
    the loop is normalised such that the physical rate at leading order
    is obtained multiplying the loop by Gf^2 fpi^2 * phi2 (space phase
    factor) which is (1-rl^2)/(16 pi mpi) where rl=ml/mpi, whereas the
    interference is obtained by the full mesolepton contrelation
    multiplied by 4*mpi*fpi*Gf^2*phi2 */
  
  //compute the meson part of the lepton contraction function
  //as usual, FIRST propagator is reverted
  THREADABLE_FUNCTION_6ARG(meson_part_leptonic_contr, spinspin*,hadr, int,iq1, int,prop1_type, int,iq2, int,prop2_type, int,irev)
  {
    crash("to be reviewed");
    
    // GET_THREAD_ID();
    
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
    // 	}
    // THREAD_BARRIER();
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
	safe_spinspin_prod_dirac(promu[i],promu[i],base_gamma+igamma_of_mu[0]);
    
    //compute the right part of the leptonic loop: G0 G^dag
    dirac_matr meslep_proj_gamma[nmeslep_proj];
    for(int ig_proj=0;ig_proj<nmeslep_proj;ig_proj++)
      {
	int ig=meslep_projs[ig_proj];
	dirac_matr temp_gamma;
	dirac_herm(&temp_gamma,base_gamma+ig);
	dirac_prod(meslep_proj_gamma+ig_proj,base_gamma+igamma_of_mu[0],&temp_gamma);
      }
    //insert gamma5 on the sink-hadron-gamma: S1^dag G5 GW S2 (G5 G5) - no dagger, no commutator because it's on the LO leptonic part
    dirac_matr weak_ins_hadr_gamma[nmeslep_weak_ins];
    for(int ins=0;ins<nmeslep_weak_ins;ins++) dirac_prod(weak_ins_hadr_gamma+ins,base_gamma+5,base_gamma+list_weak_insq[ins]);
    
    //define the combined weak projectors (see below)
    dirac_matr neutr_1m_g5_proj;
    dirac_subt(&neutr_1m_g5_proj,base_gamma+0,base_gamma+5);
    
    for(int ins=0;ins<nmeslep_weak_ins;ins++)
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
	      int i=glb_t+glb_size[0]*(ig_proj+nmeslep_proj*(list_weak_ind_contr[ins]+nindep_meslep_weak*ext_ind));
	      complex_summ_the_prod_double(meslep_contr[i],mesolep,1.0/glb_spat_vol); //here to remove the statistical average on xw
	    }
	if(IS_MASTER_THREAD) nmeslep_contr_made+=nmeslep_proj;
	THREAD_BARRIER();
      }
    
    nissa_free(loc_contr);
  }
  THREADABLE_FUNCTION_END
  
  //return the index of the combination of r, orientation, etc
  int meslep_contrpack_ind(int rl,int orie,int irev,int qins,int ilepton)
  {return rl+nr_lep*(orie+norie*(irev+nrev*(qins+nins*ilepton)));}
  
  //compute the total meslep contraction functions
  void compute_meslep_contr()
  {
    crash("to be revirewd");
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
