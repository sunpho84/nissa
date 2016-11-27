#include <nissa.hpp>

#define EXTERN_CONTR
 #include "contr.hpp"

#include "prop.hpp"

namespace nissa
{
  //allocate mesonic contractions
  void allocate_mes2pts_contr()
  {
    mes2pts_contr_size=glb_size[0]*mes_gamma_list.size()*mes2pts_contr_map.size();
    mes2pts_contr=nissa_malloc("mes2pts_contr",mes2pts_contr_size,complex);
  }
  
  //free mesonic contractions
  void free_mes2pts_contr()
  {nissa_free(mes2pts_contr);}
  
  //compute all the meson contractions
  THREADABLE_FUNCTION_0ARG(compute_mes2pts_contr)
  {
    GET_THREAD_ID();
    
    master_printf("Computing meson 2pts_contractions\n");
    
    // Tr [ G1 G5 S1^+ G5 G2 S2 ]      G2 is on the sink
    // (G1)_{ij(i)} (G5)_{j(i)} (S1*)^{ab}_{kj(i)} (G5)_k (G2)_{kl(k)} (S2)^{ab}_{l(k)i}
    //
    // A(i)=(G1)_{ij(i)} (G5)_{j(i)}
    // B(k)=(G5)_k (G2)_{kl(k)}
    //
    // A(i) (S1*)^{ab}_{kj(i)} B(k) (S2)^{ab}_{l(k)i}
    
    //allocate loc storage
    complex *loc_contr=new complex[mes2pts_contr_size];
    memset(loc_contr,0,sizeof(complex)*mes2pts_contr_size);
    
    if(IS_MASTER_THREAD) mes2pts_contr_time-=take_time();
    
    for(size_t icombo=0;icombo<mes2pts_contr_map.size();icombo++)
      {
	qprop_t &Q1=Q[mes2pts_contr_map[icombo].a];
	qprop_t &Q2=Q[mes2pts_contr_map[icombo].b];
	double norm=12/sqrt(Q1.ori_source_norm2*Q2.ori_source_norm2); //12 in case of a point source
	for(size_t ihadr_contr=0;ihadr_contr<mes_gamma_list.size();ihadr_contr++)
	  {
	    int ig1=mes_gamma_list[ihadr_contr].so;
	    int ig2=mes_gamma_list[ihadr_contr].si;
	    if(nso_spi==1 && ig1!=5) crash("implemented only g5 contraction on the source for non-diluted source");
	    
	    for(int i=0;i<nso_spi;i++)
	      {
		int j=(base_gamma+ig1)->pos[i];
		
		complex A;
		unsafe_complex_prod(A,(base_gamma+ig1)->entr[i],(base_gamma+5)->entr[j]);
		
		for(int b=0;b<nso_col;b++)
		  {
		    int isou=so_sp_col_ind(i,b);
		    spincolor *q1=Q1[isou];
		    spincolor *q2=Q2[isou];
		    
		    for(int k=0;k<NDIRAC;k++)
		      {
			int l=(base_gamma+ig2)->pos[k];
			
			//compute AB*norm
			complex B;
			unsafe_complex_prod(B,(base_gamma+5)->entr[k],(base_gamma+ig2)->entr[l]);
			complex AB;
			unsafe_complex_prod(AB,A,B);
			complex_prodassign_double(AB,norm);
			
			NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
			  {
			    complex c={0,0};
			    for(int a=0;a<NCOL;a++)
			      complex_summ_the_conj1_prod(c,q1[ivol][k][a],q2[ivol][l][a]);
			    complex_summ_the_prod(loc_contr[ind_mes2pts_contr(icombo,ihadr_contr,glb_coord_of_loclx[ivol][0])],c,AB);
			  }
		      }
		  }
	      }
	  }
      }
    THREAD_BARRIER();
    
    //reduce between threads and summ
    complex *red_contr=glb_threads_reduce_complex_vect(loc_contr,mes2pts_contr_size);
    NISSA_PARALLEL_LOOP(i,0,mes2pts_contr_size) complex_summassign(mes2pts_contr[i],red_contr[i]);
    //disallocate after all threads finished
    THREAD_BARRIER();
    delete[] loc_contr;
    
    //stats
    if(IS_MASTER_THREAD)
      {
	nmes2pts_contr_made+=mes2pts_contr_map.size()*mes_gamma_list.size();
	mes2pts_contr_time+=take_time();
      }
  }
  THREADABLE_FUNCTION_END
  
  //print all mesonic 2pts contractions
  void print_mes2pts_contr()
  {
    contr_print_time-=take_time();
    
    //reduce and normalise
    glb_nodes_reduce_complex_vect(mes2pts_contr,mes2pts_contr_size);
    for(int i=0;i<mes2pts_contr_size;i++) complex_prodassign_double(mes2pts_contr[i],1.0);
    
    for(size_t icombo=0;icombo<mes2pts_contr_map.size();icombo++)
      {
	FILE *fout=open_file(combine("%s/mes_contr_%s",outfolder,mes2pts_contr_map[icombo].name.c_str()),"w");
	
	print_contractions_to_file(fout,mes_gamma_list,mes2pts_contr+ind_mes2pts_contr(icombo,0,0),0,"",1.0/nhits);
	master_fprintf(fout,"\n");
	
	//close the file
	close_file(fout);
      }
    
    contr_print_time+=take_time();
  }
  
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
  
  //////////////////////////////////////// baryonic contractions //////////////////////////////////////////////////////////
  
  /*
    We follow eq.6.21 of Gattringer, pag 131 and compute the two Wick
    contractions separately, as in
    
    eps_{a,b,c} eps_{a',b',c'} (Cg5)_{al',be'} (Cg5)_{al,be}
    (P+-)_{ga,ga'} S_{be',be}{b',b} (
    S_{al',al}{a',a} S_{ga',ga}{c',c} -
    S_{al',ga}{a',c} S_{ga',al}{c',a})
    
     a',al'---------a,al           a',al'--@   @--a,al
       |             |		    |       \ /    |
     b',be'---------b,be           b',be'---------b,be
       |             |		    |       / \    |
     c',ga'---------c,ga	   c',ga'--@   @--c,ga
     
     insertions are labelled as abc on the source (left) side
   
  */
  
  //set Cg5=ig2g4g5
  void set_Cg5()
  {
    dirac_matr g2g4,C;
    dirac_prod(&g2g4,base_gamma+2,base_gamma+4);
    dirac_prod_idouble(&C,&g2g4,1);
    dirac_prod(&Cg5,&C,base_gamma+5);
  }
  
  //allocate baryionic contr
  void allocate_bar2pts_contr()
  {
    bar2pts_contr_size=ind_bar2pts_contr(bar2pts_contr_map.size()-1,2-1,glb_size[0]-1)+1;
    bar2pts_contr=nissa_malloc("bar2pts_contr",bar2pts_contr_size,complex);
  }
  
  //free them
  void free_bar2pts_contr()
  {nissa_free(bar2pts_contr);}
  
  //compute all contractions
  THREADABLE_FUNCTION_0ARG(compute_bar2pts_contr)
  {
    GET_THREAD_ID();
    master_printf("Computing baryon 2pts contractions\n");
    
    //allocate loc storage
    complex *loc_contr=new complex[bar2pts_contr_size];
    memset(loc_contr,0,sizeof(complex)*bar2pts_contr_size);
    
    const int eps[3][2]={{1,2},{2,0},{0,1}},sign[2]={1,-1};
    
    void (*list_fun[2])(complex,complex,complex)={complex_summ_the_prod,complex_subt_the_prod};
    UNPAUSE_TIMING(bar2pts_contr_time);
    for(size_t icombo=0;icombo<bar2pts_contr_map.size();icombo++)
      {
	qprop_t &Q1=Q[bar2pts_contr_map[icombo].a];
	qprop_t &Q2=Q[bar2pts_contr_map[icombo].b];
	qprop_t &Q3=Q[bar2pts_contr_map[icombo].c];
	double norm=pow(12,1.5)/sqrt(Q1.ori_source_norm2*Q2.ori_source_norm2*Q3.ori_source_norm2); //12 is even in case of a point source
	
	for(int al=0;al<NDIRAC;al++)
	  for(int ga=0;ga<NDIRAC;ga++)
	    for(int b=0;b<NCOL;b++)
	      for(int iperm=0;iperm<2;iperm++)
		{
		  int c=eps[b][iperm],a=eps[b][!iperm];
		  int be=Cg5.pos[al];
		  
		  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
		    {
		      int t=glb_coord_of_loclx[ivol][0];
		      
		      int ga1_l[2][NDIRAC]={{0,1,2,3},{2,3,0,1}}; //ga1 index for 1 or gamma0 matrix
    			  int sign_idg0[2]={(t<(glb_size[0]/2))?1:-1,-1}; //gamma0 is -1 always
    			  for(int al1=0;al1<NDIRAC;al1++)
    			    for(int b1=0;b1<NCOL;b1++)
    			      {
    				complex diquark_dir={0,0},diquark_exc={0,0};
				
    				//build the diquark
    				for(int iperm1=0;iperm1<2;iperm1++)
    				  {
    				    int c1=eps[b1][iperm1],a1=eps[b1][!iperm1];
				    
    				    for(int idg0=0;idg0<2;idg0++)
    				      {
    					int isign=((sign[iperm]*sign[iperm1]*sign_idg0[idg0])==1);
    					int ga1=ga1_l[idg0][ga];
					
    					list_fun[isign](diquark_dir,Q1[so_sp_col_ind(al,a)][ivol][al1][a1],Q3[so_sp_col_ind(ga,c)][ivol][ga1][c1]); //direct
    					list_fun[isign](diquark_exc,Q1[so_sp_col_ind(ga,c)][ivol][al1][a1],Q3[so_sp_col_ind(al,a)][ivol][ga1][c1]); //exchange
    				      }
    				  }
				
    				//close it
    				complex w;
    				unsafe_complex_prod(w,Cg5.entr[al1],Cg5.entr[al]);
    				int be1=Cg5.pos[al1];
    				complex_prodassign_double(diquark_dir,w[RE]*norm);
    				complex_prodassign_double(diquark_exc,w[RE]*norm);
    				complex_summ_the_prod(loc_contr[ind_bar2pts_contr(icombo,0,t)],Q2[so_sp_col_ind(be,b)][ivol][be1][b1],diquark_dir);
    				complex_summ_the_prod(loc_contr[ind_bar2pts_contr(icombo,1,t)],Q2[so_sp_col_ind(be,b)][ivol][be1][b1],diquark_exc);
    			      }
		    }
		}
      }
    STOP_TIMING(bar2pts_contr_time);
    
    //reduce
    complex *master_reduced_contr=glb_threads_reduce_complex_vect(loc_contr,bar2pts_contr_size);
    NISSA_PARALLEL_LOOP(i,0,bar2pts_contr_size) complex_summassign(bar2pts_contr[i],master_reduced_contr[i]);
    THREAD_BARRIER();
    delete[] loc_contr;
    
    //stats
    if(IS_MASTER_THREAD) nbar2pts_contr_made+=bar2pts_contr_map.size();
  }
  THREADABLE_FUNCTION_END
  
  //print all contractions
  void print_bar2pts_contr()
  {
    //reduce
    glb_nodes_reduce_complex_vect(bar2pts_contr,bar2pts_contr_size);
    
    for(size_t icombo=0;icombo<bar2pts_contr_map.size();icombo++)
      for(int dir_exc=0;dir_exc<2;dir_exc++)
	{
	  //open output
	  FILE *fout=open_file(combine("%s/bar_contr_%s_%s",outfolder,(dir_exc==0)?"dir":"exc",bar2pts_contr_map[icombo].name.c_str()),"w");
	  for(int t=0;t<glb_size[0];t++)
	    {
	      //remove border phase
	      double arg=3*QUARK_BOUND_COND*M_PI*t/glb_size[0];
	      complex phase={cos(arg),sin(arg)};
	      
	      //normalize for nsources and 1+g0
	      complex c;
	      complex_prod_double(c,bar2pts_contr[ind_bar2pts_contr(icombo,dir_exc,t)],1.0/(2*nhits));
	      
	      //put the phase and print
	      safe_complex_prod(c,c,phase);
	      master_fprintf(fout,"%+16.16lg %+16.16lg\n",c[RE],c[IM]);
	    }
	  
	  close_file(fout);
	}
  }
}
