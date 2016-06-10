#include <nissa.hpp>

#define EXTERN_CONTR
 #include "contr.hpp"

#include "prop.hpp"

namespace nissa
{
  //allocate mesonic contractions
  void allocate_mes_contr()
  {
    mes_contr_size=glb_size[0]*mes_gamma_list.size()*prop_mes_contr_map.size()*nqmass*nqmass*nr;
    mes_contr=nissa_malloc("mes_contr",mes_contr_size,complex);
  }
  
  //free mesonic contractions
  void free_mes_contr()
  {nissa_free(mes_contr);}
  
  //set all the mesonic contractions
  void set_mes_contr_list()
  {
    prop_mes_contr_map.clear();
    
    //non-perturbed
    prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_0));
    
    //mass corrections
    if(compute_mass_corrections) prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_S));
    
    //QED corrections
    if(compute_QED_corrections)
      {
	prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_P));
	prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_T));
	prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_PHOTON_AB));
	prop_mes_contr_map.push_back(mes_doublet_t(PROP_PHOTON_A,PROP_PHOTON_B));
      }
    //prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_VECTOR));
  }
  
  //compute all the meson contractions
  THREADABLE_FUNCTION_0ARG(compute_mes_contr)
  {
    GET_THREAD_ID();
    
    master_printf("Computing meson contractions\n");
    
    // Tr [ G1 G5 S1^+ G5 G2 S2 ]      G2 is on the sink
    // (G1)_{ij(i)} (G5)_{j(i)} (S1*)^{ab}_{kj(i)} (G5)_k (G2)_{kl(k)} (S2)^{ab}_{l(k)i}
    //
    // A(i)=(G1)_{ij(i)} (G5)_{j(i)}
    // B(k)=(G5)_k (G2)_{kl(k)}
    //
    // A(i) (S1*)^{ab}_{kj(i)} B(k) (S2)^{ab}_{l(k)i}
    
    //allocate loc storage
    complex *loc_contr=new complex[mes_contr_size];
    memset(loc_contr,0,sizeof(complex)*mes_contr_size);
    
    mes_contr_time-=take_time();
    for(size_t icombo=0;icombo<prop_mes_contr_map.size();icombo++)
      for(int imass=0;imass<nqmass;imass++)
	for(int jmass=0;jmass<nqmass;jmass++)
	  for(int r=0;r<nr;r++)
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
			int ip1=iqprop(imass,prop_mes_contr_map[icombo].a,r,j,b);
			int ip2=iqprop(jmass,prop_mes_contr_map[icombo].b,r,i,b);
			
			for(int k=0;k<NDIRAC;k++)
			  {
			    int l=(base_gamma+ig2)->pos[k];
			    
			    //compute AB
			    complex B;
			    unsafe_complex_prod(B,(base_gamma+5)->entr[k],(base_gamma+ig2)->entr[k]);
			    complex AB;
			    unsafe_complex_prod(AB,A,B);
			    
			    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
			      {
				complex c={0,0};
				for(int a=0;a<NCOL;a++)
				  complex_summ_the_conj1_prod(c,Q[ip1][ivol][k][a],Q[ip2][ivol][l][a]);
				complex_summ_the_prod(loc_contr[ind_mes_contr(icombo,imass,jmass,r,ihadr_contr,glb_coord_of_loclx[ivol][0])],c,AB);
			      }
			  }
		      }
		  }
	      }
    THREAD_BARRIER();
    
    //reduce between threads and summ
    complex *red_contr=glb_threads_reduce_complex_vect(loc_contr,mes_contr_size);
    NISSA_PARALLEL_LOOP(i,0,mes_contr_size) complex_summassign(mes_contr[i],red_contr[i]);
    //disallocate after all threads finished
    THREAD_BARRIER();
    delete[] loc_contr;
    
    //stats
    if(IS_MASTER_THREAD)
      {
	nmes_contr+=prop_mes_contr_map.size()*nqmass*nqmass*nr*mes_gamma_list.size();
	mes_contr_time+=take_time();
      }
  }
  THREADABLE_FUNCTION_END
  
  //print all contractions averaging
  void print_mes_contr()
  {
    contr_print_time-=take_time();
    
    //reduce and normalise
    double norm=1.0/nsources;
    glb_nodes_reduce_complex_vect(mes_contr,mes_contr_size);
    for(int i=0;i<mes_contr_size;i++) complex_prodassign_double(mes_contr[i],norm);
    
    int ind=0;
    for(size_t icombo=0;icombo<prop_mes_contr_map.size();icombo++)
      {
	FILE *fout=open_file(combine("%s/mes_contr_%c%c",outfolder,qprop_list[prop_mes_contr_map[icombo].a].shortname,qprop_list[prop_mes_contr_map[icombo].b].shortname).c_str(),"w");
	
	for(int imass=0;imass<nqmass;imass++)
	  for(int jmass=0;jmass<nqmass;jmass++)
	    for(int r=0;r<nr;r++)
	      {
		if(!pure_wilson) master_fprintf(fout," # m1(rev)=%lg m2(ins)=%lg r=%d\n",qmass[imass],qmass[jmass],r);
		else             master_fprintf(fout," # kappa1(rev)=%lg kappa2(ins)=%lg\n",qkappa[imass],qkappa[jmass]);
		print_contractions_to_file(fout,mes_gamma_list,mes_contr+ind*glb_size[0],0,"",1.0);
		master_fprintf(fout,"\n");
		ind+=mes_gamma_list.size();
	      }
	
	//close the file
	close_file(fout);
      }
    
    contr_print_time+=take_time();
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //set Cg5=ig2g4g5
  void set_Cg5()
  {
    dirac_matr g2g4,C;
    dirac_prod(&g2g4,base_gamma+2,base_gamma+4);
    dirac_prod_idouble(&C,&g2g4,1);
    dirac_prod(&Cg5,&C,base_gamma+5);
  }
  
  //allocate baryionic contr
  void allocate_bar_contr()
  {
    bar_contr_size=ind_bar_contr(prop_bar_contr_map.size()-1,nsm_sink-1,nqmass-1,nr-1,nqmass-1,nr-1,nqmass-1,nr-1,2-1,glb_size[0]-1)+1;
    bar_contr=nissa_malloc("bar_contr",bar_contr_size,complex);
  }
  
  //free them
  void free_bar_contr()
  {nissa_free(bar_contr);}
  
  //set all the baryonic contractions
  void set_bar_contr_list()
  {
    //clear the list
    prop_bar_contr_map.clear();
    
    //non-perturbed
    prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_0));
    
    //mass corrections
    if(compute_mass_corrections)
      {
	//scalar insertion on one of the three lines
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_S,PROP_0,PROP_0));
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_S,PROP_0));
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_S));
      }
    
    //QED corrections
    if(compute_QED_corrections)
      {
	//pseudoscalar insertion on one of the three lines
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_P,PROP_0,PROP_0));
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_P,PROP_0));
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_P));
	//tadpole insertion on one of the three lines
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_T,PROP_0,PROP_0));
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_T,PROP_0));
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_T));
	//self-energy insertion on one of the three lines
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_PHOTON_AB,PROP_0,PROP_0));
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_PHOTON_AB,PROP_0));
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_PHOTON_AB));
	//photon exchange between one of the three lines
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_PHOTON_A,PROP_PHOTON_B));
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_PHOTON_A,PROP_0,PROP_PHOTON_B));
	prop_bar_contr_map.push_back(bar_triplet_t(PROP_PHOTON_A,PROP_PHOTON_B,PROP_0));
      }
  }
  
  //compute all contractions
  THREADABLE_FUNCTION_0ARG(compute_bar_contr)
  {
    GET_THREAD_ID();
    
    master_printf("Computing baryon contractions\n");
    
    //local thread/node contractions
    complex *loc_contr=new complex[bar_contr_size];
    memset(loc_contr,0,sizeof(complex)*bar_contr_size);
    
    for(int ism_sink=0;ism_sink<nsm_sink;ism_sink++)
      {
	//smear all sinks
	START_TIMING(smear_oper_time,nsmear_oper);
	if(ism_sink)
	  for(int ip=0;ip<nqprop;ip++)
	    gaussian_smearing(Q[ip],Q[ip],ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
	STOP_TIMING(smear_oper_time);
	
	const int eps[3][2]={{1,2},{2,0},{0,1}},sign[2]={1,-1};
	
	void (*list_fun[2])(complex,complex,complex)={complex_summ_the_prod,complex_subt_the_prod};
	UNPAUSE_TIMING(bar_contr_time);
	for(size_t icombo=0;icombo<prop_bar_contr_map.size();icombo++)
	  for(int ima=0;ima<nqmass;ima++)
	    for(int imb=0;imb<nqmass;imb++)
	      for(int imc=0;imc<nqmass;imc++)
		for(int ra=0;ra<nr;ra++)
		  for(int rb=0;rb<nr;rb++)
		    for(int rc=0;rc<nr;rc++)
		      {
			if(IS_MASTER_THREAD) nbar_contr++;
			
			for(int al=0;al<NDIRAC;al++)
			  for(int ga=0;ga<NDIRAC;ga++)
			    for(int b=0;b<NCOL;b++)
			      for(int iperm=0;iperm<2;iperm++)
				{
				  int c=eps[b][iperm],a=eps[b][!iperm];
				  int be=Cg5.pos[al];
				  
				  int ipa_al_a=iqprop(ima,prop_bar_contr_map[icombo].a,ra,al,a);
				  int ipa_ga_c=iqprop(ima,prop_bar_contr_map[icombo].a,ra,ga,c);
				  int ipb=iqprop(imb,prop_bar_contr_map[icombo].b,rb,be,b);
				  int ipc_ga_c=iqprop(imc,prop_bar_contr_map[icombo].c,rc,ga,c);
				  int ipc_al_a=iqprop(imc,prop_bar_contr_map[icombo].c,rc,al,a);
				  
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
						    
						    list_fun[isign](diquark_dir,Q[ipa_al_a][ivol][al1][a1],Q[ipc_ga_c][ivol][ga1][c1]); //direct
						    list_fun[isign](diquark_exc,Q[ipa_ga_c][ivol][al1][a1],Q[ipc_al_a][ivol][ga1][c1]); //exchange
						  }
					      }
					    
					    //close it
					    complex w;
					    unsafe_complex_prod(w,Cg5.entr[al1],Cg5.entr[al]);
					    int be1=Cg5.pos[al1];
					    complex_prodassign_double(diquark_dir,w[RE]);
					    complex_prodassign_double(diquark_exc,w[RE]);
					   complex_summ_the_prod(loc_contr[ind_bar_contr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,0,t)],Q[ipb][ivol][be1][b1],diquark_dir);
					   complex_summ_the_prod(loc_contr[ind_bar_contr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,1,t)],Q[ipb][ivol][be1][b1],diquark_exc);
					  }
				    }
				}
		      }
	STOP_TIMING(bar_contr_time);
      }
    
    //reduce
    complex *master_reduced_contr=(complex*)glb_threads_reduce_double_vect((double*)loc_contr,bar_contr_size*2);
    NISSA_PARALLEL_LOOP(i,0,bar_contr_size) complex_summassign(bar_contr[i],master_reduced_contr[i]);
    THREAD_BARRIER();
    delete[] loc_contr;
  }
  THREADABLE_FUNCTION_END
  
  //print all contractions averaging
  void print_bar_contr()
  {
    //reduce
    glb_nodes_reduce_complex_vect(bar_contr,bar_contr_size);
    
    //open output
    FILE *fout=open_file(combine("%s/bar_contr",outfolder),"w");
    
    for(size_t icombo=0;icombo<prop_bar_contr_map.size();icombo++)
      for(int ism_sink=0;ism_sink<nsm_sink;ism_sink++)
	for(int ima=0;ima<nqmass;ima++)
	  for(int imb=0;imb<nqmass;imb++)
	    for(int imc=0;imc<nqmass;imc++)
	      for(int ra=0;ra<nr;ra++)
		for(int rb=0;rb<nr;rb++)
		  for(int rc=0;rc<nr;rc++)
		    for(int dir_exc=0;dir_exc<2;dir_exc++)
		      {
			master_fprintf(fout,"\n # icombo %c%c%c , sink_smeared %d , ma %lg , mb %lg , mc %lg , ra %d , rb %d , rc %d , dir_exc %d\n\n",
				       qprop_list[prop_bar_contr_map[icombo].a].shortname,
				       qprop_list[prop_bar_contr_map[icombo].b].shortname,
				       qprop_list[prop_bar_contr_map[icombo].c].shortname,
				       ism_sink,qmass[ima],qmass[imb],qmass[imc],ra,rb,rc,dir_exc);
			for(int t=0;t<glb_size[0];t++)
			  {
			    //remove border phase
			    double arg=3*QUARK_BOUND_COND*M_PI*t/glb_size[0];
			    complex phase={cos(arg),sin(arg)};
			    
			    //normalize for nsources and 1+g0
			    complex c;
			    complex_prod_double(c,bar_contr[ind_bar_contr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,dir_exc,t)],1.0/(2*nsources));
			    
			    //put the phase and print
			    safe_complex_prod(c,c,phase);
			    master_fprintf(fout,"%+016.016lg %+016.016lg\n",c[RE],c[IM]);
			  }
			master_fprintf(fout,"\n");
		      }
    
    close_file(fout);
  }
}
