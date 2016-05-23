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
  {
    nissa_free(mes_contr);
  }
  
  //set all the mesonic contractions
  void set_mes_contr_list()
  {
    prop_mes_contr_map.clear();
    prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_0));
    prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_S));
    prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_P));
    prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_T));
    prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_PHOTON_AB));
    prop_mes_contr_map.push_back(mes_doublet_t(PROP_PHOTON_A,PROP_PHOTON_B));
    //prop_mes_contr_map.push_back(mes_doublet_t(PROP_0,PROP_VECTOR));
  }
  
  //compute all the meson contractions
  void compute_mes_contr()
  {
    contr_print_time-=take_time();
    
    master_printf("Computing meson contractions\n");
    
    complex *glb_contr=nissa_malloc("glb_contr",glb_size[0]*mes_gamma_list.size(),complex);
    complex *loc_contr=nissa_malloc("loc_contr",glb_size[0]*mes_gamma_list.size(),complex);
    
    mes_contr_time-=take_time();
    for(size_t icombo=0;icombo<prop_mes_contr_map.size();icombo++)
      for(int imass=0;imass<nqmass;imass++)
	for(int jmass=0;jmass<nqmass;jmass++)
	  for(int r=0;r<nr;r++)
	    {
	      //compute the contraction function
	      int ip1=iqprop(imass,prop_mes_contr_map[icombo].a,r);
	      int ip2=iqprop(jmass,prop_mes_contr_map[icombo].b,r);
	      
	      meson_two_points_Wilson_prop(glb_contr,loc_contr,Q[ip1],Q[ip2],mes_gamma_list);
	      nmes_contr+=mes_gamma_list.size();
	      
	      //save to the total stack
	      for(size_t ihadr_contr=0;ihadr_contr<mes_gamma_list.size();ihadr_contr++)
		for(int t=0;t<glb_size[0];t++)
		  complex_summassign(mes_contr[ind_mes_contr(icombo,imass,jmass,r,ihadr_contr,t)],glb_contr[t+glb_size[0]*ihadr_contr]);
	    }
    mes_contr_time+=take_time();
    
    nissa_free(glb_contr);
    nissa_free(loc_contr);
    
    contr_print_time+=take_time();
  }
  
  //print all contractions averaging
  void print_mes_contr()
  {
    //normalise
    double n=1.0/nsources;
    for(int i=0;i<mes_contr_size;i++) complex_prodassign_double(mes_contr[i],n);
    
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
  {
    nissa_free(bar_contr);
  }
  
  //set all the baryonic contractions
  void set_bar_contr_list()
  {
    //add the contraction combination
    prop_bar_contr_map.clear();
    prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_0));
    //scalar insertion on one of the three lines
    prop_bar_contr_map.push_back(bar_triplet_t(PROP_S,PROP_0,PROP_0));
    prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_S,PROP_0));
    prop_bar_contr_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_S));
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
  
#ifdef POINT_SOURCE_VERSION
  
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
			int ipa=iqprop(ima,prop_bar_contr_map[icombo].a,ra);
			int ipb=iqprop(imb,prop_bar_contr_map[icombo].b,rb);
			int ipc=iqprop(imc,prop_bar_contr_map[icombo].c,rc);
			
			if(IS_MASTER_THREAD) nbar_contr++;
			
			NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
			  {
			    int t=glb_coord_of_loclx[ivol][0];
			    
			    int ga1_l[2][NDIRAC]={{0,1,2,3},{2,3,0,1}}; //ga1 index for 1 or gamma0 matrix
			    int sign_idg0[2]={(t<(glb_size[0]/2))?1:-1,-1}; //gamma0 is -1 always
			    for(int al1=0;al1<NDIRAC;al1++)
			      for(int al=0;al<NDIRAC;al++)
				for(int b1=0;b1<NCOL;b1++)
				  for(int b=0;b<NCOL;b++)
				    {
				      complex diquark_dir={0,0},diquark_exc={0,0};
				      
				      //build the diquark
				      for(int iperm1=0;iperm1<2;iperm1++)
					for(int iperm=0;iperm<2;iperm++)
					  {
					    int c=eps[b][iperm],a=eps[b][!iperm];
					    int c1=eps[b1][iperm1],a1=eps[b1][!iperm1];
					    
					    for(int ga=0;ga<NDIRAC;ga++)
					      for(int idg0=0;idg0<2;idg0++)
						{
						  int isign=((sign[iperm]*sign[iperm1]*sign_idg0[idg0])==1);
						  int ga1=ga1_l[idg0][ga];
						  
						  list_fun[isign](diquark_dir,Q[ipa][ivol][a1][a][al1][al],Q[ipc][ivol][c1][c][ga1][ga]); //direct
						  list_fun[isign](diquark_exc,Q[ipa][ivol][a1][c][al1][ga],Q[ipc][ivol][c1][a][ga1][al]); //exchange
						}
					  }
				      
				      //close it
				      complex w;
				      unsafe_complex_prod(w,Cg5.entr[al1],Cg5.entr[al]);
				      int be1=Cg5.pos[al1],be=Cg5.pos[al];
				      complex_prodassign_double(diquark_dir,w[RE]);
				      complex_prodassign_double(diquark_exc,w[RE]);
				      complex_summ_the_prod(loc_contr[ind_bar_contr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,0,t)],Q[ipb][ivol][b1][b][be1][be],
							    diquark_dir);
				      complex_summ_the_prod(loc_contr[ind_bar_contr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,1,t)],Q[ipb][ivol][b1][b][be1][be],
							    diquark_exc);
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
  
#endif
  
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
