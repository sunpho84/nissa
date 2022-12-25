#include <nissa.hpp>

#define EXTERN_CONTR
# include "contr.hpp"

#include <set>

#include "prop.hpp"

namespace nissa
{
  //class to open or append a path, depending if it was already used
  class open_or_append_t
  {
    //list of already opened file
    std::set<std::string> opened;
    
  public:
    
    //open for write or append, depending
    FILE *open(const std::string &path,
	       const bool& force_append=false)
    {
      //detect the mode
      std::string mode;
      if((not force_append) and opened.find(path)==opened.end())
	{
	  mode="w";
	  opened.insert(path);
	}
      else
	mode="a";
      
      //open
      FILE *fout=
	open_file(path.c_str(),mode.c_str());
      
      return fout;
    }
  };
  
  //allocate mesonic contractions
  void allocate_mes2pts_contr()
  {
    mes2pts_contr_size=glbSize[0]*mes_gamma_list.size()*mes2pts_contr_map.size();
    mes2pts_contr=nissa_malloc("mes2pts_contr",mes2pts_contr_size,complex);
  }
  
  //free mesonic contractions
  void free_mes2pts_contr()
  {
    nissa_free(mes2pts_contr);
  }
  
  //compute a single scalar product
  void compute_prop_scalprod(complex& res,
			     const std::string& pr_dag,
			     const std::string& pr)
  {
    master_printf("Computing the scalar product between %s and %s\n",pr_dag.c_str(),pr.c_str());
    
    LxField<complex> loc("loc");
    loc.reset();
    
    for(int idc_so=0;idc_so<nso_spi*nso_col;idc_so++)
      {
	LxField<spincolor>& q_dag=Q[pr_dag][idc_so];
	LxField<spincolor>& q=Q[pr][idc_so];
	
	PAR(0,locVol,
	    CAPTURE(TO_WRITE(loc),
		    TO_READ(q_dag),
		    TO_READ(q)),
	    ivol,
	    {
	      complex t;
	      spincolor_scalar_prod(t,q_dag[ivol],q[ivol]);
	      complex_summassign(loc[ivol],t);
	    });
      }
    
    glb_reduce(&res,loc,locVol);
  }
  
  //compute all the meson contractions
  void compute_mes2pts_contr(int normalize)
  //void compute_mes2pts_contr(int normalize)
  {
    master_printf("Computing meson 2pts_contractions\n");
    
    // Tr [ GSO G5 S1^+ G5 GSI S2 ]      GSI is on the sink
    // (GSO)_{ij(i)} (G5)_{j(i)} (S1*)^{ab}_{kj(i)} (G5)_k (GSI)_{kl(k)} (S2)^{ab}_{l(k)i}
    //
    // A(i)=(GSO)_{ij(i)} (G5)_{j(i)}
    // B(k)=(G5)_k (GSI)_{kl(k)}
    //
    // A(i) (S1*)^{ab}_{kj(i)} B(k) (S2)^{ab}_{l(k)i}
    
    if(IS_MASTER_THREAD) mes2pts_contr_time-=take_time();
    
    for(size_t icombo=0;icombo<mes2pts_contr_map.size();icombo++)
      {
	master_printf("icombo %d/%d\n",icombo,mes2pts_contr_map.size());
	
	const qprop_t &Q1=Q[mes2pts_contr_map[icombo].a];
	const qprop_t &Q2=Q[mes2pts_contr_map[icombo].b];
	
	const double norm=12/sqrt(Q1.ori_source_norm2*Q2.ori_source_norm2); //12 in case of a point source
	for(size_t ihadr_contr=0;ihadr_contr<mes_gamma_list.size();ihadr_contr++)
	  {
	    const int& ig_so=mes_gamma_list[ihadr_contr].so;
	    const int& ig_si=mes_gamma_list[ihadr_contr].si;
	    
	    if(nso_spi==1 and ig_so!=5)
	      crash("implemented only g5 contraction on the source for non-diluted source");
	    
	    loc_contr->reset();
	    
	    for(int i=0;i<nso_spi;i++)
	      {
		const int& j=base_gamma[ig_so].pos[i];
		
		complex A;
		unsafe_complex_prod(A,base_gamma[ig_so].entr[i],base_gamma[5].entr[j]);
		
		complex AB[NDIRAC];
		for(int k=0;k<NDIRAC;k++)
		  {
		    //compute AB*norm
		    complex B;
		    unsafe_complex_prod(B,base_gamma[5].entr[k],base_gamma[ig_si].entr[k]);
		    unsafe_complex_prod(AB[k],A,B);
		    
		    if(normalize)
		      complex_prodassign_double(AB[k],norm);
		  }
		
		for(int b=0;b<nso_col;b++)
		  {
		    const LxField<spincolor>& q1=Q1[so_sp_col_ind(j,b)];
		    const LxField<spincolor>& q2=Q2[so_sp_col_ind(i,b)];
		    
		    NISSA_PARALLEL_LOOP(ivol,0,locVol)
		      {
			for(int k=0;k<NDIRAC;k++)
			  {
			    const int l=base_gamma[ig_si].pos[k];
			    
			    complex c={0,0};
			    for(int a=0;a<NCOL;a++)
			      complex_summ_the_conj1_prod(c,q1[ivol][k][a],q2[ivol][l][a]);
			    complex_summ_the_prod((*loc_contr)[ivol],c,AB[k]);
			  }
		      }
		    NISSA_PARALLEL_LOOP_END;
		  }
	      }
	    
	    THREAD_BARRIER();
	    complex temp_contr[glbSize[0]];
	    glb_reduce(temp_contr,*loc_contr,locVol,glbSize[0],locSize[0],glbCoordOfLoclx[0][0]);
	    
	    for(int t=0;t<glbSize[0];t++)
	      complex_summassign(mes2pts_contr[ind_mes2pts_contr(icombo,ihadr_contr,rel_time_of_glb_time(t))],temp_contr[t]);
	  }
      }
    
    //stats
    if(IS_MASTER_THREAD)
      {
	nmes2pts_contr_made+=mes2pts_contr_map.size()*mes_gamma_list.size();
	mes2pts_contr_time+=take_time();
      }
  }
  
  //print all mesonic 2pts contractions
  void print_mes2pts_contr(int n,int force_append,int skip_inner_header,const std::string &alternative_header_template)
  {
    //set the header template
    std::string header_template;
    if(alternative_header_template=="") header_template="\n # Contraction of %s ^ \\dag and %s\n\n";
    else header_template=alternative_header_template;
    
    contr_print_time-=take_time();
    
    //reduce and normalise
    // glb_nodes_reduce_complex_vect(mes2pts_contr,mes2pts_contr_size);
    // for(int i=0;i<mes2pts_contr_size;i++) complex_prodassign_double(mes2pts_contr[i],1.0);
    
    //list to open or append
    open_or_append_t list;
    
    for(size_t icombo=0;icombo<mes2pts_contr_map.size();icombo++)
      {
	auto& combo=mes2pts_contr_map[icombo];
	
	//path to use
	FILE *fout=list.open(combine("%s/%s_%s",outfolder,mes2pts_prefix.c_str(),combo.name.c_str()),force_append);
	
	master_fprintf(fout,header_template.c_str(),combo.a.c_str(),combo.b.c_str());
	
	print_contractions_to_file(fout,mes_gamma_list,mes2pts_contr+ind_mes2pts_contr(icombo,0,0),0,"",1.0/n,skip_inner_header);
	master_fprintf(fout,"\n");
	
	//close the file
	close_file(fout);
      }
    
    contr_print_time+=take_time();
  }
  
  //////////////////////////////////////// barionic contractions //////////////////////////////////////////////////////////
  
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
    const dirac_matr g2g4=base_gamma[2]*base_gamma[4];
    const dirac_matr C=dirac_prod_idouble(g2g4,1.0);
    Cg5=dirac_prod(C,base_gamma[5]);
  }
  
  //allocate bariionic contr
  void allocate_bar2pts_contr()
  {
    bar2pts_contr_size=ind_bar2pts_contr(bar2pts_contr_map.size(),0,0);
    bar2pts_alt_contr_size=ind_bar2pts_alt_contr(bar2pts_contr_map.size(),0,0,0);
    bar2pts_contr=nissa_malloc("bar2pts_contr",bar2pts_contr_size,complex);
    bar2pts_alt_contr=nissa_malloc("bar2pts_alt_contr",bar2pts_alt_contr_size,complex);
  }
  
  //free them
  void free_bar2pts_contr()
  {
    nissa_free(bar2pts_contr);
    nissa_free(bar2pts_alt_contr);
  }
  
  //Compute all contractions for 1/2 and 3/2 barions.  The calculation
  //is organized in three blocks, one corresponding to g5 g5 matrix
  //elements, the two others to the gi gi combination, finally to the
  //gi gj one.  The three possible Wick contractions are obtained by
  //simply permutating the sink Dirac index of the three propagators,
  //such that the "direct" contraction correspond to iWick=0, exchange
  //to one of the two others (should they be simply related by time
  //reversal?)
  //
  // O_ga=          eps_{a,b,c} q1_{al,a}(CG)_{al,be}q2_{be,b}q3_{ga,c}
  //
  // O^\dagger_ga = eps_{a,b,c}} q3^*_{ga,c}q2^*_{be,b}(CG)^\dag_{al,be}q3_{ga,c}
  void compute_bar2pts_alt_contr()
  {
    crash("reimplement");
    
    //     master_printf("Computing barion 2pts contractions alternative\n");
    
//     const int nWicks=2;
    
//     UNPAUSE_TIMING(bar2pts_alt_contr_time);
//     for(size_t icombo=0;icombo<bar2pts_contr_map.size();icombo++)
//       {
// 	qprop_t &Q1=Q[bar2pts_contr_map[icombo].a];
// 	qprop_t &Q2=Q[bar2pts_contr_map[icombo].b];
// 	qprop_t &Q3=Q[bar2pts_contr_map[icombo].c];
// 	const double norm=pow(12,1.5)/sqrt(Q1.ori_source_norm2*Q2.ori_source_norm2*Q3.ori_source_norm2); //12 is even in case of a point source
	
// 	for(auto &iProjGroup : std::array<std::array<int,3>,NBAR_ALT_SINGLE_PROJ>
// #ifdef BAR_ALT_LIMITED_PROJ
// 	       {{{5,5, 0},
// 		 {1,1, 1},{2,2, 1},{3,3, 1},
// 		 {1,2, 2},{1,3, 2},{2,1, 2},{2,3, 2},{3,1, 2},{3,2, 2}}})
// #else
// 	      {{{5,5, 0},
// 		{1,1, 1},{2,2, 1},{3,3, 1},
// 		{1,2, 2},{1,3, 2},{2,1, 2},{2,3, 2},{3,1, 2},{3,2, 2},
// 		{4,4, 3},
// 		{4,1, 4},{4,2, 4},{4,3, 4},
// 		{1,4, 5},{2,4, 5},{3,4, 5}}})
// #endif
// 	  {
// 	    const int igSo=iProjGroup[0];
// 	    const int igSi=iProjGroup[1];
	    
// 	    complex contr[glbSize[0]*nWicks];
// 	    tm_corr_op::compute_baryon_2pts_proj_contr(contr,igSo,igSi,Q1.sp,Q2.sp,Q3.sp,source_coord[0],temporal_bc);
	    
// 	    for(int dt=0;dt<glbSize[0];dt++)
// 	      for(int iWick=0;iWick<nWicks;iWick++)
// 		{
// 		  /// Input index
// 		  const int iin=iWick+nWicks*dt;
		  
// 		  /// Out index
// 		  const int iout=ind_bar2pts_alt_contr(icombo,iWick,iProjGroup[2],dt);
		  
// 		  complex_summ_the_prod_double(bar2pts_alt_contr[iout],contr[iin],norm);
// 		}
// 	  }
//   }
//   STOP_TIMING(bar2pts_alt_contr_time);
  
//   //stats
//   if(IS_MASTER_THREAD) nbar2pts_alt_contr_made+=bar2pts_contr_map.size();
  }
  
  // //compute all contractions
  // void compute_bar2pts_contr_free_theory()
  // {
    
  //   master_printf("Computing barion contractions\n");
    
  //   //local thread/node contractions
  //   complex *loc_contr=new complex[bar2pts_contr_size];
  //   memset(loc_contr,0,sizeof(complex)*bar2pts_contr_size);
    
  //   int eps[3][3][3]={};
  //   for(int i=0;i<3;i++)
  //     {
  // 	eps[i][(i+1)%3][(i+2)%3]=+1;
  // 	eps[i][(i+2)%3][(i+1)%3]=-1;
  //     }
    
  //   for(int i=0;i<3;i++)
  //     for(int j=0;j<3;j++)
  // 	for(int k=0;k<3;k++)
  // 	  master_printf("%d %d %d %d\n",i,j,k,eps[i][j][k]);
    
  //   std::vector<std::array<int,3>> list({{{5,5,0},
  // 	  {1,1,1},{2,2,1},{3,3,1},
  // 					 {1,2,2},{1,3,2},{2,1,2},{2,3,2},{3,1,2},{3,2,2}}});
  //   for(int iP=0;iP<10;iP++)
  //     {
  // 	int g_si=list[iP][0];
  // 	int g_so=list[iP][1];
        
  // 	spinspin opg0[2];
  // 	for(int it=0;it<2;it++)
  // 	  {
  // 	    dirac_matr temp=base_gamma[g_si]*base_gamma[g_so];
  // 	    dirac_matr temp0=temp*base_gamma[4];
	    
  // 	    spinspin_dirac_prod_double(opg0[it],&temp,1.0);
  // 	    if(it==0)
  // 	      spinspin_dirac_summ_the_prod_double(opg0[it],&temp0,1.0);
  // 	    else
  // 	      spinspin_dirac_summ_the_prod_double(opg0[it],&temp0,-1.0);
  // 	  }
	
  // 	spinspin Cg_so;
  // 	dirac_matr temp=base_gamma[4]*herm(set_CgX(g_so))*base_gamma[4];
  // 	spinspin_dirac_prod_double(Cg_so,&temp,1.0);
  // 	//
  // 	spinspin Cg_si;
  // 	temp=set_CgX(g_si);
  // 	spinspin_dirac_prod_double(Cg_si,&temp,1.0);
	
  // 	for(size_t icombo=0;icombo<bar2pts_contr_map.size();icombo++)
  // 	  {
  // 	    qprop_t &Q1=Q[bar2pts_contr_map[icombo].a];
  // 	    qprop_t &Q2=Q[bar2pts_contr_map[icombo].b];
  // 	    qprop_t &Q3=Q[bar2pts_contr_map[icombo].c];
  // 	    double norm=6*pow(12,1.5)/sqrt(Q1.ori_source_norm2*Q2.ori_source_norm2*Q3.ori_source_norm2); //12 is even in case of a point source
	    
  // 	    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
  // 	      {
  // 		const int t=rel_time_of_loclx(ivol);
  // 		const int it=0;//ANNA (t<=glb_size[0]/2)?0:1;
		
  // 		su3spinspin p1,p2,p3;
		
  // 		//Takes a slice
  // 		for(auto &k : std::vector<std::pair<su3spinspin&,qprop_t&>>{{p1,Q1},{p2,Q2},{p3,Q3}})
  // 		  for(int sp_so=0;sp_so<NDIRAC;sp_so++)
  // 		    for(int sp_si=0;sp_si<NDIRAC;sp_si++)
  // 		      for(int co_so=0;co_so<NCOL;co_so++)
  // 			for(int co_si=0;co_si<NCOL;co_si++)
  // 			  complex_copy(k.first[co_si][co_so][sp_si][sp_so],k.second[so_sp_col_ind(sp_so,co_so)][ivol][sp_si][co_si]);
		
  // 		const int a_so=0,a_si=0,b_so=0,b_si=0,c_so=0,c_si=0;
  // 		// for(int a_so=0;a_so<NCOL;a_so++)
  // 		//   for(int b_so=0;b_so<NCOL;b_so++)
  // 		//     for(int c_so=0;c_so<NCOL;c_so++)
  // 		//       for(int a_si=0;a_si<NCOL;a_si++)
  // 		// 	for(int b_si=0;b_si<NCOL;b_si++)
  // 		// 	  for(int c_si=0;c_si<NCOL;c_si++)
  // 			    for(int al_so=0;al_so<NDIRAC;al_so++)
  // 			      for(int be_so=0;be_so<NDIRAC;be_so++)
  // 				for(int ga_so=0;ga_so<NDIRAC;ga_so++)
  // 				  for(int al_si=0;al_si<NDIRAC;al_si++)
  // 				    for(int be_si=0;be_si<NDIRAC;be_si++)
  // 				      for(int ga_si=0;ga_si<NDIRAC;ga_si++)
  // 					{
  // 					  complex temp={1,0};
  // 					  complex_prodassign(temp,p1[a_si][a_so][al_si][al_so]);
  // 					  complex_prodassign(temp,p2[b_si][b_so][be_si][be_so]);
  // 					  complex_prodassign(temp,p3[c_si][c_so][ga_si][ga_so]);
  // 					  complex_prodassign(temp,Cg_si[al_si][be_si]);
  // 					  complex_prodassign(temp,Cg_so[al_so][be_so]);
  // 					  complex_prodassign(temp,opg0[it][ga_so][ga_si]);
					  
  // 					  complex_summ_the_prod_double(loc_contr[ind_bar2pts_contr(icombo,0,iP,t)],temp,-// eps[a_so][b_so][c_so]*eps[a_si][b_si][c_si]*
  // 								       norm);
					  
  // 					  complex temp2={1,0};
  // 					  complex_prodassign(temp2,p1[a_si][c_so][al_si][ga_so]);
  // 					  complex_prodassign(temp2,p2[b_si][b_so][be_si][be_so]);
  // 					  complex_prodassign(temp2,p3[c_si][a_so][ga_si][al_so]);
  // 					  complex_prodassign(temp2,Cg_si[al_si][be_si]);
  // 					  complex_prodassign(temp2,Cg_so[al_so][be_so]);
  // 					  complex_prodassign(temp2,opg0[it][ga_so][ga_si]);
					  
  // 					  complex_subt_the_prod_double(loc_contr[ind_bar2pts_contr(icombo,1,iP,t)],temp2,-// eps[a_so][b_so][c_so]*eps[a_si][b_si][c_si]*
  // 								       norm);
					  
  // 					}
  // 	      }
  //NISSA_PARALLEL_LOOP_END;
  // 	  }
  //     }
    
  //   //reduce
  //   complex *master_reduced_contr=(complex*)glb_threads_reduce_double_vect((double*)loc_contr,bar2pts_contr_size*2);
  //   NISSA_PARALLEL_LOOP(i,0,bar2pts_contr_size)
  //     {
  // 	//remove border phase
  // 	int t=i%glb_size[0];
  // 	double arg=3*temporal_bc*M_PI*t/glb_size[0];
  // 	complex phase={cos(arg),sin(arg)};
  // 	complex_summ_the_prod(bar2pts_contr[i],master_reduced_contr[i],phase);
  //     }
  //NISSA_PARALLEL_LOOP_END;
  //   THREAD_BARRIER();
  //   delete[] loc_contr;
  // }
  
  //compute all contractions
  void compute_bar2pts_contr()
  {
    master_printf("Computing barion 2pts contractions\n");
    
    crash("#warning reimplement");
    //allocate loc storage
    complex *loc_contr=new complex[bar2pts_contr_size];
    memset(loc_contr,0,sizeof(complex)*bar2pts_contr_size);
    
    //void (*list_fun[2])(complex,const complex,const complex)={complex_summ_the_prod,complex_subt_the_prod};
    UNPAUSE_TIMING(bar2pts_contr_time);
    for(size_t icombo=0;icombo<bar2pts_contr_map.size();icombo++)
      {
	// qprop_t &Q1=Q[bar2pts_contr_map[icombo].a];
	// qprop_t &Q2=Q[bar2pts_contr_map[icombo].b];
	// qprop_t &Q3=Q[bar2pts_contr_map[icombo].c];
	//double norm=pow(12,1.5)/sqrt(Q1.ori_source_norm2*Q2.ori_source_norm2*Q3.ori_source_norm2); //12 is even in case of a point source
	
	for(int al=0;al<NDIRAC;al++)
	  for(int ga=0;ga<NDIRAC;ga++)
	    for(int b=0;b<NCOL;b++)
	      for(int iperm=0;iperm<2;iperm++)
		{
		  NISSA_PARALLEL_LOOP(ivol,0,locVol)
		    {
		      // const int eps[3][2]={{1,2},{2,0},{0,1}};
		      // int sign[2]={1,-1};
    
		      // int t=rel_time_of_loclx(ivol);
		      
		      // int c=eps[b][iperm],a=eps[b][!iperm];
		      // int be=Cg5.pos[al];
		  
		      // int ga1_l[2][NDIRAC]={{0,1,2,3},{2,3,0,1}}; //ga1 index for 1 or gamma0 matrix
		      // int sign_idg0[2]={1,(t<(glb_size[0]/2))?-1:+1}; //gamma0 is -1 always
    			  for(int al1=0;al1<NDIRAC;al1++)
    			    for(int b1=0;b1<NCOL;b1++)
    			      {
    				// complex diquark_dir={0,0},diquark_exc={0,0};
				
    				//build the diquark
    				for(int iperm1=0;iperm1<2;iperm1++)
    				  {
    				    // int c1=eps[b1][iperm1],a1=eps[b1][!iperm1];
				    
    				    for(int idg0=0;idg0<2;idg0++)
    				      {
    					// int isign=((sign[iperm]*sign[iperm1]*sign_idg0[idg0])==1);
    					// int ga1=ga1_l[idg0][ga];
					
    					// list_fun[isign](diquark_dir,Q1[so_sp_col_ind(al,a)][ivol][al1][a1],Q3[so_sp_col_ind(ga,c)][ivol][ga1][c1]); //direct
    					// list_fun[isign](diquark_exc,Q1[so_sp_col_ind(ga,c)][ivol][al1][a1],Q3[so_sp_col_ind(al,a)][ivol][ga1][c1]); //exchange
    				      }
    				  }
				
    				//close it
    				// complex w;
    				// unsafe_complex_prod(w,Cg5.entr[al1],Cg5.entr[al]);
				// int be1=Cg5.pos[al1];
    				// complex_prodassign_double(diquark_dir,w[RE]*norm);
    				// complex_prodassign_double(diquark_exc,w[RE]*norm);
				// complex_summ_the_prod(loc_contr[ind_bar2pts_contr(icombo,0,t)],Q2[so_sp_col_ind(be,b)][ivol][be1][b1],diquark_dir);
 				// complex_summ_the_prod(loc_contr[ind_bar2pts_contr(icombo,1,t)],Q2[so_sp_col_ind(be,b)][ivol][be1][b1],diquark_exc);
    			      }
		    }
		  NISSA_PARALLEL_LOOP_END;
		}
      }
    STOP_TIMING(bar2pts_contr_time);
    
    //reduce
    crash("#warning");
    // complex *master_reduced_contr=glb_threads_reduce_complex_vect(loc_contr,bar2pts_contr_size);
    // NISSA_PARALLEL_LOOP(i,0,bar2pts_contr_size)
    //   {
    // 	//remove border phase
    // 	int t=i%glb_size[0];
    // 	double arg=3*temporal_bc*M_PI*t/glb_size[0];
    // 	complex phase={cos(arg),sin(arg)};
    // 	complex_summ_the_prod(bar2pts_contr[i],master_reduced_contr[i],phase);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    delete[] loc_contr;
    
    //stats
    if(IS_MASTER_THREAD) nbar2pts_contr_made+=bar2pts_contr_map.size();
  }
  
  //print all contractions
  void print_bar2pts_contr()
  {
    //list to open or append
    open_or_append_t list;
    
    //reduce
    crash("#warning reimplement");
    //glb_nodes_reduce_complex_vect(bar2pts_contr,bar2pts_contr_size);
    
    for(size_t icombo=0;icombo<bar2pts_contr_map.size();icombo++)
	for(int dir_exc=0;dir_exc<2;dir_exc++)
	{
	  //open output
	  FILE *fout=list.open(combine("%s/bar_contr_%s_%s",outfolder,(dir_exc==0)?"dir":"exc",bar2pts_contr_map[icombo].name.c_str()));
	  for(int t=0;t<glbSize[0];t++)
	    {
	      //normalize for nsources and 1+g0
	      complex c;
	      complex_prod_double(c,bar2pts_contr[ind_bar2pts_contr(icombo,dir_exc,t)],1.0/(2*nhits));
	      master_fprintf(fout,"%+16.16lg %+16.16lg\n",c[RE],c[IM]);
	    }
	  
	  close_file(fout);
	}
  }
  
  //print all contractions
  void print_bar2pts_alt_contr()
  {
    //list to open or append
    open_or_append_t list;
    
    for(size_t icombo=0;icombo<bar2pts_contr_map.size();icombo++)
      for(int iProj=0;iProj<NBAR_ALT_PROJ;iProj++)
	for(int iWick=0;iWick<2;iWick++)
	  {
	    //open output
	    FILE *fout=list.open(combine("%s/bar_alt_contr_%s_proj_%d_Wick_%d",outfolder,bar2pts_contr_map[icombo].name.c_str(),iProj,iWick));
	    
	    master_fprintf(fout,"# proj %d\n",iProj);
	    for(int t=0;t<glbSize[0];t++)
	      {
		//normalize for nsources and 1+g0
		complex c;
		complex_prod_double(c,bar2pts_alt_contr[ind_bar2pts_alt_contr(icombo,iWick,iProj,t)],1.0/(2*nhits));
		master_fprintf(fout,"%+16.16lg %+16.16lg\n",c[RE],c[IM]);
	      }
	    
	    close_file(fout);
	  }
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //compute the matrix element of the conserved current between two propagators. If asking to revert, g5 is inserted between the two propagators
  void conserved_vector_current_mel(const LxField<quad_su3>& conf,
				    LxField<spin1field>& si,
				    const dirac_matr& ext_g,
				    const int& r,
				    const char* name_bw,
				    const char* name_fw,
				    const bool& revert)
  {
    si.reset();
    
    //compute the gammas
    dirac_matr GAMMA[5];
    dirac_matr temp_gamma;
    if(twisted_run>0)	temp_gamma=dirac_prod_idouble(base_gamma[5],-tau3[r]);
    else                temp_gamma=base_gamma[0];
    
    //Add g5 on the gamma, only if asked to revert
    if(revert)
      {
	GAMMA[4]=base_gamma[5]*temp_gamma;
	for(int mu=0;mu<NDIM;mu++)
	 GAMMA[mu]=base_gamma[5]*base_gamma[igamma_of_mu[mu]];
      }
    else
      {
	GAMMA[4]=temp_gamma;
	for(int mu=0;mu<NDIM;mu++) GAMMA[mu]=base_gamma[igamma_of_mu[mu]];
      }
    
    dirac_matr g;
    if(revert) g=dirac_prod(ext_g,base_gamma[5]);
    else       g=ext_g;
    
    for(int iso_spi_bw=0;iso_spi_bw<nso_spi;iso_spi_bw++)
      for(int iso_col=0;iso_col<nso_col;iso_col++)
	{
	  int iso_spi_fw=g.pos[iso_spi_bw];
	  int idc_fw=so_sp_col_ind(iso_spi_fw,iso_col);
	  int idc_bw=so_sp_col_ind(iso_spi_bw,iso_col);
	  
	  //get componentes
	  LxField<spincolor>& Qfw=Q[name_fw][idc_fw];
	  LxField<spincolor>& Qbw=Q[name_bw][idc_bw];
	  
	  Qfw.updateHalo();
	  Qbw.updateHalo();
	  conf->updateHalo();
	  
	  PAR(0,locVol,
	      CAPTURE(g,iso_spi_bw,GAMMA,
		      TO_WRITE(si),
		      TO_READ(Qfw),
		      TO_READ(Qbw),
		      TO_READ(conf)),
	      ivol,
	      {
		for(int mu=0;mu<NDIM;mu++)
		  {
		    const int ivol_fw=loclxNeighup[ivol][mu];
		    spincolor f,Gf;
		    complex c;
		    
		    //piece psi_ivol U_ivol psi_fw
		    unsafe_su3_prod_spincolor(f,conf[ivol][mu],Qfw[ivol_fw]);
		    unsafe_dirac_prod_spincolor(Gf,GAMMA[4],f);
		    dirac_subt_the_prod_spincolor(Gf,GAMMA[mu],f);
		    spincolor_scalar_prod(c,Qbw[ivol],Gf);
		    complex_prodassign(c,g.entr[iso_spi_bw]);
		    complex_summ_the_prod_idouble(si[ivol][mu],c,-0.5);
		    
		    //piece psi_fw U_ivol^dag psi_ivol
		    unsafe_su3_dag_prod_spincolor(f,conf[ivol][mu],Qfw[ivol]);
		    unsafe_dirac_prod_spincolor(Gf,GAMMA[4],f);
		    dirac_summ_the_prod_spincolor(Gf,GAMMA[mu],f);
		    spincolor_scalar_prod(c,Qbw[ivol_fw],Gf);
		    complex_prodassign(c,g.entr[iso_spi_bw]);
		    complex_summ_the_prod_idouble(si[ivol][mu],c,+0.5);
		  }
	      });
	}
  }
  
  /// Compute the matrix element of the current between two propagators
  void vector_current_mel(LxField<spin1field>& si,
			  const dirac_matr& ext_g,
			  const int& r,
			  const char* id_Qbw,
			  const char* id_Qfw,
			  const bool& revert)
  {
    si.reset();
    
    dirac_matr GAMMA[NDIM];
    for(int mu=0;mu<NDIM;mu++)
      if(revert) GAMMA[mu]=base_gamma[5]*base_gamma[igamma_of_mu[mu]];
      else       GAMMA[mu]=base_gamma[igamma_of_mu[mu]];
    
    dirac_matr g;
    if(revert) g=dirac_prod(ext_g,base_gamma[5]);
    else       g=ext_g;
    
    for(int iso_spi_fw=0;iso_spi_fw<nso_spi;iso_spi_fw++)
      for(int iso_col=0;iso_col<nso_col;iso_col++)
	{
	  int iso_spi_bw=g.pos[iso_spi_fw];
	  
	  //get componentes
	  LxField<spincolor>& Qbw=Q[id_Qbw][so_sp_col_ind(iso_spi_bw,iso_col)];
	  LxField<spincolor>& Qfw=Q[id_Qfw][so_sp_col_ind(iso_spi_fw,iso_col)];
	  
	  PAR(0,locVol,
	      CAPTURE(g,GAMMA,iso_spi_fw,
		      TO_WRITE(si),
		      TO_READ(Qfw),
		      TO_READ(Qbw)),
	      ivol,
	      {
		for(int mu=0;mu<NDIM;mu++)
		  {
		    spincolor temp;
		    complex c;
		    
		    unsafe_dirac_prod_spincolor(temp,GAMMA[mu],Qfw[ivol]);
		    spincolor_scalar_prod(c,Qbw[ivol],temp);
		    complex_summ_the_prod(si[ivol][mu],c,g.entr[iso_spi_fw]);
		  }
	      });
	}
  }
  
  //compute local or conserved vector current matrix element
  void local_or_conserved_vector_current_mel(LxField<spin1field>& si,
					     const dirac_matr &g,
					     const std::string &prop_name_bw,
					     const std::string &prop_name_fw,
					     const bool& revert)
  {
    if(loc_hadr_curr)
      vector_current_mel(si,g,Q[prop_name_fw].r,prop_name_bw.c_str(),prop_name_fw.c_str(),revert);
    else
      {
	momentum_t plain_bc;
	plain_bc[0]=temporal_bc;
	for(int mu=1;mu<NDIM;mu++)
	  plain_bc[mu]=0.0;
	
	const LxField<quad_su3>* conf=get_updated_conf(Q[prop_name_fw].charge,plain_bc,*glb_conf);
	const int r=Q[prop_name_fw].r;
	conserved_vector_current_mel(*conf,si,g,r,prop_name_bw.c_str(),prop_name_fw.c_str(),revert);
      }
  }
  
  //                                                          handcuffs
  
  void compute_handcuffs_contr()
  {
    master_printf("Computing handcuffs contractions\n");
    
    //allocate all sides
    std::map<std::string,LxField<spin1field>*> sides;
    
    //store source norm for all sides
    std::map<std::string,double> normSides;
    
    //loop over sides
    for(size_t iside=0;iside<handcuffs_side_map.size();iside++)
      {
	//allocate
	handcuffs_side_map_t &h=handcuffs_side_map[iside];
	std::string side_name=h.name;
	
        normSides.insert({side_name,Q[h.fw].ori_source_norm2*Q[h.bw].ori_source_norm2});
	
	sides[side_name]=new LxField<spin1field>(side_name.c_str());
	LxField<spin1field>& si=*sides[side_name];
	si.reset();
	
	//check r are the same (that is, opposite!)
	if(twisted_run and (not loc_hadr_curr)) {
	  if(Q[h.fw].r==Q[h.bw].r and (not Q[h.bw].is_source))
	    crash("conserved current needs opposite r (before reverting), but quarks %s and %s have the same",h.fw.c_str(),h.bw.c_str());
	}
	
	//compute dirac combo
	dirac_matr g;
	const int ig=::abs(handcuffs_side_map[iside].igamma);
	const int revert=(handcuffs_side_map[iside].igamma>=0); //reverting only if positive ig asked
	if(ig!=5 and not diluted_spi_source)
	  crash("ig %d not available if not diluting in spin",ig);
	//dirac_prod(&g,base_gamma+5,base_gamma+ig);
	g=base_gamma[ig];
	
	//compute the matrix element
	local_or_conserved_vector_current_mel(si,g,h.bw,h.fw,revert);
	
	//if(h.store) store_spin1field(combine("%s/handcuff_side_%s",outfolder,h.name.c_str()),si);
      }
    
    //add the photon
    for(size_t ihand=0;ihand<handcuffs_map.size();ihand++)
      {
	handcuffs_map_t &h=handcuffs_map[ihand];
	std::string right_with_photon=h.right+"_photon";
	master_printf("Inserting photon to compute %s\n",right_with_photon.c_str());
	
	if(sides.find(right_with_photon)==sides.end())
	  {
	    auto rp=sides[right_with_photon]=new LxField<spin1field>(right_with_photon.c_str());
	    multiply_by_tlSym_gauge_propagator(*rp,*sides[h.right],photon);
	  }
      }
    
    //compute the hands
    for(int ihand=0;ihand<handcuffs_map.size();ihand++)
      if(sides.find(handcuffs_map[ihand].left)==sides.end() or
	 sides.find(handcuffs_map[ihand].right)==sides.end())
	crash("Unable to find sides: %s or %s",handcuffs_map[ihand].left.c_str(),handcuffs_map[ihand].right.c_str());
      else
	{
	  const auto& left=*sides[handcuffs_map[ihand].left];
	  const auto& right=*sides[handcuffs_map[ihand].right+"_photon"];
	  
	  loc_contr->reset();
	  NISSA_PARALLEL_LOOP(ivol,0,locVol)
	    for(int mu=0;mu<NDIM;mu++)
	      complex_subt_the_prod((*loc_contr)[ivol],
				    left[ivol][mu],
				    right[ivol][mu]);
	  NISSA_PARALLEL_LOOP_END;
	  
          double normalization=glbSpatVol*144.0;
          normalization/=sqrt((double)normSides[handcuffs_map[ihand].left]*normSides[handcuffs_map[ihand].right]);
	  
	  complex temp[1];
	  glb_reduce(temp,*loc_contr,locVol);
	  complex_summ_the_prod_double(handcuffs_contr[ind_handcuffs_contr(ihand)],temp[0],normalization);
	}
    
    //free
    for(std::map<std::string,LxField<spin1field>*>::iterator it=sides.begin();it!=sides.end();it++)
      delete it->second;
  }
  
  //allocate handcuff contractions
  void allocate_handcuffs_contr()
  {
    handcuffs_contr_size=ind_handcuffs_contr(handcuffs_map.size()-1)+1;
    handcuffs_contr=nissa_malloc("handcuffs_contr",handcuffs_contr_size,complex);
  }
  
  //free handcuff contractions
  void free_handcuffs_contr()
  {
    nissa_free(handcuffs_contr);
  }
  
  //print handcuffs contractions
  void print_handcuffs_contr()
  {
    //list to open or append
    open_or_append_t list;
    
    contr_print_time-=take_time();
    
    //Open if size different from zero
    FILE *fout=nullptr;
    if(handcuffs_map.size()) fout=list.open(combine("%s/handcuffs",outfolder));
    
    for(size_t icombo=0;icombo<handcuffs_map.size();icombo++)
      {
	//normalize for nsources
	complex c;
	complex_prod_double(c,handcuffs_contr[ind_handcuffs_contr(icombo)],1.0/nhits);
	master_fprintf(fout,"%s %+16.16lg %+16.16lg\n",handcuffs_map[icombo].name.c_str(),c[RE],c[IM]);
      }
    //close the file
    if(fout) close_file(fout);
    
    contr_print_time+=take_time();
  }
}
