#include <nissa.hpp>

#include "conf.hpp"
#include "pars.hpp"
#include "prop.hpp"

// #include <immintrin.h>

using namespace nissa;

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
     
     we store them separately
   
 */

int nsmear,ntranspose,ncontract;
double smear_time=0,transpose_time=0,contract_time=0;

//hadron contractions
struct bar_triplet_t : std::vector<double>
{
  int a,b,c;
  bar_triplet_t(int a,int b,int c) : a(a),b(b),c(c) {}
};
std::vector<bar_triplet_t> prop_hadr_combo_map;

bool compare_el(std::pair<bar_triplet_t,int> a,std::pair<bar_triplet_t,int> b)
{
  return ((a.first.a<b.first.a)||
	  ((a.first.a==b.first.a)&&
	   ((a.first.b<b.first.b)||
	    ((a.first.b==b.first.b)&&
	     (a.first.c<b.first.c)))));
}

//index inside a colorspinspin
inline int dirspin_ind(int ic_si,int ic_so,int id_si,int id_so,int ri)
{return
    (ri+2*
     (id_so+4*
      (id_si+4*
       (ic_so+NCOL*ic_si))));
}
inline int ind_rot_prop(int iprop,int ci,int ivol)
{return
    (ivol+loc_vol*
     (ci+sizeof(su3spinspin)/sizeof(double)*iprop));
}
inline int ind_rot_prop(int iprop,int ic_si,int ic_so,int id_si,int id_so,int ri,int ivol)
{return ind_rot_prop(iprop,dirspin_ind(ic_si,ic_so,id_si,id_so,ri),ivol);}

//contains how to bind the indices of the three propagators to make a Wick contraction
const int neps_perm=6;
int eps_perm[neps_perm][4]={{0,1,2,+1},{1,2,0,+1},{2,0,1,+1},
			    {0,2,1,-1},{2,1,0,-1},{1,0,2,-1}};
dirac_matr Cg5;
int nwick;
int ind_wick(int dir_exc,int ri,int id_g0)
{return id_g0+2*(ri+2*dir_exc);}
std::vector<std::vector<std::pair<bar_triplet_t,double> > > wickes;
void set_wickes_contractions()
{
  //set total wick contractions and resize
  nwick=ind_wick(2-1,2-1,2-1)+1;
  wickes.resize(nwick);
  
  //set C*g5
  dirac_matr g2g4,C;
  dirac_prod(&g2g4,base_gamma+2,base_gamma+4);
  dirac_prod_idouble(&C,&g2g4,1);
  dirac_prod(&Cg5,&C,base_gamma+5);
  
  int id_g0[2]={0,4};
  for(int iid_g0=0;iid_g0<2;iid_g0++)
    for(int iperm_si=0;iperm_si<neps_perm;iperm_si++)
      for(int iperm_so=0;iperm_so<neps_perm;iperm_so++)
	for(int al1=0;al1<4;al1++)
	  for(int al=0;al<4;al++)
	    for(int ga=0;ga<4;ga++)
	      {
		//find other indices
		int a =eps_perm[iperm_so][0],b =eps_perm[iperm_so][1],c =eps_perm[iperm_so][2];
		int a1=eps_perm[iperm_si][0],b1=eps_perm[iperm_si][1],c1=eps_perm[iperm_si][2];
		int be1=Cg5.pos[al1],be=Cg5.pos[al];
		int ga1=base_gamma[id_g0[iid_g0]].pos[ga];
		//build dirac structure
		complex Cg5Cg5,PCg5Cg5,PCg5Cg5_eps;
		unsafe_complex_prod(Cg5Cg5,Cg5.entr[al1],Cg5.entr[al]);
		unsafe_complex_prod(PCg5Cg5,base_gamma[id_g0[iid_g0]].entr[al],Cg5Cg5);
		//include eps
		complex_prod_double(PCg5Cg5_eps,PCg5Cg5,eps_perm[iperm_si][3]*eps_perm[iperm_so][3]);
		
		//loop over real or imaginary part
		for(int ria=0;ria<2;ria++)
		  for(int rib=0;rib<2;rib++)
		    for(int ric=0;ric<2;ric++)
		      {
			complex wa={0,0},wb={0,0},wc={0,0};
			wa[ria]=wb[rib]=wc[ric]=1;
			complex wab,wabc;
			unsafe_complex_prod(wab,wa,wb);
			unsafe_complex_prod(wabc,wab,wc);
			complex w;
			unsafe_complex_prod(w,wabc,PCg5Cg5_eps);
			
			//loop over direct or exchange
			for(int idir_exc=0;idir_exc<2;idir_exc++)
			  {
			    const int itb=dirspin_ind(b1,b,be1,be,rib);
			    const int ita=dirspin_ind(a1,(idir_exc==0)?a:c,al1,(idir_exc==0)?al:ga,ria);
			    const int itc=dirspin_ind(c1,(idir_exc==0)?c:a,ga1,(idir_exc==0)?ga:al,ric);
			    
			    //find whether the real or imaginary part contributes
			    for(int ri=0;ri<2;ri++)
			      if(w[ri])
				wickes[ind_wick(idir_exc,ri,iid_g0)].push_back(std::make_pair(bar_triplet_t(ita,itb,itc),w[ri]));
			  }
		      }
	      }
  //sort the term (possibly to achieve faster access)
  for(int iwick=0;iwick<nwick;iwick++)
    std::sort(wickes[iwick].begin(),wickes[iwick].end(),compare_el);
}

//store correlation function
int nsm_sink=2;
double *corr;
int ind_corr(int icombo,int ism_sink,int ima,int ra,int imb,int rb,int imc,int rc,int iwick,int t)
{return
    (t+glb_size[0]*
     (iwick+nwick*
      (rc+nr*
       (imc+nqmass*
	(rb+nr*
	 (imb+nqmass*
	  (ra+nr*
	   (ima+nqmass*
	    (ism_sink+nsm_sink*icombo)))))))));
}
int corr_size;

//set all the combinations
void set_combinations()
{
  //add the contraction combination
  prop_hadr_combo_map.clear();
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_0));
  //scalar insertion on one of the three lines
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_S,PROP_0,PROP_0));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_S,PROP_0));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_S));
  //pseudoscalar insertion on one of the three lines
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_P,PROP_0,PROP_0));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_P,PROP_0));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_P));
  //tadpole insertion on one of the three lines
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_T,PROP_0,PROP_0));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_T,PROP_0));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_T));
  //self-energy insertion on one of the three lines
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_PHOTON2,PROP_0,PROP_0));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_PHOTON2,PROP_0));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_PHOTON2));
  //photon exchange between one of the three lines
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_PHOTON,PROP_PHOTON));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_PHOTON,PROP_0,PROP_PHOTON));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_PHOTON,PROP_PHOTON,PROP_0));
}
//init everything
void init_simulation(char *path)
{
  //set how to compute propagators, how to make barions and how to
  //combine the different kind of propagators
  set_inversions();
  set_wickes_contractions();
  set_combinations();
  
  //open input file and read it
  open_input(path);
  read_input_preamble();
  read_ape_smearing_pars();
  read_gaussian_smearing_pars();
  read_photon_pars();
  read_seed_start_random();
  read_free_theory_flag();
  read_random_gauge_transform();
  read_nsources();
  read_ngauge_conf();
  
  //allocate
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  ape_smeared_conf=nissa_malloc("ape_smeared_conf",loc_vol+bord_vol,quad_su3);
  photon_field=nissa_malloc("photon_phield",loc_vol+bord_vol,spin1field);
  source=nissa_malloc("source",loc_vol,su3spinspin);
  original_source=nissa_malloc("original_source",loc_vol,su3spinspin);
  corr_size=ind_corr(prop_hadr_combo_map.size()-1,nsm_sink-1,nqmass-1,nr-1,nqmass-1,nr-1,nqmass-1,nr-1,nwick-1,glb_size[0]-1)+1;
  corr=nissa_malloc("corr",corr_size,double);
  nqprop=iqprop(nqmass-1,nqprop_kind()-1,nr-1)+1;
  Q=nissa_malloc("Q*",nqprop,PROP_TYPE*);
  for(int iprop=0;iprop<nqprop;iprop++) Q[iprop]=nissa_malloc("Q",loc_vol+bord_vol,PROP_TYPE);
}

//init a new conf
void start_new_conf()
{
  setup_conf(conf,old_theta,put_theta,conf_path,rnd_gauge_transform,free_theory);
  //put periodic
  put_theta[0]=0;
  adapt_theta(conf,old_theta,put_theta,0,0);
  //spatial smearing
  ape_spatial_smear_conf(ape_smeared_conf,conf,ape_smearing_alpha,ape_smearing_niters);
  master_printf("Smeared plaquette: %.16lg\n",global_plaquette_lx_conf(ape_smeared_conf));
  //put back anti-periodic
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,0,0);
  
  //reset correlations
  vector_reset(corr);
}

//handle to discard the source
void skip_conf()
{
  for(int isource=0;isource<nsources;isource++)
    {
      coords coord;
      generate_random_coord(coord);
    }
}

//compute all correlations
THREADABLE_FUNCTION_0ARG(compute_correlations)
{
  GET_THREAD_ID();
  
  //local thread/node correlations
  double *loc_corr=new double[corr_size];
  memset(loc_corr,0,sizeof(double)*corr_size);
  
  //transposed propagators
  double *p=nissa_malloc("rotated_prop",nqprop*sizeof(su3spinspin)*loc_vol,double);
  
  for(int ism_sink=0;ism_sink<nsm_sink;ism_sink++)
    {
      //smear all sinks
      START_TIMING(smear_time,nsmear);
      if(ism_sink==1)
	for(int ip=0;ip<nqprop;ip++)
	gaussian_smearing(Q[ip],Q[ip],ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
      STOP_TIMING(smear_time);
      
      //transpose all propagators
      START_TIMING(transpose_time,ntranspose);
      for(int ip=0;ip<nqprop;ip++)
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  for(int ic_si=0;ic_si<NCOL;ic_si++)
	    for(int ic_so=0;ic_so<NCOL;ic_so++)
	      for(int id_si=0;id_si<4;id_si++)
		for(int id_so=0;id_so<4;id_so++)
		  for(int ri=0;ri<2;ri++)
		    p[ind_rot_prop(ip,ic_si,ic_so,id_si,id_so,ri,ivol)]=Q[ip][ivol][ic_si][ic_so][id_si][id_so][ri];
      THREAD_BARRIER();
      STOP_TIMING(transpose_time);
      
      UNPAUSE_TIMING(contract_time);
      for(size_t icombo=0;icombo<prop_hadr_combo_map.size();icombo++)
	for(int ima=0;ima<nqmass;ima++)
	  for(int imb=0;imb<nqmass;imb++)
	    for(int imc=0;imc<nqmass;imc++)
	      for(int ra=0;ra<nr;ra++)
		for(int rb=0;rb<nr;rb++)
		  for(int rc=0;rc<nr;rc++)
		    {
		      int ipa=iqprop(ima,prop_hadr_combo_map[icombo].a,ra);
		      int ipb=iqprop(imb,prop_hadr_combo_map[icombo].b,rb);
		      int ipc=iqprop(imc,prop_hadr_combo_map[icombo].c,rc);
		      
		      for(int iwick=0;iwick<nwick;iwick++)
			for(size_t iterm=0;iterm<wickes[iwick].size();iterm++)
			  {
			    int iwa=wickes[iwick][iterm].first.a;
			    int iwb=wickes[iwick][iterm].first.b;
			    int iwc=wickes[iwick][iterm].first.c;
			    
			    if(IS_MASTER_THREAD) ncontract++;
#ifndef BGQ
			    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
			      loc_corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,iwick,glb_coord_of_loclx[ivol][0])]+=
			      wickes[iwick][iterm].second*
			      p[ind_rot_prop(ipa,iwa,ivol)]*
			      p[ind_rot_prop(ipb,iwb,ivol)]*
			      p[ind_rot_prop(ipc,iwc,ivol)];
#else
			    for(int t=0;t<loc_size[0];t++)
			      {
				DECLARE_REG_BI_COMPLEX(reg_tot);
				REG_SPLAT_BI_COMPLEX(reg_tot,0);
				
				NISSA_CHUNK_WORKLOAD(START,CHUNK_LOAD,END,0,loc_spat_vol/4,thread_id,NACTIVE_THREADS);
				
				complex *a_ptr=(complex*)(p+ind_rot_prop(ipa,iwa,t*loc_spat_vol+START*4));
				complex *b_ptr=(complex*)(p+ind_rot_prop(ipb,iwb,t*loc_spat_vol+START*4));
				complex *c_ptr=(complex*)(p+ind_rot_prop(ipc,iwc,t*loc_spat_vol+START*4));
				
				DECLARE_REG_BI_COMPLEX(reg_a);
				DECLARE_REG_BI_COMPLEX(reg_b);
				DECLARE_REG_BI_COMPLEX(reg_c);
				for(int i=START;i<END;i++)
				  {
				    REG_LOAD_BI_COMPLEX(reg_a,a_ptr);
				    REG_LOAD_BI_COMPLEX(reg_b,b_ptr);
				    REG_LOAD_BI_COMPLEX(reg_c,c_ptr);
				    
				    a_ptr+=2;
				    b_ptr+=2;
				    c_ptr+=2;
				    
				    BI_COMPLEX_PREFETCH(a_ptr);
				    BI_COMPLEX_PREFETCH(b_ptr);
				    BI_COMPLEX_PREFETCH(c_ptr);
				    
				    DECLARE_REG_BI_COMPLEX(reg_temp);
				    REG_BI_COMPLEX_PROD_4DOUBLE(reg_temp,reg_a,reg_b);
				    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(reg_tot,reg_tot,reg_temp,reg_c);
				  }
				double tot[4];
				STORE_REG_BI_COMPLEX(tot,reg_tot);
				for(int i=0;i<4;i++) loc_corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,iwick,t+loc_size[0]*rank_coord[0])]+=tot[i]*wickes[iwick][iterm].second;
			      }
#endif			    
			    // for(int t=0;t<loc_size[0];t++)
			    //   {
			    // 	__m256d tot=_mm256_setzero_pd();
			    // 	NISSA_PARALLEL_LOOP(ispat_vol_quarter,0,loc_spat_vol/4)
			    // 	  {
			    // 	    int ivol_base=t*loc_spat_vol+ispat_vol_quarter*4;
			    // 	    __m256d wa=_mm256_load_pd(p+ind_rot_prop(ipa,iwa,ivol_base));
			    // 	    __m256d wb=_mm256_load_pd(p+ind_rot_prop(ipb,iwb,ivol_base));
			    // 	    __m256d wc=_mm256_load_pd(p+ind_rot_prop(ipc,iwc,ivol_base));
			    // 	    tot=_mm256_fmadd_pd(_mm256_mul_pd(wa,wb),wc,tot);
			    // 	  }
			    // 	double tott[4] __attribute__((aligned(32)));
			    // 	_mm256_store_pd(tott,tot);
			    // 	for(int i=0;i<4;i++) loc_corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,iwick,t+loc_size[0]*rank_coord[0])]+=tott[i]*wickes[iwick][iterm].second;
			    //   }
			  }
		    }
      STOP_TIMING(contract_time);
    }
  nissa_free(p);
  
  //reduce
  double *master_reduced_corr=glb_threads_reduce_double_vect(loc_corr,corr_size);
  NISSA_PARALLEL_LOOP(i,0,corr_size) corr[i]+=master_reduced_corr[i];
  THREAD_BARRIER();
  delete[] loc_corr;
}
THREADABLE_FUNCTION_END

//print all correlations averaging
void print_correlations()
{
  //reduce
  glb_nodes_reduce_double_vect(corr,corr_size);
  
  //open output
  FILE *fout=open_file(combine("%s/corr",outfolder),"w");
  
  for(size_t icombo=0;icombo<prop_hadr_combo_map.size();icombo++)
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
				     qprop_list[prop_hadr_combo_map[icombo].a].shortname,
				     qprop_list[prop_hadr_combo_map[icombo].b].shortname,
				     qprop_list[prop_hadr_combo_map[icombo].c].shortname,
				     ism_sink,qmass[ima],qmass[imb],qmass[imc],ra,rb,rc,dir_exc);
		      for(int t=0;t<glb_size[0];t++)
			{
			  //reassembling real/imaginary
			  //reassembling (1+-g0)/2 depending on t<>T/2
			  //putting factor 1/nsources
			  int sign=(t<glb_size[0]/2)?1:-1;
			  complex c;
			  
			  //remove anti-periodic condition phase
			  double arg=3*M_PI*t/glb_size[0];
			  complex phase={cos(arg),sin(arg)};
			  
			  for(int ri=0;ri<2;ri++)
			    {
			      //take id and g0
			      int iwick_id=ind_wick(dir_exc,ri,0/*id*/);
			      int iwick_g0=ind_wick(dir_exc,ri,1/*g0*/);
			      
			      //the sign should be on g0, but this way the correlator is periodic
			      c[ri]=
				(corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,iwick_id,t)]*sign+
				 corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,iwick_g0,t)])/
				(2*nsources);
			    }
			  
			  //put the phase and print
			  safe_complex_prod(c,c,phase);
			  master_fprintf(fout,"%+016.016lg %+016.016lg\n",c[RE],c[IM]);
			}
		      master_fprintf(fout,"\n");
		    }
  
  close_file(fout);
}

//close deallocating everything
void close()
{
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g, of which:\n",tot_prog_time);
  master_printf(" - %02.2f%s to prepare %d photon stochastic propagators (%2.2gs avg)\n",photon_prop_time/tot_prog_time*100,"%",nphoton_prop_tot,photon_prop_time/nphoton_prop_tot);
  master_printf(" - %02.2f%s to prepare %d generalized sources (%2.2gs avg)\n",source_time/tot_prog_time*100,"%",nsource_tot,source_time/nsource_tot);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf("    of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d contractions (%2.2gs avg), %02.2f MFlops/rank\n",contract_time/tot_prog_time*100,"%",ncontract,contract_time/ncontract,
		(long long int)ncontract*4*loc_vol/(contract_time*1e6));
  master_printf(" - %02.2f%s to perform %d transpositions (%2.2gs avg)\n",transpose_time/tot_prog_time*100,"%",ntranspose,transpose_time/ntranspose);
  master_printf(" - %02.2f%s to perform %d smearing (%2.2gs avg)\n",smear_time/tot_prog_time*100,"%",nsmear,smear_time/nsmear);
  
  nissa_free(photon_field);
  nissa_free(original_source);
  nissa_free(source);
  nissa_free(conf);
  nissa_free(ape_smeared_conf);
  for(int iprop=0;iprop<nqprop;iprop++) nissa_free(Q[iprop]);
  nissa_free(Q);
  nissa_free(corr);
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
  while(iconf<ngauge_conf && enough_time && !file_exists("stop") && read_conf_parameters(iconf,skip_conf,finish_file_present))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      for(int isource=0;isource<nsources;isource++)
	{
	  master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
	  //shift the conf and create the stochastic photon field
	  random_shift_gauge_conf(conf,old_theta,put_theta);
	  generate_photon_stochastic_propagator();
	  //generate source and smear it
	  generate_original_source();
	  smear_time-=take_time();
	  gaussian_smearing(original_source,original_source,ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
	  smear_time+=take_time();
	  //compute prop and correlators
	  generate_quark_propagators();
	  compute_correlations();
	}
      
      //print out correlations
      print_correlations();
      
      //pass to the next conf if there is enough time
      char fin_file[1024];
      sprintf(fin_file,"%s/finished",outfolder);
      file_touch(fin_file);
      
      nanalyzed_conf++;
      enough_time=check_remaining_time();
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
