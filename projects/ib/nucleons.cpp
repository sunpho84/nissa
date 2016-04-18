#include <nissa.hpp>

#include "conf.hpp"
#include "pars.hpp"
#include "prop.hpp"

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

//hadron contractions
struct bar_triplet_t
{
  int a,b,c;
  bar_triplet_t(int a,int b,int c) : a(a),b(b),c(c) {}
};
std::vector<bar_triplet_t> prop_hadr_combo_map;

//index inside a colorspinspin
int dirspin_ind(int ic_si,int ic_so,int id_si,int id_so,int ri)
{return
    (ri+2*
     (id_so+4*
      (id_si+4*
       (ic_so+NCOL*ic_si))));
}
int ind_rot_prop(int iprop,int ci,int ivol)
{return
    (ivol+loc_vol*
     (ci+sizeof(su3spinspin)/sizeof(double)*iprop));
}
int ind_rot_prop(int iprop,int ic_si,int ic_so,int id_si,int id_so,int ri,int ivol)
{return ind_rot_prop(iprop,dirspin_ind(ic_si,ic_so,id_si,id_so,ri),ivol);}

//contains how to bind the indices of the three propagators to make a Wick contraction
int eps_perm[6][4]={{0,1,2,+1},{1,2,0,+1},{2,0,1,+1},
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
    for(int iperm_si=0;iperm_si<6;iperm_si++)
      for(int iperm_so=0;iperm_so<6;iperm_so++)
	for(int al1=0;al1<4;al1++)
	  for(int al=0;al<4;al++)
	    for(int ga=0;ga<4;ga++)
	      {
		//find other indices
		int a =eps_perm[iperm_si][0],b =eps_perm[iperm_si][1],c =eps_perm[iperm_si][2];
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
			for(int dir_exc=0;dir_exc<2;dir_exc++)
			  {
			    const int itb=dirspin_ind(b1,b,be1,be,rib);
			    const int ita=dirspin_ind(a1,(dir_exc==0)?a:c,al1,(dir_exc==0)?al:ga,ria);
			    const int itc=dirspin_ind(c1,(dir_exc==0)?c:a,ga1,(dir_exc==0)?ga:al,ric);
			    
			    //find whether the real or imaginary part contributes
			    for(int ri=0;ri<2;ri++)
			      if(w[ri])
				wickes[ind_wick(dir_exc,ri,iid_g0)].push_back(std::make_pair(bar_triplet_t(ita,itb,itc),w[ri]));
			  }
		      }
	      }
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
  //prop_hadr_combo_map.push_back(std::make_pair(PROP_0,PROP_S));
  //prop_hadr_combo_map.push_back(std::make_pair(PROP_0,PROP_P));
  //prop_hadr_combo_map.push_back(std::make_pair(PROP_0,PROP_T));
  //prop_hadr_combo_map.push_back(std::make_pair(PROP_0,PROP_PHOTON2));
  //prop_hadr_combo_map.push_back(std::make_pair(PROP_PHOTON,PROP_PHOTON));
  
  //prop_hadr_combo_map.push_back(std::make_pair(PROP_0,PROP_VECTOR));
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
      if(ism_sink==1)
	for(int ip=0;ip<nqprop;ip++)
	gaussian_smearing(Q[ip],Q[ip],ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
      
      //transpose all propagators
      for(int ip=0;ip<nqprop;ip++)
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  for(int ic_si=0;ic_si<NCOL;ic_si++)
	    for(int ic_so=0;ic_so<NCOL;ic_so++)
	      for(int id_si=0;id_si<4;id_si++)
		for(int id_so=0;id_so<4;id_so++)
		  for(int ri=0;ri<2;ri++)
		    p[ind_rot_prop(ip,ic_si,ic_so,id_si,id_so,ri,ivol)]=Q[ip][ivol][ic_si][ic_so][id_si][id_so][ri];
      THREAD_BARRIER();
      
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
			    
			    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
			    loc_corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,iwick,glb_coord_of_loclx[ivol][0])]+=
			      wickes[iwick][iterm].second*
			      p[ind_rot_prop(ipa,iwa,ivol)]*
			      p[ind_rot_prop(ipb,iwb,ivol)]*
			      p[ind_rot_prop(ipc,iwc,ivol)];
			  }
		    }
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
  glb_nodes_reduce_double_vect(corr,corr_size);
  
  //phase
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
		      master_printf("\n # icombo %lu , sink_smeared %d , im %d %d %d , r %d %d %d , dir_exc %d\n",icombo,ism_sink,ima,imb,imc,ra,rb,rc,dir_exc);
		      for(int t=0;t<glb_size[0];t++)
			{
			  //reassembling real/imaginary
			  //reassembling (1+-g0)/2 depending on t<>T/2
			  //putting factor 1/nsources
			  int sign=(t<glb_size[0]/2)?1:-1;
			  complex c;
			  
			  for(int ri=0;ri<2;ri++)
			    {
			      //take id and g0
			      int iwick_id=ind_wick(dir_exc,ri,0/*id*/);
			      int iwick_g0=ind_wick(dir_exc,ri,1/*g0*/);
			      
			      c[ri]=
				(corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,iwick_id,t)]+sign*
				 corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,iwick_g0,t)])/
				(2*nsources);
			    }
			  
			  master_printf("%+016.016lg %+016.016lg\n",c[RE],c[IM]);
			}
		      master_printf("\n");
		    }
}

//close deallocating everything
void close()
{
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
	  
	  random_shift_gauge_conf(conf,old_theta,put_theta);
	  generate_photon_stochastic_propagator();
	  generate_original_source();
	  gaussian_smearing(original_source,original_source,ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
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
