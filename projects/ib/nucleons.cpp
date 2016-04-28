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

int nsmear,ncontract;
double smear_time=0,contract_time=0;

//hadron contractions
struct bar_triplet_t
{
  int a,b,c;
  bar_triplet_t(int a,int b,int c) : a(a),b(b),c(c) {}
};
std::vector<bar_triplet_t> prop_hadr_combo_map;

//index inside a colorspinspin
inline int dirspin_ind(int ic_si,int ic_so,int id_si,int id_so,int ri)
{return
    (ri+2*
     (id_so+NDIRAC*
      (id_si+NDIRAC*
       (ic_so+NCOL*ic_si))));
}
inline int ind_rot_prop(int iprop,int ci,int ivol)
{return
    (ivol+loc_vol*
     (ci+sizeof(su3spinspin)/sizeof(double)*iprop));
}
inline int ind_rot_prop(int iprop,int ic_si,int ic_so,int id_si,int id_so,int ri,int ivol)
{return ind_rot_prop(iprop,dirspin_ind(ic_si,ic_so,id_si,id_so,ri),ivol);}

dirac_matr Cg5;
void set_wickes_contractions()
{
  //set C*g5
  dirac_matr g2g4,C;
  dirac_prod(&g2g4,base_gamma+2,base_gamma+4);
  dirac_prod_idouble(&C,&g2g4,1);
  dirac_prod(&Cg5,&C,base_gamma+5);
}

//store correlation function
int nsm_sink=2;
complex *corr;
int ind_corr(int icombo,int ism_sink,int ima,int ra,int imb,int rb,int imc,int rc,int dir_exc,int t)
{return
    (t+glb_size[0]*
     (dir_exc+2*
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
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_PHOTON_AB,PROP_0,PROP_0));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_PHOTON_AB,PROP_0));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_0,PROP_PHOTON_AB));
  //photon exchange between one of the three lines
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_0,PROP_PHOTON_A,PROP_PHOTON_B));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_PHOTON_A,PROP_0,PROP_PHOTON_B));
  prop_hadr_combo_map.push_back(bar_triplet_t(PROP_PHOTON_A,PROP_PHOTON_B,PROP_0));
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
  read_use_photon_field();
  read_use_photon_field();
  read_loc_hadr_curr();
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
  photon_eta=nissa_malloc("photon_eta",loc_vol+bord_vol,spin1field);
  photon_field=nissa_malloc("photon_field",loc_vol+bord_vol,spin1field);
  photon_phi=nissa_malloc("photon_phi",loc_vol+bord_vol,spin1field);
  source=nissa_malloc("source",loc_vol,su3spinspin);
  original_source=nissa_malloc("original_source",loc_vol,su3spinspin);
  corr_size=ind_corr(prop_hadr_combo_map.size()-1,nsm_sink-1,nqmass-1,nr-1,nqmass-1,nr-1,nqmass-1,nr-1,2-1,glb_size[0]-1)+1;
  corr=nissa_malloc("corr",corr_size,complex);
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
  complex *loc_corr=new complex[corr_size];
  memset(loc_corr,0,sizeof(complex)*corr_size);
  
  for(int ism_sink=0;ism_sink<nsm_sink;ism_sink++)
    {
      //smear all sinks
      START_TIMING(smear_time,nsmear);
      if(ism_sink)
	for(int ip=0;ip<nqprop;ip++)
	  gaussian_smearing(Q[ip],Q[ip],ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
      STOP_TIMING(smear_time);
      
      const int eps[3][2]={{1,2},{2,0},{0,1}},sign[2]={1,-1};
      
      void (*list_fun[2])(complex,complex,complex)={complex_summ_the_prod,complex_subt_the_prod};
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
		      
		      if(IS_MASTER_THREAD) ncontract++;
		      
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
				    complex_summ_the_prod(loc_corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,0,t)],Q[ipb][ivol][b1][b][be1][be],diquark_dir);
				    complex_summ_the_prod(loc_corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,1,t)],Q[ipb][ivol][b1][b][be1][be],diquark_exc);
				  }
			}
		    }
      STOP_TIMING(contract_time);
    }
  
  //reduce
  complex *master_reduced_corr=(complex*)glb_threads_reduce_double_vect((double*)loc_corr,corr_size*2);
  NISSA_PARALLEL_LOOP(i,0,corr_size) complex_summassign(corr[i],master_reduced_corr[i]);
  THREAD_BARRIER();
  delete[] loc_corr;
}
THREADABLE_FUNCTION_END

//print all correlations averaging
void print_correlations()
{
  //reduce
  glb_nodes_reduce_complex_vect(corr,corr_size);
  
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
			  //remove anti-periodic condition phase
			  double arg=3*M_PI*t/glb_size[0];
			  complex phase={cos(arg),sin(arg)};
			  
			  //normalize for nsources and 1+g0
			  complex c;
			  complex_prod_double(c,corr[ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,dir_exc,t)],1.0/(2*nsources));
			  
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
  master_printf(" - %02.2f%s to perform %d contractions (%2.2gs avg)\n",contract_time/tot_prog_time*100,"%",ncontract,contract_time/ncontract);
  master_printf(" - %02.2f%s to perform %d smearing (%2.2gs avg)\n",smear_time/tot_prog_time*100,"%",nsmear,smear_time/nsmear);
  
  nissa_free(photon_eta);
  nissa_free(photon_phi);
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
