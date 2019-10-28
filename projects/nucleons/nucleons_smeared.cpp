#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "nissa.hpp"

using namespace nissa;

typedef spinspin sss[4];
typedef sss ssss[4];
typedef ssss cssss[3];
typedef cssss diquark[3];
typedef ssss sssss[4];
typedef sssss ssssss[4];

//configuration
int nconf;
double put_theta[4]={0,0,0,0};
char **conf_path,**out_path;
quad_su3 *conf,*smea_conf;
int nmass,im_3pts;
double *mass;
double kappa;
clover_term_t *Cl;

//source
coords source_pos;
spincolor *source,*temp_source;
su3spinspin *original_source;
su3spinspin *seq_source;

//the propagators
su3spinspin ***S0_SL,***S0_SS;
su3spinspin **S1;

//cgm inverter spinors and parameters
spincolor **solDD,*sol_reco[2];
double *stopping_residues;
int niter_max=100000;

//smearing parameters
double Gauss_kappa,ape_alpha;
int Gauss_niter,ape_niter;

//insertion
int tseparation;
int tsink;

spinspin proj[3]; //projectors over N and N*, and 00 compont of N (in the spinorial representation)
int proj_ind1[8]={0,0,1,1,2,2,3,3},proj_ind2[8]={0,2,1,3,0,2,1,3}; //indices different from 0
int proj_couples[4][2]={{0,5},{2,7},{1,4},{3,6}}; //independent couples
int proj_entr[3][4]={{1,1,-1,-1},{1,1,1,1},{1,0,-1,0}};
spinspin C5; //C*gamma5

//two points contractions
int nproton_2pt_contr=32;      //VA     AV       VV       AA       TT        BB        TB        BT
int list_2pt_op1[32]={5,0,5,0, 1,2,3,4, 6,7,8,9, 1,2,3,4, 6,7,8,9, 10,11,12, 13,14,15, 10,11,12, 13,14,15};
int list_2pt_op2[32]={0,5,5,0, 6,7,8,9, 1,2,3,4, 1,2,3,4, 6,7,8,9, 10,11,12, 13,14,15, 13,14,15, 10,11,12};
int compute_also_SL_2pts;

//which 3 pts compute
int compute_3pts[2][2];

//three point contractions
int nproton_3pt_contr=16;
int list_3pt_op[16]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

//three point chromo contractions
int nproton_3pt_chromo_contr=2;
int list_3pt_chromo_op[2]={0,5};

//                      e_00x   e_01x    e_02x     e_10x    e_11x   e_12x     e_20x   e_21x    e_22x
int epsilon[3][3][3]={{{0,0,0},{0,0,1},{0,-1,0}},{{0,0,-1},{0,0,0},{1,0,0}},{{0,1,0},{-1,0,0},{0,0,0}}};
//same thing, but in a compact way
int eps_pos[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
int eps_sig[6]={1,-1,-1,1,1,-1};
//even more compact
int eps1[3]={2,0,1},eps2[3]={1,2,0};

//timings
double tinv=0,tcontr_2pt=0,tcontr_3pt=0,tot_prog_time=0;

dirac_matr gC;

void put_dirac_matr_into_spinspin(spinspin out,dirac_matr *in)
{
  memset(out,0,sizeof(spinspin));
  
  for(int id1=0;id1<4;id1++)
    {
      int id2=in->pos[id1];
      for(int ri=0;ri<2;ri++)
	out[id1][id2][ri]=in->entr[id1][ri];
    }
}

void initialize_nucleons(char *input_path)
{
  //C5
  complex ima={0,1};
  dirac_matr migC,gC5;
  dirac_prod(&migC,&(base_gamma[2]),&(base_gamma[4]));
  unsafe_dirac_compl_prod(&gC,&migC,ima);
  dirac_prod(&gC5,&gC,&(base_gamma[5]));
  
  put_dirac_matr_into_spinspin(C5,&gC5);
  
  //proj[0] and proj[1]
  for(int nns=0;nns<3;nns++) memset(proj[nns],0,sizeof(spinspin));
  for(int id1=0;id1<4;id1++)
    {
      int id2=base_gamma[4].pos[id1];
      
      proj[0][id1][id1][0]=proj[1][id1][id1][0]=0.5;
      complex_prod_double(proj[0][id1][id2],base_gamma[4].entr[id1],+0.5);
      complex_prod_double(proj[1][id1][id2],base_gamma[4].entr[id1],-0.5);
    }
  
  for(int id1=0;id1<4;id1++) if(id1==0||id1==2) proj[2][id1][id1][0]=0.5;
  proj[2][0][2][0]=proj[2][2][0][0]=-0.5;
  
  open_input(input_path);
  
  // 1) Information about the gauge conf
  
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid
  init_grid(T,L);
  //Allocate the gauge Conf
  conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  smea_conf=nissa_malloc("smea_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  Cl=nissa_malloc("Cl",loc_vol,clover_term_t);
  //Read the gauge conf
  read_str_int("NGaugeConf",&nconf);
  conf_path=(char**)malloc(sizeof(char*)*nconf);
  out_path=(char**)malloc(sizeof(char*)*nconf);
  for(int iconf=0;iconf<nconf;iconf++)
    {
      conf_path[iconf]=(char*)malloc(1024);
      out_path[iconf]=(char*)malloc(1024);
      read_str(conf_path[iconf],1024);
      read_str(out_path[iconf],1024);
    }
  
  //Kappa
  read_str_double("Kappa",&(kappa));
  
  // 2) Source position, masses and smearing parameters
  
  expect_str("SourcePosition");
  master_printf("Source position: ");
  for(int idir=0;idir<NDIM;idir++)
    {
      read_int(&(source_pos[idir]));
      master_printf("%d ",source_pos[idir]);
    }
  master_printf("\n");
  //Smearing parameters
  read_str_double("ApeAlpha",&ape_alpha);
  read_str_int("ApeNiter",&ape_niter);
  read_str_double("GaussKappa",&Gauss_kappa);
  read_str_int("GaussNiter",&Gauss_niter);
  //Mass
  read_list_of_double_pairs("MassResidues",&nmass,&mass,&stopping_residues);
  read_str_int("Ind3ptsMass",&(im_3pts));
  
  // 3) 2pts and insertion info
  //compute also SL 2pts?
  read_str_int("ComputeAlsoSL2pts",&compute_also_SL_2pts);
  //tsink-tsource
  read_str_int("TSeparation",&tseparation);
  tsink=(source_pos[0]+tseparation)%glb_size[0];
  
  // 4) three points computation infos
  read_str_int("Compute3ptsLike0Dislike0",&(compute_3pts[0][0]));
  read_str_int("Compute3ptsLike0Dislike1",&(compute_3pts[0][1]));
  read_str_int("Compute3ptsLike1Dislike0",&(compute_3pts[1][0]));
  read_str_int("Compute3ptsLike1Dislike1",&(compute_3pts[1][1]));
  
  close_input();
  
  ///////////////////// Allocate the various spinors ///////////////////////
  
  original_source=nissa_malloc("original_source",loc_vol,su3spinspin);
  
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  temp_source=nissa_malloc("temp_source",loc_vol+bord_vol,spincolor);
  
  //S0 and similars
  solDD=(spincolor**)malloc(sizeof(spincolor*)*nmass);
  S0_SL=(su3spinspin***)malloc(sizeof(su3spinspin**)*nmass);
  S0_SS=(su3spinspin***)malloc(sizeof(su3spinspin**)*nmass);
  for(int imass=0;imass<nmass;imass++)
    {
      solDD[imass]=nissa_malloc("solDD",loc_vol+bord_vol,spincolor);
      
      //smearead-local spinor
      S0_SL[imass]=(su3spinspin**)malloc(sizeof(su3spinspin*)*2);
      S0_SL[imass][0]=nissa_malloc("S0_SL[X][0]",loc_vol,su3spinspin);
      S0_SL[imass][1]=nissa_malloc("S0_SL[X][1]",loc_vol,su3spinspin);
      
      //smeared-smeared
      S0_SS[imass]=(su3spinspin**)malloc(sizeof(su3spinspin*)*2);
      S0_SS[imass][0]=nissa_malloc("S0_SS[X][0]",loc_vol,su3spinspin);
      S0_SS[imass][1]=nissa_malloc("S0_SS[X][1]",loc_vol,su3spinspin);
    }
  
  sol_reco[0]=nissa_malloc("solution_reco[0]",loc_vol,spincolor);
  sol_reco[1]=nissa_malloc("solution_reco[1]",loc_vol,spincolor);
  
  seq_source=nissa_malloc("seqential_source",loc_vol,su3spinspin);
  
  S1=nissa_malloc("S1",nmass,su3spinspin*);
  for(int imass=0;imass<nmass;imass++) S1[imass]=nissa_malloc("S1",loc_vol,su3spinspin);
}

//read a configuration and put anti-periodic condition at the slice tsource-1
void read_conf_and_put_antiperiodic(quad_su3 *conf,char *conf_path,int tsource)
{
  read_ildg_gauge_conf(conf,conf_path);
  
  //calculate plaquette of original conf
  master_printf("plaq: %+16.016lg\n",global_plaquette_lx_conf(conf));
  
  //calcolate Cl
  chromo_operator(Cl,conf);
  
  //prepared the smeared version and  calculate plaquette
  ape_spatial_smear_conf(smea_conf,conf,ape_alpha,ape_niter);
  master_printf("smeared plaq: %+16.16lg\n",global_plaquette_lx_conf(smea_conf));
  
  //Put the anti-periodic condition on the temporal border
  put_theta[0]=1;
  put_boundaries_conditions(conf,put_theta,1,1);
}

//perform the first inversion to produce the S0 for u and d
THREADABLE_FUNCTION_0ARG(calculate_S0)
{
  GET_THREAD_ID();
  
  for(int ic_sour=0;ic_sour<NCOL;ic_sour++)
    for(int id_sour=0;id_sour<4;id_sour++)
      {
	// =====================================================================
	
	// 1) prepare the source
	master_printf("\n(S0) source index: id=%d, ic=%d\n",id_sour,ic_sour);
	get_spincolor_from_su3spinspin(temp_source,original_source,id_sour,ic_sour);
	
	//smear the source
	gaussian_smearing(source,temp_source,smea_conf,Gauss_kappa,Gauss_niter);
	master_printf(" -> source smeared\n");
	
	//============================================================================
	
	// 2) peform the inversion taking time
	if(IS_MASTER_THREAD) tinv-=take_time();
	inv_tmDQ_cgm(solDD,conf,kappa,mass,nmass,niter_max,stopping_residues,source);
	if(IS_MASTER_THREAD) tinv+=take_time();
	master_printf("inversions finished\n");
	
	//============================================================================
	
	// 3) reconstruct the doublet and smear the sink
	for(int imass=0;imass<nmass;imass++)
	  {
	    //reconstruct the doublet
	    reconstruct_tm_doublet(sol_reco[0],sol_reco[1],conf,kappa,mass[imass],solDD[imass]);
	    
	    //convert the id-th spincolor into the colorspinspin and prepare the sink smeared version
	    for(int r=0;r<2;r++)
	      {
		NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
		  {
		    //put the anti-periodic condition on the propagator
		    int dt=glb_coord_of_loclx[ivol][0]-source_pos[0];
		    double arg=M_PI*dt/glb_size[0];
		    complex phase={cos(arg),sin(arg)};
		    
		    unsafe_spincolor_prod_complex(temp_source[ivol],sol_reco[r][ivol],phase);
		    put_spincolor_into_su3spinspin(S0_SL[imass][r][ivol],temp_source[ivol],id_sour,ic_sour);
		  }
		NISSA_PARALLEL_LOOP_END;
		set_borders_invalid(S0_SL[imass][r]);
		
		//smear the sink
		gaussian_smearing(source,temp_source,smea_conf,Gauss_kappa,Gauss_niter);
		NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
		  put_spincolor_into_su3spinspin(S0_SS[imass][r][ivol],source[ivol],id_sour,ic_sour);
		NISSA_PARALLEL_LOOP_END;
		set_borders_invalid(S0_SS[imass][r]);
	      }
	  }
      }
  
  //put the (1+ig5)/sqrt(2) factor
  for(int imass=0;imass<nmass;imass++)
    for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      	for(int ic1=0;ic1<NCOL;ic1++)
	  for(int ic2=0;ic2<NCOL;ic2++)
	    {
	      rotate_spinspin_to_physical_basis(S0_SL[imass][r][ivol][ic1][ic2],!r,!r);
	      rotate_spinspin_to_physical_basis(S0_SS[imass][r][ivol][ic1][ic2],!r,!r);
	    }
  NISSA_PARALLEL_LOOP_END;
  THREAD_BARRIER();
  
  master_printf(" final rotations performed\n");
}
THREADABLE_FUNCTION_END

//Calculate the proton contraction for a single point
THREADABLE_FUNCTION_2ARG(local_diquark, diquark*,diq, su3spinspin*,S)
{
  GET_THREAD_ID();
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int al=0;al<4;al++)
      for(int al1=0;al1<4;al1++)
	for(int ga=0;ga<4;ga++)
	  for(int ga1=0;ga1<4;ga1++)
	    for(int b=0;b<NCOL;b++)
	      for(int b1=0;b1<NCOL;b1++)
		{
		  int a=eps1[b],c=eps2[b],a1=eps1[b1],c1=eps2[b1];
		  //first time reset
		  unsafe_complex_prod  (diq[ivol][b][b1][al][ga][al1][ga1],S[ivol][a1][a][al1][al],S[ivol][c1][c][ga1][ga]);
		  complex_subt_the_prod(diq[ivol][b][b1][al][ga][al1][ga1],S[ivol][a1][c][al1][ga],S[ivol][c1][a][ga1][al]);
		  //both epsilon index (at fixed b) exchanged, so summ
		  complex_summ_the_prod(diq[ivol][b][b1][al][ga][al1][ga1],S[ivol][c1][c][al1][al],S[ivol][a1][a][ga1][ga]);
		  complex_subt_the_prod(diq[ivol][b][b1][al][ga][al1][ga1],S[ivol][c1][a][al1][ga],S[ivol][a1][c][ga1][al]);
		  //now only b indices (a and c) exchanged, so subt
		  complex_subt_the_prod(diq[ivol][b][b1][al][ga][al1][ga1],S[ivol][a1][c][al1][al],S[ivol][c1][a][ga1][ga]);
		  complex_summ_the_prod(diq[ivol][b][b1][al][ga][al1][ga1],S[ivol][a1][a][al1][ga],S[ivol][c1][c][ga1][al]);
		  //again only b1 indices (a1 and c1) exchanged, so subt
		  complex_subt_the_prod(diq[ivol][b][b1][al][ga][al1][ga1],S[ivol][c1][a][al1][al],S[ivol][a1][c][ga1][ga]);
		  complex_summ_the_prod(diq[ivol][b][b1][al][ga][al1][ga1],S[ivol][c1][c][al1][ga],S[ivol][a1][a][ga1][al]);
		}
  NISSA_PARALLEL_LOOP_END;
  THREAD_BARRIER();
}
THREADABLE_FUNCTION_END

//close the diquark with the third propagator
THREADABLE_FUNCTION_3ARG(close_diquark, ssssss*,prot6, diquark*,diq, su3spinspin*,S)
{
  GET_THREAD_ID();
  
  ssssss *loc_prot6=new ssssss[glb_size[0]];
  memset(loc_prot6,0,sizeof(ssssss)*glb_size[0]);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int al=0;al<4;al++)
      for(int be=0;be<4;be++)
	for(int ga=0;ga<4;ga++)
	  for(int al1=0;al1<4;al1++)
	    for(int be1=0;be1<4;be1++)
	      for(int ga1=0;ga1<4;ga1++)
		for(int b=0;b<NCOL;b++)
		  for(int b1=0;b1<NCOL;b1++)
		    complex_summ_the_prod(loc_prot6[glb_coord_of_loclx[ivol][0]][al][be][ga][al1][be1][ga1],S[ivol][b1][b][be1][be],diq[ivol][b][b1][al][ga][al1][ga1]);
  NISSA_PARALLEL_LOOP_END;
  THREAD_BARRIER();
  int ndoub=sizeof(ssssss)/sizeof(double)*glb_size[0];
  glb_threads_reduce_double_vect((double*)loc_prot6,ndoub);
  if(IS_MASTER_THREAD) MPI_Reduce(loc_prot6,prot6,ndoub,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  delete[] loc_prot6;
}
THREADABLE_FUNCTION_END

//calculate all the 2pts contractions
void calculate_all_2pts(char *path,su3spinspin ***S0)
{
  FILE *output=open_text_file_for_output(path);
  
  //take initial time
  tcontr_2pt-=take_time();
  
  //tag for the contraction file
  char pm_tag[2][2]={"+","-"};
  
  //diquark term and all-open dirac indices term
  diquark *diq=nissa_malloc("diq",loc_vol,diquark);
  ssssss prot6[glb_size[0]];
  ssss prot4[glb_size[0]];
  spin prot[glb_size[0]];
  
  //contraction
  int ncontr=nproton_2pt_contr;
  
  //list of gamma
  dirac_matr o1[ncontr],o2[ncontr],o3;
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      dirac_prod(&(o1[icontr]),&gC,&base_gamma[list_2pt_op1[icontr]]);
      o2[icontr]=base_gamma[list_2pt_op2[icontr]];
    }
  dirac_prod(&o3,&gC,&base_gamma[5]);
  
  //loop over like masses
  for(int im_like=0;im_like<nmass;im_like++)
    for(int rlike=0;rlike<2;rlike++)
      {
	//calculate the di-quark part
	local_diquark(diq,S0[im_like][rlike]);
	
	//now close with the third propagator leaving all dirac indices open, and makes global reduction
	for(int im_dislike=0;im_dislike<nmass;im_dislike++)
	  for(int rdislike=0;rdislike<2;rdislike++)
	    {
	      master_fprintf(output,"  # Two points for m_like=%g rlike=%d, m_dislike=%g rdislike=%d\n",
			     mass[im_like],rlike,mass[im_dislike],rdislike);
	      
	      close_diquark(prot6,diq,S0[im_dislike][rdislike]);
	      
	      //now we still have a structure with 6 dirac indices, and is possible to factorize 2 dirac contractions
	      //perform the proton contraction putting operators on the sink or on the source
	      for(int SS=0;SS<2;SS++)
		{
		  memset(prot4,0,sizeof(ssss)*glb_size[0]);
		  
		  //put the gamma3 (the gamma4 is identity!)
		  for(int t=0;t<glb_size[0];t++)
		    for(int al=0;al<4;al++)
		      for(int be=0;be<4;be++)
			for(int ga=0;ga<4;ga++)
			  for(int ga1=0;ga1<4;ga1++)
			    for(int al1=0;al1<4;al1++)
			      {
				int be1=o3.pos[al1];
				if(SS==0) complex_summ_the_prod(prot4[t][al][be][ga][ga1],o3.entr[al1],prot6[t][al][be][ga][al1][be1][ga1]);
				else      complex_summ_the_prod(prot4[t][al][be][ga][ga1],o3.entr[al1],prot6[t][al1][be1][ga1][al][be][ga]);
			      }
		  
		  //now put the other two dirac matrices, remaining with 2 open dirac indices.
		  //We collect them in an appropriate structures:
		  //
		  //   + 0 - 0     + 0 + 0     + 0 + 0  	This are (2 times) the 3 projectors. Only 8 entries
		  //   0 + 0 -     0 + 0 +     0 0 0 0  	are different from 0. Entries are numbered seq.lly,
		  //   - 0 + 0     + 0 + 0     + 0 + 0  	and summed in pairs.
		  //   0 - 0 +     0 + 0 +     0 0 0 0
		  
		  for(int icontr=0;icontr<nproton_2pt_contr;icontr++)
		    {
		      memset(prot,0,sizeof(spin)*glb_size[0]);
		      
		      if(SS==0) master_fprintf(output,"   # %s%s-P5S0\n",gtag[list_2pt_op1[icontr]],gtag[list_2pt_op2[icontr]]);
		      else      master_fprintf(output,"   # P5S0-%s%s\n",gtag[list_2pt_op1[icontr]],gtag[list_2pt_op2[icontr]]);
		      
		      for(int t=0;t<glb_size[0];t++)
			for(int ind_p=0;ind_p<4;ind_p++)
			  for(int ico=0;ico<2;ico++)
			    {
			      int ip=proj_couples[ind_p][ico];
			      
			      int ga=proj_ind1[ip];
			      int ga1=proj_ind2[ip];
			      
			      for(int al=0;al<4;al++)
				{
				  int be=o1[icontr].pos[al],de=o2[icontr].pos[ga];
				  complex fact;
				  
				  unsafe_complex_prod(fact,o1[icontr].entr[al],o2[icontr].entr[ga]);
				  complex_summ_the_prod(prot[t][ind_p],fact,prot4[t][al][be][de][ga1]);
				}
			    }
		      
		      //ok, we managed to have something with 2 left dirac indices (altough organized)
		      //still we need to project them, using the already defined 3 projectors
		      for(int nns=0;nns<3;nns++)
			{
			  if(nns<2) master_fprintf(output,"# Contraction with (1%sg4)/2\n",pm_tag[nns]);
			  else      master_fprintf(output,"# Contraction with (1+g4)_00(spinorial)/2\n");
			      
			  for(int tt=0;tt<glb_size[0];tt++)
			    {
			      int t=(tt+source_pos[0])%glb_size[0];
			      
			      complex contr={0,0};
			      for(int ind_p=0;ind_p<4;ind_p++)
				{
				  if(proj_entr[nns][ind_p]==+1) complex_summ(contr,contr,prot[t][ind_p]);
				  if(proj_entr[nns][ind_p]==-1) complex_subt(contr,contr,prot[t][ind_p]);
				}
			      master_fprintf(output," %+016.16g\t%+016.16g\n",contr[0]/2,contr[1]/2);
			    }
			  master_fprintf(output,"\n");
			}
		    }
		}
	    }
      }
  
  tcontr_2pt+=take_time();
  
  master_printf("contractions finished\n");
  close_file(output);
  
  nissa_free(diq);
}

void prepare_like_sequential_source(int rlike,int rdislike,int slice_to_take)
{
  int dt=(tseparation+source_pos[0])%glb_size[0]-source_pos[0];
  double arg=M_PI*dt/glb_size[0];
  complex phase={cos(arg),sin(arg)};
  
  NISSA_LOC_VOL_LOOP(ivol)
    if(glb_coord_of_loclx[ivol][0]==slice_to_take)
      {
	su3spinspin temp;
	memset(temp,0,sizeof(su3spinspin));
	
	for(int ieps1=0;ieps1<6;ieps1++)
	  {
	    int a1=eps_pos[ieps1][0];
	    int b1=eps_pos[ieps1][1];
	    int c1=eps_pos[ieps1][2];
	    
	    for(int ieps2=0;ieps2<6;ieps2++)
	      {
		int a2=eps_pos[ieps2][0];
		int b2=eps_pos[ieps2][1];
		int c2=eps_pos[ieps2][2];
		
		for(int mu1=0;mu1<4;mu1++)
		  for(int mu2=0;mu2<4;mu2++)
		    for(int be1=0;be1<4;be1++)
		      for(int be2=0;be2<4;be2++)
			for(int al1=0;al1<4;al1++)
			  for(int al2=0;al2<4;al2++)
			    if((C5[be1][mu1][0]||C5[be1][mu1][1])&&(C5[be2][mu2][0]||C5[be2][mu2][1]))
			      {
				int se=eps_sig[ieps1]*eps_sig[ieps2];
				
				complex ter;
				unsafe_complex_prod(ter,S0_SL[im_3pts][rlike][ivol][b1][b2][al1][al2],S0_SL[im_3pts][rlike][ivol][c1][c2][be1][be2]);
				complex_summ_the_prod(ter,S0_SL[im_3pts][rlike][ivol][b1][b2][al1][be2],S0_SL[im_3pts][rlike][ivol][c1][c2][be1][al2]);
				
				safe_complex_prod(ter,C5[be1][mu1],ter);
				safe_complex_prod(ter,C5[be2][mu2],ter);
				safe_complex_prod(ter,proj[2][al2][al1],ter); //create z+ polarized proton
				
				for(int ri=0;ri<2;ri++)
				  if(se==1) temp[a2][a1][mu2][mu1][ri]+=ter[ri];
				  else      temp[a2][a1][mu2][mu1][ri]-=ter[ri];
			      }
	      }	
	  }
	
	//remove the anti-periodic condition on the sequential source
	unsafe_su3spinspin_prod_complex(seq_source[ivol],temp,phase);
      }
}

void prepare_dislike_sequential_source(int rlike,int rdislike,int slice_to_take)
{
  int dt=(tseparation+source_pos[0])%glb_size[0]-source_pos[0];
  double arg=M_PI*dt/glb_size[0];
  complex phase={cos(arg),sin(arg)};
  
  NISSA_LOC_VOL_LOOP(ivol)
    if(glb_coord_of_loclx[ivol][0]==slice_to_take)
      {
	su3spinspin temp;
	memset(temp,0,sizeof(su3spinspin));
	
	for(int x=0;x<3;x++)
	  for(int z=0;z<3;z++)
	    for(int rho=0;rho<4;rho++)
	      for(int tau=0;tau<4;tau++)
		{
		  for(int ga1=0;ga1<4;ga1++)
		    for(int ga2=0;ga2<4;ga2++)
		      {
			if(proj[2][ga2][ga1][0]||proj[2][ga2][ga1][1])
			  for(int a1=0;a1<3;a1++)
			    for(int b1=0;b1<3;b1++)
			      for(int c1=0;c1<3;c1++)
				if(epsilon[a1][b1][c1])
				  for(int al1=0;al1<4;al1++)
				    for(int be1=0;be1<4;be1++)
				      if(C5[al1][be1][0]||C5[al1][be1][1])
					for(int be2=0;be2<4;be2++)
					  {
					    //first part
					    if(C5[rho][be2][0]||C5[rho][be2][1])
					      for(int b2=0;b2<3;b2++)
						for(int c2=0;c2<3;c2++)
						  if(epsilon[x][b2][c2])
						    {
						      complex part={0,0};
						      
						      if(a1==z && al1==tau) complex_summ(part,part,S0_SL[im_3pts][rlike][ivol][c1][c2][ga1][ga2]);
						      if(c1==z && ga1==tau) complex_subt(part,part,S0_SL[im_3pts][rlike][ivol][a1][c2][al1][ga2]);
						      
						      safe_complex_prod(part,C5[rho][be2],part);
						      complex_prod_double(part,part,epsilon[x][b2][c2]);
						      safe_complex_prod(part,S0_SL[im_3pts][rdislike][ivol][b1][b2][be1][be2],part);
						      safe_complex_prod(part,proj[2][ga2][ga1],part); //spin z+ polarized proton
						      
						      if(epsilon[a1][b1][c1]==1) complex_summ_the_prod(temp[x][z][rho][tau],part,C5[al1][be1]);
						      else complex_subt_the_prod(temp[x][z][rho][tau],part,C5[al1][be1]);
						    }
					      
					    //second part
					    for(int al2=0;al2<4;al2++)
					      if((C5[al2][be2][0]||C5[al2][be2][1])&&(ga2==rho))
						for(int a2=0;a2<3;a2++)
						  for(int b2=0;b2<3;b2++)
						    if(epsilon[a2][b2][x])
						      {
							complex part={0,0};
							
							if(c1==z && ga1==tau) complex_summ(part,part,S0_SL[im_3pts][rlike][ivol][a1][a2][al1][al2]);
							if(a1==z && al1==tau) complex_subt(part,part,S0_SL[im_3pts][rlike][ivol][c1][a2][ga1][al2]);
							
							safe_complex_prod(part,C5[al2][be2],part);
							complex_prod_double(part,part,epsilon[a2][b2][x]);
							safe_complex_prod(part,S0_SL[im_3pts][rdislike][ivol][b1][b2][be1][be2],part);
							safe_complex_prod(part,proj[2][ga2][ga1],part);
							
							if(epsilon[a1][b1][c1]==1) complex_summ_the_prod(temp[x][z][rho][tau],part,C5[al1][be1]);
							else complex_subt_the_prod(temp[x][z][rho][tau],part,C5[al1][be1]);
						      }
					  }
		      }
		}
	//remove the anti-periodic condition on the sequential source
	unsafe_su3spinspin_prod_complex(seq_source[ivol],temp,phase);
      }
}

//perform the sequential inversion, on like or dislike according to "ld"
void calculate_S1_like_dislike(int rlike,int rdislike,int ld)
{
  char tag[2][1024]={"like","dislike"};
  
  for(int id_sink=0;id_sink<4;id_sink++)
    for(int ic_sink=0;ic_sink<3;ic_sink++)
      { //take the source and put g5 on the source, and pre-rotate
	NISSA_LOC_VOL_LOOP(ivol)
	  {
	    for(int id_sour=0;id_sour<4;id_sour++)
	      for(int ic_sour=0;ic_sour<3;ic_sour++)
		for(int ri=0;ri<2;ri++)
		  if(id_sour<2)source[ivol][id_sour][ic_sour][ri]= seq_source[ivol][ic_sink][ic_sour][id_sink][id_sour][ri];
		  else         source[ivol][id_sour][ic_sour][ri]=-seq_source[ivol][ic_sink][ic_sour][id_sink][id_sour][ri];
	    if((ld==0 && rdislike==0)||(ld==1 && rlike==0))safe_spincolor_prod_dirac(source[ivol],source[ivol],&Pminus);//up
	    else                                           safe_spincolor_prod_dirac(source[ivol],source[ivol],&Pplus); //dw
	  }
	set_borders_invalid(source);
	
	master_printf("\n(S1) %s, rlike=%d rdislike=%d, sink index: id=%d, ic=%d\n",tag[ld],rlike,rdislike,id_sink,ic_sink);
	
	tinv-=take_time();
	inv_tmQ2_left_cgm(solDD,conf,kappa,mass,nmass,niter_max,stopping_residues,source);
	tinv+=take_time();
	
	//use sol_reco[0] as temporary storage. We solve for rdislike
	for(int imass=0;imass<nmass;imass++)
	  {
	    if((ld==0 && rdislike==0)||(ld==1 && rlike==1))
	      {
		apply_tmQ_left(sol_reco[0],conf,kappa,+mass[imass],solDD[imass]); //cancel d
		NISSA_LOC_VOL_LOOP(ivol) safe_spincolor_prod_dirac(sol_reco[0][ivol],sol_reco[0][ivol],&Pminus);
	      }
	    else
	      {
		apply_tmQ_left(sol_reco[0],conf,kappa,-mass[imass],solDD[imass]); //cancel u
		NISSA_LOC_VOL_LOOP(ivol) safe_spincolor_prod_dirac(sol_reco[0][ivol],sol_reco[0][ivol],&Pplus);
	      }
	    
	    NISSA_LOC_VOL_LOOP(ivol)
	      {
		//put the anti-periodic condition on the propagator
		int dt=glb_coord_of_loclx[ivol][0]-source_pos[0];
		double arg=-M_PI*dt/glb_size[0];
		complex phase={cos(arg),sin(arg)};
		
		spincolor temp;
		
		unsafe_spincolor_prod_complex(temp,sol_reco[0][ivol],phase);
		
		for(int id_sour=0;id_sour<4;id_sour++)
		  for(int ic_sour=0;ic_sour<3;ic_sour++)
		    for(int ri=0;ri<2;ri++)
		      S1[imass][ivol][ic_sink][ic_sour][id_sink][id_sour][ri]=temp[id_sour][ic_sour][ri];
	      }
	  }
      }
  
  master_printf("%s sequential inversions finished\n",tag[ld]);
}

//wrappers
void calculate_S1_like(int rlike,int rdislike)
{calculate_S1_like_dislike(rlike,rdislike,0);}
void calculate_S1_dislike(int rlike,int rdislike)
{calculate_S1_like_dislike(rlike,rdislike,1);}

//this is needed to check 2pts
void contract_with_source(complex *glb_contr,su3spinspin *eta,su3spinspin *S)
{
  complex *loc_contr=(complex*)malloc(sizeof(complex)*glb_size[0]);
  memset(loc_contr,0,sizeof(complex)*glb_size[0]);
  
  NISSA_LOC_VOL_LOOP(ivol)
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++) //eta+*S
	    complex_summ_the_conj1_prod(loc_contr[glb_coord_of_loclx[ivol][0]],eta[ivol][ic2][ic1][id2][id1],S[ivol][ic2][ic1][id2][id1]);
  
  MPI_Reduce(loc_contr,glb_contr,2*glb_size[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  free(loc_contr);
}

//check two points with the current sequential source
void check_2pts_with_current_sequential_source(char *path)
{
  FILE *fout=open_text_file_for_output(path);
  complex *contr_2pts=(complex*)malloc(sizeof(complex)*glb_size[0]);
  
  for(int imass=0;imass<nmass;imass++)
    {
      contract_with_source(contr_2pts,original_source,S1[imass]);
      if(rank==0) fprintf(fout,"mseq=%g %+016.16g\t%+016.16g\n",
			  mass[imass],contr_2pts[source_pos[0]][0],contr_2pts[source_pos[0]][1]);
    }
  
  if(rank==0) fclose(fout);
  
  free(contr_2pts);
}

//Calculate the proton contraction with the inserction of a gamma
void point_proton_sequential_contraction(complex contr,su3spinspin S0,dirac_matr g,su3spinspin S1)
{
  contr[0]=contr[1]=0;
  
  complex temp;
  
  for(int ic1=0;ic1<3;ic1++)
    for(int ic2=0;ic2<3;ic2++)
      for(int id1=0;id1<4;id1++)
	for(int id2=0;id2<4;id2++)
	  {
	    int idg=g.pos[id2];
	    unsafe_complex_prod(temp,        S0[ic2][ic1][id2][id1],g.entr[id2]);
	    complex_summ_the_prod(contr,temp,S1[ic1][ic2][id1][idg]);
	  }
}

//Apply the dipole operator on a su3spinspin
void apply_dipole_operator(su3spinspin *S_out,su3spinspin *S_in,int dir)
{
  NISSA_LOC_VOL_LOOP(ivol)
    {
      int coor=(glb_coord_of_loclx[ivol][dir]-source_pos[dir]+glb_size[dir])%glb_size[dir];
      
      if(coor> glb_size[dir]/2) coor-=glb_size[dir]; //minor distance
      if(coor==glb_size[dir]/2) coor=0; //take away the border (Simpson)
      
      for(int icso=0;icso<3;icso++)
	for(int icsi=0;icsi<3;icsi++)
	  for(int idso=0;idso<4;idso++)
	    for(int idsi=0;idsi<4;idsi++)
	      for(int ri=0;ri<2;ri++)
		S_out[ivol][icsi][icso][idsi][idso][ri]=S_in[ivol][icsi][icso][idsi][idso][ri]*coor;
    }
}

//calculate all the 3pts contractions
void calculate_all_3pts_with_current_sequential(int rlike,int rdislike,int rS0,char *path,const char *tag)
{
  int ncontr;
  FILE *fout=open_text_file_for_output(path);
  su3spinspin *supp_S=nissa_malloc("suppS",loc_vol,su3spinspin);
  
  tcontr_3pt-=take_time();
  
  complex *loc_contr=(complex*)malloc(sizeof(complex)*glb_size[0]);
  complex *glb_contr=(complex*)malloc(sizeof(complex)*glb_size[0]);
  
  for(int im_seq=0;im_seq<nmass;im_seq++)
    for(int im_close=0;im_close<nmass;im_close++)
      {
	
	//this loop is made in order to avoid duplicating the routine
	for(int norm_chro_EDM=0;norm_chro_EDM<3;norm_chro_EDM++)
	  {
	    //switch the contraction
	    switch(norm_chro_EDM)
	      {
	      case 0:ncontr=nproton_3pt_contr;break;
	      case 1:
		ncontr=nproton_3pt_chromo_contr;
		unsafe_apply_chromo_operator_to_su3spinspin(supp_S,Cl,S0_SL[im_close][rS0]);
		break;
	      case 2:
		ncontr=3; //three spatial directions
		//the Dipole Operator is applied inside the contraction (that is, direction) loop
		break;
	      }
	    
	    //loop over contraction
	    for(int icontr=0;icontr<ncontr;icontr++)
	      {
		//reset output
		memset(loc_contr,0,sizeof(complex)*glb_size[0]);
		
		//apply the EDM in the dir=icontr+1
		if(norm_chro_EDM==2) apply_dipole_operator(supp_S,S0_SL[im_close][rS0],icontr+1);
		
		//loop over local node volume
		NISSA_LOC_VOL_LOOP(ivol)
		  {
		    complex point_contr;
		    int glb_t=glb_coord_of_loclx[ivol][0];
		    
		    //contract the single point
		    switch(norm_chro_EDM)
		      {
		      case 0:point_proton_sequential_contraction(point_contr,S0_SL[im_close][rS0][ivol],base_gamma[list_3pt_op[icontr]],S1[im_seq][ivol]);break;
		      case 1:point_proton_sequential_contraction(point_contr,supp_S[ivol],base_gamma[list_3pt_chromo_op[icontr]],S1[im_seq][ivol]);break;
		      case 2:point_proton_sequential_contraction(point_contr,supp_S[ivol],base_gamma[4],S1[im_seq][ivol]);break;
		      }
		    
		    complex_summ(loc_contr[glb_t],loc_contr[glb_t],point_contr);
		  }
		
		//final reduction
		MPI_Reduce(loc_contr,glb_contr,2*glb_size[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		
		if(rank==0)
		  {
		    fprintf(fout," # Three point for rlike=%d, rdislike=%d, %s source, m_spe=%g m_seq=%g m_close=%g\n",
			    rlike,rdislike,tag,mass[im_3pts],mass[im_seq],mass[im_close]);
		    
		    switch(norm_chro_EDM)
		      {
		      case 0:fprintf(fout," # Proton-%s-Proton\n",gtag[list_3pt_op[icontr]]);break;
		      case 1:fprintf(fout," # Proton-CHROMO_%s-Proton\n",gtag[list_3pt_op[icontr]]);break;
		      case 2:fprintf(fout," # Proton-EDM_%d-Proton\n",icontr+1);break;
		      }
		    
		    //print the contraction
		    for(int tt=0;tt<glb_size[0];tt++)
		      {
			int t=(tt+source_pos[0])%glb_size[0];
			fprintf(fout," %+016.16g\t%+016.16g\n",glb_contr[t][0],glb_contr[t][1]);
		      }
		    fprintf(fout,"\n");
		  }
	      }
	  }
      }
  
  tcontr_3pt+=take_time();
  
  if(rank==0)
    {
      fclose(fout);
      printf("contractions ultimated\n");
    }
  
  free(loc_contr);
  free(glb_contr);
  nissa_free(supp_S);
}

void close_nucleons()
{
  nissa_free(conf);nissa_free(smea_conf);nissa_free(Cl);
  nissa_free(original_source);nissa_free(source);
  nissa_free(temp_source);nissa_free(seq_source);
  for(int imass=0;imass<nmass;imass++)
    {
      nissa_free(solDD[imass]);
      nissa_free(S0_SL[imass][0]);nissa_free(S0_SL[imass][1]);
      nissa_free(S0_SS[imass][0]);nissa_free(S0_SS[imass][1]);
      nissa_free(S1[imass]);
    }
  nissa_free(S1);
  nissa_free(sol_reco[0]);nissa_free(sol_reco[1]);
}

void in_main(int narg,char **arg)
{
  tot_prog_time-=take_time();
  
  if(narg<2) crash("Use: %s input_file\n",arg[0]);
  
  initialize_nucleons(arg[1]);
  
  ///////////////////////////////////////////
  
  generate_delta_source(original_source,source_pos);
  
  //loop over configurations
  for(int iconf=0;iconf<nconf;iconf++)
    {
      //prepare the SL output path
      char path_SL[1024],path_SS[1024];
      sprintf(path_SL,"%s_SL",out_path[iconf]);
      sprintf(path_SS,"%s_SS",out_path[iconf]);
      
      read_conf_and_put_antiperiodic(conf,conf_path[iconf],source_pos[0]);
      
      calculate_S0();
      
      //now print both SL and SS
      if(compute_also_SL_2pts) calculate_all_2pts(path_SL,S0_SL);
      calculate_all_2pts(path_SS,S0_SS);
      
      //time for the three points
      for(int rlike=0;rlike<2;rlike++)
	for(int rdislike=0;rdislike<2;rdislike++)
	  if(compute_3pts[rlike][rdislike])
	    {
	      char out2pts_check_like[1024],out3pts_like[1024];
	      char out2pts_check_dislike[1024],out3pts_dislike[1024];
	      
	      sprintf(out2pts_check_like,"2pts_check_like%d%d%d",rlike,rlike,rdislike);
	      sprintf(out2pts_check_dislike,"2pts_check_dislike%d%d%d",rlike,rlike,rdislike);
	      sprintf(out3pts_like,"3pts_like%d%d%d",rlike,rlike,rdislike);
	      sprintf(out3pts_dislike,"3pts_dislike%d%d%d",rlike,rlike,rdislike);
	      
	      //like three points
	      prepare_like_sequential_source(rlike,rdislike,tsink);
	      calculate_S1_like(rlike,rdislike);
	      check_2pts_with_current_sequential_source(out2pts_check_like);
	      calculate_all_3pts_with_current_sequential(rlike,rdislike,rdislike,out3pts_like,"like");
	      
	      //dislike three points
	      prepare_dislike_sequential_source(rlike,rdislike,tsink);
	      calculate_S1_dislike(rlike,rdislike);
	      check_2pts_with_current_sequential_source(out2pts_check_dislike);
	      calculate_all_3pts_with_current_sequential(rlike,rdislike,rlike,out3pts_dislike,"dislike");
	    }
      
    }
  
  tot_prog_time+=take_time();
  
  ///////////////////////////////////////////
  
  master_printf("Total time: %g s\n",tot_prog_time);
  master_printf("-inversion time: %g%s avg: %d s\n",tinv/tot_prog_time*100,"%",(int)(tinv/nconf/12));
  master_printf("-contraction time for 2pts: %g%s\n",tcontr_2pt/tot_prog_time*100,"%");
  master_printf("-contraction time for 3pts: %g%s\n",tcontr_3pt/tot_prog_time*100,"%");
  
  close_nucleons();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
