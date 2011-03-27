#include "appretto.h"

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
as2t_su3 *Pmunu;

//source
int source_pos[4];
spincolor *source,*temp_source;
su3spinspin *original_source;
su3spinspin *seq_source;

//the propagators
su3spinspin ***S0_SL,***S0_SS;
su3spinspin *S1;

//cgmms inverter spinors and parameters
spincolor **solDD,*sol_reco[2];
double stopping_residue;
double minimal_residue;
int stopping_criterion;
int niter_max;

//smearing parameters
double jacobi_kappa,ape_alpha;
int jacobi_niter,ape_niter;

//insertion
int tseparation;
int tsink;

spinspin Proj[3]; //projectors over N and N*, and 00 compont of N (in the spinorial representation)
spinspin C5; //C*gamma5

//two points contractions
int nproton_2pt_contr=32;      //VA     AV       VV       AA       TT        BB        TB        BT
int list_2pt_op1[32]={5,0,5,0, 1,2,3,4, 6,7,8,9, 1,2,3,4, 6,7,8,9, 10,11,12, 13,14,15, 10,11,12, 13,14,15};
int list_2pt_op2[32]={0,5,5,0, 6,7,8,9, 1,2,3,4, 1,2,3,4, 6,7,8,9, 10,11,12, 13,14,15, 13,14,15, 10,11,12};

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
int eps_pos[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
int eps_sig[6]={1,-1,-1,1,1,-1};
int eps1[3]={2,0,1},eps2[3]={1,2,0};

//timings
double tinv=0,tcontr=0,tot_time=0;

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
  
  //Proj[0] and Proj[1]
  for(int nns=0;nns<3;nns++) memset(Proj[nns],0,sizeof(spinspin));
  for(int id1=0;id1<4;id1++)
    {
      int id2=base_gamma[4].pos[id1];
      
      Proj[0][id1][id1][0]=Proj[1][id1][id1][0]=0.5;
      complex_prod_with_real(Proj[0][id1][id2],base_gamma[4].entr[id1],+0.5);
      complex_prod_with_real(Proj[1][id1][id2],base_gamma[4].entr[id1],-0.5);
    }

  for(int id1=0;id1<4;id1++) if(id1==0||id1==2) Proj[2][id1][id1][0]=Proj[2][id1][id1][0]=0.5;
  Proj[2][0][2][0]=Proj[2][2][0][0]=-0.5; 

  open_input(input_path);

  // 1) Information about the gauge conf
  
  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  //Init the MPI grid 
  init_grid();
  //Allocate the gauge Conf
  conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"conf");
  smea_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"conf");
  Pmunu=allocate_as2t_su3(loc_vol,"Pmunu");
  //Read the gauge conf
  read_str_int("NGaugeConf",&nconf);
  conf_path=(char**)malloc(sizeof(char*)*nconf);
  out_path=(char**)malloc(sizeof(char*)*nconf);
  for(int iconf=0;iconf<nconf;iconf++)
    {
      conf_path[iconf]=(char*)malloc(1024);
      out_path[iconf]=(char*)malloc(1024);
      read_str_str("GaugeConfPath",conf_path[iconf],1024);
      read_str(out_path[iconf],1024);
    }

  //Kappa
  read_str_double("Kappa",&(kappa));
  
  // 2) Source position, masses and smerding parameters

  expect_str("SourcePosition");
  if(rank==0) printf("Source position: ");
  for(int idir=0;idir<4;idir++) 
    {
      read_int(&(source_pos[idir]));
      if(rank==0) printf("%d ",source_pos[idir]);
    }
  if(rank==0) printf("\n");
  //Smearing parameters
  read_str_double("ApeAlpha",&ape_alpha);
  read_str_int("ApeNiter",&ape_niter);
  read_str_double("JacobiKappa",&jacobi_kappa);
  read_str_int("JacobiNiter",&jacobi_niter);
  //Mass
  read_list_of_doubles("NMass",&nmass,&mass);
  read_str_int("Ind3ptsMass",&(im_3pts));
  
  // 3) inverter

  //Residue
  read_str_double("Residue",&stopping_residue);
  //Stopping criterion
  stopping_criterion=numb_known_stopping_criterion;
  char str_stopping_criterion[1024];
  read_str_str("StoppingCriterion",str_stopping_criterion,1024);
  int isc=0;
  do
    {
      if(strcasecmp(list_known_stopping_criterion[isc],str_stopping_criterion)==0) stopping_criterion=isc;
      isc++;
    }
  while(isc<numb_known_stopping_criterion && stopping_criterion==numb_known_stopping_criterion);
  
  if(stopping_criterion==numb_known_stopping_criterion && rank==0)
    {
      fprintf(stderr,"Unknown stopping criterion: %s\n",str_stopping_criterion);
      fprintf(stderr,"List of known stopping criterions:\n");
      for(int isc=0;isc<numb_known_stopping_criterion;isc++) fprintf(stderr," %s\n",list_known_stopping_criterion[isc]);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  if(stopping_criterion==sc_standard) read_str_double("MinimalResidue",&minimal_residue);
      
  //Number of iterations
  read_str_int("NiterMax",&niter_max);
  
  // 4) insertion info
  //tsink-tsource
  read_str_int("TSeparation",&tseparation);
  tsink=(source_pos[0]+tseparation)%glb_size[0];
  
  // 5) three points computation infos
  read_str_int("Compute3ptsLike0Dislike0",&(compute_3pts[0][0]));
  read_str_int("Compute3ptsLike0Dislike1",&(compute_3pts[0][1]));
  read_str_int("Compute3ptsLike1Dislike0",&(compute_3pts[1][0]));
  read_str_int("Compute3ptsLike1Dislike1",&(compute_3pts[1][1]));
  
  close_input();

  ///////////////////// Allocate the various spinors ///////////////////////
  
  original_source=allocate_su3spinspin(loc_vol,"original_source");
  
  source=allocate_spincolor(loc_vol+loc_bord,"source");
  temp_source=allocate_spincolor(loc_vol+loc_bord,"temp_source");

  //S0 and similars
  solDD=(spincolor**)malloc(sizeof(spincolor*)*nmass);
  S0_SL=(su3spinspin***)malloc(sizeof(su3spinspin**)*nmass);
  S0_SS=(su3spinspin***)malloc(sizeof(su3spinspin**)*nmass);
  for(int imass=0;imass<nmass;imass++)
    {
      solDD[imass]=allocate_spincolor(loc_vol+loc_bord,"solDD");

      //smearead-local spinor
      S0_SL[imass]=(su3spinspin**)malloc(sizeof(su3spinspin*)*2);
      S0_SL[imass][0]=allocate_su3spinspin(loc_vol,"S0_SL[X][0]");
      S0_SL[imass][1]=allocate_su3spinspin(loc_vol,"S0_SL[X][1]");

      //smeared-smeared
      S0_SS[imass]=(su3spinspin**)malloc(sizeof(su3spinspin*)*2);
      S0_SS[imass][0]=allocate_su3spinspin(loc_vol,"S0_SS[X][0]");
      S0_SS[imass][1]=allocate_su3spinspin(loc_vol,"S0_SS[X][1]");
    }

  sol_reco[0]=allocate_spincolor(loc_vol,"solution_reco[0]");
  sol_reco[1]=allocate_spincolor(loc_vol,"solution_reco[1]");

  seq_source=allocate_su3spinspin(loc_vol,"seqential_source");

  S1=allocate_su3spinspin(loc_vol,"S1");
}

//read a configuration and put anti-periodic condition at the slice tsource-1
void read_conf_and_put_antiperiodic(quad_su3 *conf,char *conf_path,int tsource)
{
  read_local_gauge_conf(conf,conf_path);
  
  //commmunicate borders
  communicate_gauge_borders(conf);  
  communicate_gauge_edges(conf);
  
  //calculate plaquette of original conf
  double gplaq=global_plaquette(conf);
  if(rank==0) printf("plaq: %.18g\n",gplaq);

  //calcolate Pmunu
  Pmunu_term(Pmunu,conf);
  
  //prepared the smerded version and  calculate plaquette
  ape_smearing(smea_conf,conf,ape_alpha,ape_niter);
  gplaq=global_plaquette(smea_conf);
  if(rank==0) printf("smerded plaq: %.18g\n",gplaq);
  
  //Put the anti-periodic condition on the temporal border
  put_theta[0]=1;
  put_boundaries_conditions(conf,put_theta,1,1);

  //re-communicate borders
  communicate_gauge_borders(conf);
  communicate_gauge_borders(smea_conf);
  communicate_gauge_edges(conf);
}

//create a 12 index point source
void prepare_source()
{
  int isloc=1;
  int lx[4];

  //reset the source
  memset(original_source,0,sizeof(spincolor)*loc_vol);

  //check if the source position is associated to the rank and calculate its local position
  for(int idir=0;idir<4;idir++)
    {
      lx[idir]=source_pos[idir]-proc_coord[idir]*loc_size[idir];
      isloc=isloc && (lx[idir]>=0 && lx[idir]<loc_size[idir]);
    }
  
  //if the source is in the rank, put it
  int ivol=loclx_of_coord(lx);
  for(int ic=0;ic<3;ic++)
    for(int id=0;id<4;id++)
      if(isloc) original_source[ivol][ic][ic][id][id][0]=1;
}      

//perform the first inversion to produce the S0 for u and d
void calculate_S0()
{
  for(int ic_sour=0;ic_sour<3;ic_sour++)
    for(int id_sour=0;id_sour<4;id_sour++)
      {
	
	// =====================================================================

	// 1) prepare the source
	
	if(rank==0) printf("\n(S0) source index: id=%d, ic=%d\n",id_sour,ic_sour);
	
	//take the source and put g5
	for(int ivol=0;ivol<loc_vol;ivol++)
	  {
	    get_spincolor_from_su3spinspin(temp_source[ivol],original_source[ivol],id_sour,ic_sour);
	    for(int id_sink=2;id_sink<4;id_sink++)
	      for(int ic_sink=0;ic_sink<3;ic_sink++)
		for(int ri=0;ri<2;ri++)
		  temp_source[ivol][id_sink][ic_sink][ri]=-temp_source[ivol][id_sink][ic_sink][ri];
	  }
	
	//smerd the source
	jacobi_smearing(source,temp_source,smea_conf,jacobi_kappa,jacobi_niter);

	//print the denisity profile
        if(ic_sour==0 && id_sour==0)
          {
            int L=glb_size[1];
            double rho[L];
            density_profile(rho,source,source_pos);
            FILE *fout=open_text_file_for_output("profile");
            if(rank==0)
              {
                for(int d=0;d<L;d++) fprintf(fout,"%d %g\n",d,rho[d]);
                fclose(fout);
              }
          }     

	if(rank==0) printf(" -> source smeared\n");
	
	//============================================================================
	
	// 2) peform the inversion taking time

	tinv-=take_time();
	inv_Q2_cgmms(solDD,source,NULL,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	tinv+=take_time();

	if(rank==0) printf("inversions finished\n");
	
	//============================================================================

	// 3) reconstruct the doublet and smerd the sink

	for(int imass=0;imass<nmass;imass++)
	  {
	    //reconstruct the doublet
	    reconstruct_doublet(sol_reco[0],sol_reco[1],solDD[imass],conf,kappa,mass[imass]);
	    
	    //convert the id-th spincolor into the colorspinspin and prepare the sink smerded version
	    for(int r=0;r<2;r++)
	      {
		for(int ivol=0;ivol<loc_vol;ivol++)
		  {
		    //put the anti-periodic condition on the propagator
		    int dt=glb_coord_of_loclx[ivol][0]-source_pos[0];
		    double arg=M_PI*dt/glb_size[0];
		    complex phase={cos(arg),sin(arg)};
		    
		    unsafe_spincolor_prod_complex(temp_source[ivol],sol_reco[r][ivol],phase);
		    put_spincolor_into_su3spinspin(S0_SL[imass][r][ivol],temp_source[ivol],id_sour,ic_sour);
		  }
		
		//smerd the sink
		jacobi_smearing(source,temp_source,smea_conf,jacobi_kappa,jacobi_niter);
		for(int ivol=0;ivol<loc_vol;ivol++)
		  put_spincolor_into_su3spinspin(S0_SS[imass][r][ivol],source[ivol],id_sour,ic_sour);
	      }
	  }
      }

  //put the (1+ig5)/sqrt(2) factor
  for(int imass=0;imass<nmass;imass++)
    for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
      for(int ivol=0;ivol<loc_vol;ivol++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    {
	      rotate_spinspin_to_physical_basis(S0_SL[imass][r][ivol][ic1][ic2],!r,!r);
	      rotate_spinspin_to_physical_basis(S0_SS[imass][r][ivol][ic1][ic2],!r,!r);
	    }
  
  if(rank==0) printf(" final rotations performed\n");
}

//Calculate the proton contraction for a single point
void local_diquark(diquark *diq,su3spinspin *S)
{
  for(int l=0;l<loc_vol;l++)
    for(int al=0;al<4;al++)
      for(int al1=0;al1<4;al1++)
	for(int ga=0;ga<4;ga++)
	  for(int ga1=0;ga1<4;ga1++)
	    for(int b=0;b<3;b++)
	      for(int b1=0;b1<3;b1++)
		{
		  int a=eps1[b],c=eps2[b],a1=eps1[b1],c1=eps2[b1];

		  //first time reset
		  unsafe_complex_prod  (diq[l][b][b1][al][ga][al1][ga1],S[l][a1][a][al1][al],S[l][c1][c][ga1][ga]);
		  complex_subt_the_prod(diq[l][b][b1][al][ga][al1][ga1],S[l][a1][c][al1][ga],S[l][c1][a][ga1][al]);
		  //both epsilon index (at fixed b) exchanged, so summ
		  complex_summ_the_prod(diq[l][b][b1][al][ga][al1][ga1],S[l][c1][c][al1][al],S[l][a1][a][ga1][ga]);
		  complex_subt_the_prod(diq[l][b][b1][al][ga][al1][ga1],S[l][c1][a][al1][ga],S[l][a1][c][ga1][al]);
		  //now only b indexes (a and c) exchanged, so subt
		  complex_subt_the_prod(diq[l][b][b1][al][ga][al1][ga1],S[l][a1][c][al1][al],S[l][c1][a][ga1][ga]);
		  complex_summ_the_prod(diq[l][b][b1][al][ga][al1][ga1],S[l][a1][a][al1][ga],S[l][c1][c][ga1][al]);
		  //again only b1 indexes (a1 and c1) exchanged, so subt
		  complex_subt_the_prod(diq[l][b][b1][al][ga][al1][ga1],S[l][c1][a][al1][al],S[l][a1][c][ga1][ga]);
		  complex_summ_the_prod(diq[l][b][b1][al][ga][al1][ga1],S[l][c1][c][al1][ga],S[l][a1][a][ga1][al]);
		}
}

void close_diquark(ssssss *prot6,diquark *diq,su3spinspin* S)
{
  ssssss *loc_prot6=(ssssss*)malloc(sizeof(ssssss)*glb_size[0]);
  memset(loc_prot6,0,sizeof(ssssss)*glb_size[0]);

  for(int l=0;l<loc_vol;l++)
    {
      int t=glb_coord_of_loclx[l][0];
      
      for(int al=0;al<4;al++)
	for(int be=0;be<4;be++)
	  for(int ga=0;ga<4;ga++)
	    for(int al1=0;al1<4;al1++)
	      for(int be1=0;be1<4;be1++)
		for(int ga1=0;ga1<4;ga1++)
		  for(int b=0;b<3;b++)
		    for(int b1=0;b1<3;b1++)
		      complex_summ_the_prod(loc_prot6[t][al][be][ga][al1][be1][ga1],S[l][b1][b][be1][be],diq[l][b][b1][al][ga][al1][ga1]);
    }
  
  MPI_Reduce(loc_prot6,prot6,glb_size[0]*8192,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  free(loc_prot6);
}

//calculate all the 2pts contractions
void calculate_all_2pts(char *path,su3spinspin ***S0)
{
  //output file
  FILE *output=open_text_file_for_output(path);

  //take initial time
  tcontr-=take_time();
  
  //tag for the contraction file
  char pm_tag[2][2]={"+","-"};
  
  //diquark term and all-open dirac indexes term
  diquark *diq=(diquark*)malloc(sizeof(diquark)*loc_vol);
  ssssss *prot6=(ssssss*)malloc(sizeof(ssssss)*glb_size[0]);
  ssss *prot4=(ssss*)malloc(sizeof(ssss)*glb_size[0]);
  spinspin *prot2=(spinspin*)malloc(sizeof(spinspin)*glb_size[0]);

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
	
	//now close with the third propagator leaving all dirac indexes open, and makes global reduction
	for(int im_dislike=0;im_dislike<nmass;im_dislike++)
	  for(int rdislike=0;rdislike<2;rdislike++)
	    {
	      if(rank==0) fprintf(output,"  # Two points for m_like=%g rlike=%d, m_dislike=%g rdislike=%d\n",
				  mass[im_like],rlike,mass[im_dislike],rdislike);
	      
	      close_diquark(prot6,diq,S0[im_dislike][rdislike]);
	      
	      //now we still have a structure with 6 dirac indexes, and is possible to factorize 2 dirac contractions
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
		  
		  //now put the other two dirac matrices, remaining with 2 open dirac indexes
		  for(int icontr=0;icontr<nproton_2pt_contr;icontr++)
		    {
		      memset(prot2,0,sizeof(spinspin)*glb_size[0]);
		      
		      if(rank==0) //header
			{
			  if(SS==0)fprintf(output,"   # %s%s-P5S0\n",gtag[list_2pt_op1[icontr]],gtag[list_2pt_op2[icontr]]);
			  else     fprintf(output,"   # P5S0-%s%s\n",gtag[list_2pt_op1[icontr]],gtag[list_2pt_op2[icontr]]);
			}

		      for(int t=0;t<glb_size[0];t++)
			for(int ga=0;ga<4;ga++)
			  for(int ga1=0;ga1<4;ga1++)
			    for(int al=0;al<4;al++)
			      {
				int be=o1[icontr].pos[al],de=o2[icontr].pos[ga];
				complex fact;
				
				unsafe_complex_prod(fact,o1[icontr].entr[al],o2[icontr].entr[ga]);
				complex_summ_the_prod(prot2[t][ga][ga1],fact,prot4[t][al][be][de][ga1]);
			      }
		    
		      //ok, we managed to have something with 2 left dirac indexes.
		      //still we need to project them, using the already defined 3 projectors
		      for(int nns=0;nns<3;nns++)
			{
			  if(rank==0)
			    {
			      if(nns<2) fprintf(output,"# Contraction with (1%sg4)/2\n",pm_tag[nns]);
			      else      fprintf(output,"# Contraction with (1+g4)_00(spinorial)/2\n");
			      
			      for(int tt=0;tt<glb_size[0];tt++)
				{
				  int t=(tt+source_pos[0])%glb_size[0];

				  complex contr;
				  trace_prod_spinspins(contr,prot2[t],Proj[nns]);
				  fprintf(output," %+016.16g\t%+016.16g\n",contr[0],contr[1]);
				}
			      fprintf(output,"\n");
			    }
			}
		    }
		}
	    }
      }
  
  tcontr+=take_time();

  if(rank==0)
    {
      printf("contractions finished\n");
      fclose(output);
    }

  free(diq);
  free(prot6);
  free(prot4);
  free(prot2);
}

void prepare_like_sequential_source(int rlike,int rdislike,int slice_to_take)
{
  int dt=(tseparation+source_pos[0])%glb_size[0]-source_pos[0];
  double arg=M_PI*dt/glb_size[0];
  complex phase={cos(arg),sin(arg)};
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    if(glb_coord_of_loclx[ivol][0]==slice_to_take)
      {
	su3spinspin temp;
	memset(temp,0,sizeof(su3spinspin));
	
	for(int a1=0;a1<3;a1++)
	  for(int b1=0;b1<3;b1++)
	    for(int c1=0;c1<3;c1++)
	      if(epsilon[a1][b1][c1])
		for(int a2=0;a2<3;a2++)
		  for(int b2=0;b2<3;b2++)
		    for(int c2=0;c2<3;c2++)
		      if(epsilon[a2][b2][c2])
			for(int mu1=0;mu1<4;mu1++)
			  for(int mu2=0;mu2<4;mu2++)
			    for(int be1=0;be1<4;be1++)
			      for(int be2=0;be2<4;be2++)
				for(int al1=0;al1<4;al1++)
				  for(int al2=0;al2<4;al2++)
				    if((C5[be1][mu1][0]||C5[be1][mu1][1])&&(C5[be2][mu2][0]||C5[be2][mu2][1]))
				      {
					int se=epsilon[a1][b1][c1]*epsilon[a2][b2][c2];
					
					complex ter;
					unsafe_complex_prod(ter,S0_SL[im_3pts][rlike][ivol][b1][b2][al1][al2],S0_SL[im_3pts][rlike][ivol][c1][c2][be1][be2]);
					complex_summ_the_prod(ter,S0_SL[im_3pts][rlike][ivol][b1][b2][al1][be2],S0_SL[im_3pts][rlike][ivol][c1][c2][be1][al2]);
					
					safe_complex_prod(ter,C5[be1][mu1],ter);
					safe_complex_prod(ter,C5[be2][mu2],ter);
					safe_complex_prod(ter,Proj[2][al2][al1],ter); //create z+ polarized proton
					
					for(int ri=0;ri<2;ri++)
					  if(se==1) temp[a2][a1][mu2][mu1][ri]+=ter[ri];
					  else      temp[a2][a1][mu2][mu1][ri]-=ter[ri];
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
  
  for(int ivol=0;ivol<loc_vol;ivol++)
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
			if(Proj[2][ga2][ga1][0]||Proj[2][ga2][ga1][1])
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
							complex_prod_with_real(part,part,epsilon[x][b2][c2]);
							safe_complex_prod(part,S0_SL[im_3pts][rdislike][ivol][b1][b2][be1][be2],part);
							safe_complex_prod(part,Proj[2][ga2][ga1],part); //spin z+ polarized proton
							
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
							  complex_prod_with_real(part,part,epsilon[a2][b2][x]);
							  safe_complex_prod(part,S0_SL[im_3pts][rdislike][ivol][b1][b2][be1][be2],part);
							  safe_complex_prod(part,Proj[2][ga2][ga1],part);
							  
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

void calculate_S1_like(int rlike,int rdislike)
{
  for(int id_sink=0;id_sink<4;id_sink++)
    for(int ic_sink=0;ic_sink<3;ic_sink++)
      { //take the source and put g5 on the source, and pre-rotate
	for(int ivol=0;ivol<loc_vol;ivol++)
	  {
	    for(int id_sour=0;id_sour<4;id_sour++)
	      for(int ic_sour=0;ic_sour<3;ic_sour++)
		for(int ri=0;ri<2;ri++)
		  if(id_sour<2)source[ivol][id_sour][ic_sour][ri]= seq_source[ivol][ic_sink][ic_sour][id_sink][id_sour][ri];
		  else         source[ivol][id_sour][ic_sour][ri]=-seq_source[ivol][ic_sink][ic_sour][id_sink][id_sour][ri];
	    if(rdislike==0) safe_spincolor_prod_dirac(source[ivol],source[ivol],&Pminus); //up
	    else            safe_spincolor_prod_dirac(source[ivol],source[ivol],&Pplus);  //dw
	  }

	if(rank==0) printf("\n(S1) like, rlike=%d rdislike=%d, sink index: id=%d, ic=%d\n",rlike,rdislike,id_sink,ic_sink);
	
	tinv-=take_time();
	inv_Q2_cg_left(solDD[0],source,NULL,conf,kappa,mass[im_3pts],niter_max,1,stopping_residue);
	tinv+=take_time();
	
	//use sol_reco[0] as temporary storage. We solve for rdislike
	if(rdislike==0)
	  {
	    apply_Q_left(sol_reco[0],solDD[0],conf,kappa,mass[im_3pts]); //cancel d
	    for(int ivol=0;ivol<loc_vol;ivol++) safe_spincolor_prod_dirac(sol_reco[0][ivol],sol_reco[0][ivol],&Pminus);
	  }
	else
	  {
            apply_Q_left(sol_reco[0],solDD[0],conf,kappa,-mass[im_3pts]); //cancel u
	    for(int ivol=0;ivol<loc_vol;ivol++) safe_spincolor_prod_dirac(sol_reco[0][ivol],sol_reco[0][ivol],&Pplus);
	  }
	
	for(int ivol=0;ivol<loc_vol;ivol++)
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
		  S1[ivol][ic_sink][ic_sour][id_sink][id_sour][ri]=temp[id_sour][ic_sour][ri];
	  }
      }
  
  if(rank==0) printf("like sequential inversions finished\n");
}

void calculate_S1_dislike(int rlike,int rdislike)
{
  for(int id_sink=0;id_sink<4;id_sink++)
    for(int ic_sink=0;ic_sink<3;ic_sink++)
      { //take the source and put g5 on the source
	for(int ivol=0;ivol<loc_vol;ivol++)
	  {
	    for(int id_sour=0;id_sour<4;id_sour++)
	      for(int ic_sour=0;ic_sour<3;ic_sour++)
		for(int ri=0;ri<2;ri++)
		  if(id_sour<2)source[ivol][id_sour][ic_sour][ri]= seq_source[ivol][ic_sink][ic_sour][id_sink][id_sour][ri];
		  else         source[ivol][id_sour][ic_sour][ri]=-seq_source[ivol][ic_sink][ic_sour][id_sink][id_sour][ri];
	    
	    if(rlike==0) safe_spincolor_prod_dirac(source[ivol],source[ivol],&Pminus); //up
	    else         safe_spincolor_prod_dirac(source[ivol],source[ivol],&Pplus);  //dw
	  }

	if(rank==0) printf("\n(S1) dislike, rlike=%d rdislike=%d, sink index: id=%d, ic=%d\n",rlike,rdislike,id_sink,ic_sink);
	
	tinv-=take_time();
	inv_Q2_cg_left(solDD[0],source,NULL,conf,kappa,mass[im_3pts],niter_max,1,stopping_residue);
	tinv+=take_time();
	
	//use sol_reco[0] as temporary storage. We solve for rlike
	if(rlike==0)
	  {
	    apply_Q_left(sol_reco[0],solDD[0],conf,kappa,mass[im_3pts]); //cancel d
	    for(int ivol=0;ivol<loc_vol;ivol++) safe_spincolor_prod_dirac(sol_reco[0][ivol],sol_reco[0][ivol],&Pminus); //up
	  }
	else
	  {
	    apply_Q_left(sol_reco[0],solDD[0],conf,kappa,-mass[im_3pts]); //cancel u
	    for(int ivol=0;ivol<loc_vol;ivol++) safe_spincolor_prod_dirac(sol_reco[0][ivol],sol_reco[0][ivol],&Pplus);  //dw
	  }
	
	for(int ivol=0;ivol<loc_vol;ivol++)
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
		  S1[ivol][ic_sink][ic_sour][id_sink][id_sour][ri]=temp[id_sour][ic_sour][ri];
	  }
      }
  
  if(rank==0) printf("dislike sequential inversions finished\n");
}

//this is needed to check 2pts
void contract_with_source(complex *glb_contr,su3spinspin *eta,su3spinspin *S)
{
  complex *loc_contr=malloc(sizeof(complex)*glb_size[0]);
  memset(loc_contr,0,sizeof(complex)*glb_size[0]);

  for(int ivol=0;ivol<loc_vol;ivol++)
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
  
  complex *contr_2pts=malloc(sizeof(complex)*glb_size[0]);
  
  contract_with_source(contr_2pts,original_source,S1);
  
  if(rank==0)
    {
      fprintf(fout," %+016.16g\t%+016.16g\n",contr_2pts[source_pos[0]][0],contr_2pts[source_pos[0]][1]);
      fclose(fout);
    }
  
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
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    {
      int coor=(glb_coord_of_loclx[loc_site][dir]-source_pos[dir]+glb_size[dir])%glb_size[dir];

      if(coor> glb_size[dir]/2) coor-=glb_size[dir]; //minor distance
      if(coor==glb_size[dir]/2) coor=0; //take away the border (Simpson)
      
      for(int icso=0;icso<3;icso++)
	for(int icsi=0;icsi<3;icsi++)
	  for(int idso=0;idso<4;idso++)
	    for(int idsi=0;idsi<4;idsi++)
	      for(int ri=0;ri<2;ri++)
		S_out[loc_site][icsi][icso][idsi][idso][ri]=S_in[loc_site][icsi][icso][idsi][idso][ri]*coor;
    }
}

//calculate all the 3pts contractions
void calculate_all_3pts_with_current_sequential(int rlike,int rdislike,int rS0,char *path,const char *tag)
{
  int ncontr;
  FILE *fout=open_text_file_for_output(path);
  su3spinspin *supp_S=allocate_su3spinspin(loc_vol,"suppS");;
  
  tcontr-=take_time();
  
  complex *loc_contr=(complex*)malloc(sizeof(complex)*glb_size[0]);
  complex *glb_contr=(complex*)malloc(sizeof(complex)*glb_size[0]);

  //this loop is made in order to avoid duplicating the routine
  for(int norm_chro_EDM=0;norm_chro_EDM<3;norm_chro_EDM++)
    {

      //switch the contraction
      switch(norm_chro_EDM)
	{
	case 0:
	  ncontr=nproton_3pt_contr;
	  break;
	case 1:
	  ncontr=nproton_3pt_chromo_contr;
	  unsafe_apply_chromo_operator_to_su3spinspin(supp_S,Pmunu,S0_SL[im_3pts][rS0]);
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
	  if(norm_chro_EDM==2) apply_dipole_operator(supp_S,S0_SL[im_3pts][rS0],icontr+1);
	  
	  //loop over local node volume
	  for(int loc_site=0;loc_site<loc_vol;loc_site++)
	    {
	      complex point_contr;
	      int glb_t=glb_coord_of_loclx[loc_site][0];
	      
	      //contract the single point
	      switch(norm_chro_EDM)
		{
		case 0:
		  point_proton_sequential_contraction(point_contr,S0_SL[im_3pts][rS0][loc_site],base_gamma[list_3pt_op[icontr]],S1[loc_site]);
		  break;
		case 1:
		  point_proton_sequential_contraction(point_contr,supp_S[loc_site],base_gamma[list_3pt_chromo_op[icontr]],S1[loc_site]);
		  break;
		case 2:
		  point_proton_sequential_contraction(point_contr,supp_S[loc_site],base_gamma[4],S1[loc_site]);
		  break;
		}
	      
	      complex_summ(loc_contr[glb_t],loc_contr[glb_t],point_contr);
	    }
	  
	  //final reduction
	  MPI_Reduce(loc_contr,glb_contr,2*glb_size[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	  
	  if(rank==0)
	    {
	      fprintf(fout," # Three point for rlike=%d, rdislike=%d, %s source\n",rlike,rdislike,tag);

	      switch(norm_chro_EDM)
		{
		case 0:
		  fprintf(fout," # Proton-%s-Proton\n",gtag[list_3pt_op[icontr]]);
		  break;
		case 1:
		  fprintf(fout," # Proton-CHROMO_%s-Proton\n",gtag[list_3pt_op[icontr]]);
		  break;
		case 2:
		  fprintf(fout," # Proton-EDM_%d-Proton\n",icontr+1);
		  break;
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

  tcontr+=take_time();
  
  if(rank==0)
    {
      fclose(fout);
      printf("contractions ultimated\n");
    }

  free(loc_contr);
  free(glb_contr);
  free(supp_S);
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();

  tot_time-=take_time();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  initialize_nucleons(arg[1]);
  
  ///////////////////////////////////////////
  
  prepare_source();
  
  for(int iconf=0;iconf<nconf;iconf++)
    {
      read_conf_and_put_antiperiodic(conf,conf_path[iconf],source_pos[0]);
      
      calculate_S0();
      calculate_all_2pts(out_path[iconf],S0_SS);

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

  tot_time+=take_time();

  ///////////////////////////////////////////

  if(rank==0)
    {
      printf("Total time: %g s\n",tot_time);
      printf("-inversion time: %g%s avg: %d s\n",tinv/tot_time*100,"%",(int)(tinv/nconf/12));
      printf("-contraction time: %g%s\n",tcontr/tot_time*100,"%");
    }

  close_appretto();

  return 0;
}
