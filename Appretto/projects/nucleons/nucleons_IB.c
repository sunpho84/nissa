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
int nmass;
double *mass;
double kappa;
as2t_su3 *Pmunu;

//source
int source_pos[4];
spincolor *source,*temp_source;
su3spinspin *original_source;

//the propagators
int nprops;
int get_iprop(int imass,int r,int sme_sink,int SING_DOUB){return ((imass*2+r)*2+sme_sink)*2+SING_DOUB;}
su3spinspin **S;

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
int tsink;

spinspin Proj[3]; //projectors over N and N*, and 00 compont of N (in the spinorial representation)
int Proj_ind1[8]={0,0,1,1,2,2,3,3},Proj_ind2[8]={0,2,1,3,0,2,1,3}; //indices different from 0
int Proj_couples[4][2]={{0,5},{2,7},{1,4},{3,6}}; //independent couples
int Proj_entr[3][4]={{1,1,-1,-1},{1,1,1,1},{1,0,-1,0}};
spinspin C5; //C*gamma5

//two points contractions
int nproton_2pt_contr=1;
int list_2pt_op1[32]={5};
int list_2pt_op2[32]={0};
int compute_also_SL_2pts;

//                      e_00x   e_01x    e_02x     e_10x    e_11x   e_12x     e_20x   e_21x    e_22x
int epsilon[3][3][3]={{{0,0,0},{0,0,1},{0,-1,0}},{{0,0,-1},{0,0,0},{1,0,0}},{{0,1,0},{-1,0,0},{0,0,0}}};
//same thing, but in a compact way
int eps_pos[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
int eps_sig[6]={1,-1,-1,1,1,-1};
//even more compact
int eps1[3]={2,0,1},eps2[3]={1,2,0};

//timings
double tinv=0,tcontr_2pt=0,tot_time=0;

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
  
  // 4) 2pts and insertion info
  //compute also SL 2pts?
  read_str_int("ComputeAlsoSL2pts",&compute_also_SL_2pts);
  
  close_input();
  
  ///////////////////// Allocate the various spinors ///////////////////////
  
  original_source=allocate_su3spinspin(loc_vol,"original_source");
  
  source=allocate_spincolor(loc_vol+loc_bord,"source");
  temp_source=allocate_spincolor(loc_vol+loc_bord,"temp_source");
  
  //allocate solDD
  solDD=(spincolor**)malloc(sizeof(spincolor*)*nmass);
  for(int imass=0;imass<nmass;imass++) solDD[imass]=allocate_spincolor(loc_vol+loc_bord,"solDD");
  
  //allocate S
  nprops=get_iprop(nmass,1,1,1);
  S=(su3spinspin**)malloc(sizeof(su3spinspin***)*nprops);
  for(int iprop=0;iprop<nprops;iprop++) S[iprop]=allocate_su3spinspin(loc_vol+loc_bord,"S");
  
  sol_reco[0]=allocate_spincolor(loc_vol,"solution_reco[0]");
  sol_reco[1]=allocate_spincolor(loc_vol,"solution_reco[1]");
}

//read a configuration and put anti-periodic condition at the slice tsource-1
void read_conf_and_put_antiperiodic(quad_su3 *conf,char *conf_path,int tsource)
{
  read_gauge_conf(conf,conf_path);
  
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

//perform the first inversion to produce the propagator for u and d
void calculate_S()
{
  for(int ic_sour=0;ic_sour<3;ic_sour++)
    for(int id_sour=0;id_sour<4;id_sour++)
      {
	
	// =====================================================================
	
	// 1) prepare the source
	
	if(rank==0) printf("\n(S1) source index: id=%d, ic=%d\n",id_sour,ic_sour);
	
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
		int iprop_SL=get_iprop(imass,r,0,0);
		int iprop_SS=get_iprop(imass,r,1,0);
		
		for(int ivol=0;ivol<loc_vol;ivol++)
		  {
		    //put the anti-periodic condition on the propagator
		    int dt=glb_coord_of_loclx[ivol][0]-source_pos[0];
		    double arg=M_PI*dt/glb_size[0];
		    complex phase={cos(arg),sin(arg)};
		    
		    unsafe_spincolor_prod_complex(temp_source[ivol],sol_reco[r][ivol],phase);
		    put_spincolor_into_su3spinspin(S[iprop_SL][ivol],temp_source[ivol],id_sour,ic_sour);
		  }
		
		//smerd the sink
		jacobi_smearing(source,temp_source,smea_conf,jacobi_kappa,jacobi_niter);
		for(int ivol=0;ivol<loc_vol;ivol++)
		  put_spincolor_into_su3spinspin(S[iprop_SS][ivol],source[ivol],id_sour,ic_sour);
	      }
	  }
      }

  //put the (1+ig5)/sqrt(2) factor
  for(int imass=0;imass<nmass;imass++)
    for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
      {
	int iprop_SL=get_iprop(imass,r,0,0);
	int iprop_SS=get_iprop(imass,r,1,0);
	
	for(int ivol=0;ivol<loc_vol;ivol++)
	  for(int ic1=0;ic1<3;ic1++)
	    for(int ic2=0;ic2<3;ic2++)
	      {
		rotate_spinspin_to_physical_basis(S[iprop_SL][ivol][ic1][ic2],!r,!r);
		rotate_spinspin_to_physical_basis(S[iprop_SS][ivol][ic1][ic2],!r,!r);
	      }
      }

  if(rank==0) printf(" final rotations performed\n");
}

//perform the double integrated inversion
void calculate_SS()
{
  //loop over mass index
  for(int imass=0;imass<nmass;imass++)
    for(int r=0;r<2;r++)
      {
	
	for(int ic_sour=0;ic_sour<3;ic_sour++)
	  for(int id_sour=0;id_sour<4;id_sour++)
	    {
	      
	      // =====================================================================
	      
	      // 1) prepare the source
	      
	      if(rank==0) printf("\n(S2) source index: id=%d, ic=%d\n",id_sour,ic_sour);
	
	      //take the source and put the scalar in twisted basis (so that it ignore g5), remove anti-per.
	      int iprop=get_iprop(imass,r,0,0);
	      for(int ivol=0;ivol<loc_vol;ivol++)
		{
		  int dt=(glb_coord_of_loclx[ivol][0]+source_pos[0])%glb_size[0]-source_pos[0];
		  double arg=M_PI*dt/glb_size[0];
		  complex phase={cos(arg),sin(arg)};
		  spincolor buf;
		  
		  get_spincolor_from_su3spinspin(buf,S[iprop][ivol],id_sour,ic_sour);
		  unsafe_spincolor_prod_complex(source[ivol],buf,phase);
		}
	      
	      //============================================================================
	      
	      // 2) peform the inversion taking time
	      
	      tinv-=take_time();
	      inv_Q2_cg(solDD[imass],source,source,conf,kappa,mass[imass],niter_max,1,stopping_residue);
	      tinv+=take_time();
	      
	      if(rank==0) printf("inversions finished\n");
	      
	      //============================================================================
	      
	      // 3) reconstruct the doublet and smerd the sink
	      
	      //reconstruct the doublet
	      if(r==0) apply_Q(sol_reco[0],solDD[imass],conf,kappa,+mass[imass]);
	      else     apply_Q(sol_reco[0],solDD[imass],conf,kappa,-mass[imass]);
	      
	      //convert the id-th spincolor into the colorspinspin and prepare the sink smerded version
	      int iprop_SL=get_iprop(imass,r,0,1);
	      int iprop_SS=get_iprop(imass,r,1,1);
	      
	      for(int ivol=0;ivol<loc_vol;ivol++)
		{
		  //put the anti-periodic condition on the propagator
		  int dt=glb_coord_of_loclx[ivol][0]-source_pos[0];
		  double arg=M_PI*dt/glb_size[0];
		  complex phase={cos(arg),sin(arg)};
		  spincolor buf;
		  
		  unsafe_spincolor_prod_complex(buf,sol_reco[r][ivol],phase);
		  
		  //pass to phys basis
		  if(r==0) safe_dirac_prod_spincolor(buf,&Pminus,buf);//up
		  else safe_dirac_prod_spincolor(buf,&Pplus,buf); //dw
		  
		  put_spincolor_into_su3spinspin(S[iprop_SL][ivol],buf,id_sour,ic_sour);
		}
	      
	      //smerd the sink
	      jacobi_smearing(source,temp_source,smea_conf,jacobi_kappa,jacobi_niter);
	      for(int ivol=0;ivol<loc_vol;ivol++)
		put_spincolor_into_su3spinspin(S[iprop_SS][ivol],source[ivol],id_sour,ic_sour);
	    }
      }
}

//Calculate the proton contraction for a single point
void local_diquark(diquark *diq,su3spinspin *Sa,su3spinspin *Sb)
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
		  unsafe_complex_prod  (diq[l][b][b1][al][ga][al1][ga1],Sa[l][a1][a][al1][al],Sb[l][c1][c][ga1][ga]);
		  complex_subt_the_prod(diq[l][b][b1][al][ga][al1][ga1],Sa[l][a1][c][al1][ga],Sb[l][c1][a][ga1][al]);
		  //both epsilon index (at fixed b) exchanged, so summ
		  complex_summ_the_prod(diq[l][b][b1][al][ga][al1][ga1],Sa[l][c1][c][al1][al],Sb[l][a1][a][ga1][ga]);
		  complex_subt_the_prod(diq[l][b][b1][al][ga][al1][ga1],Sa[l][c1][a][al1][ga],Sb[l][a1][c][ga1][al]);
		  //now only b indices (a and c) exchanged, so subt
		  complex_subt_the_prod(diq[l][b][b1][al][ga][al1][ga1],Sa[l][a1][c][al1][al],Sb[l][c1][a][ga1][ga]);
		  complex_summ_the_prod(diq[l][b][b1][al][ga][al1][ga1],Sa[l][a1][a][al1][ga],Sb[l][c1][c][ga1][al]);
		  //again only b1 indices (a1 and c1) exchanged, so subt
		  complex_subt_the_prod(diq[l][b][b1][al][ga][al1][ga1],Sa[l][c1][a][al1][al],Sb[l][a1][c][ga1][ga]);
		  complex_summ_the_prod(diq[l][b][b1][al][ga][al1][ga1],Sa[l][c1][c][al1][ga],Sb[l][a1][a][ga1][al]);
		}
}

//close the diquark with the third propagator
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
void calculate_all_2pts(char *path,int sme_sink,int doub1,int doub2,int doub3)
{
  //output file
  FILE *output=open_text_file_for_output(path);

  //take initial time
  tcontr_2pt-=take_time();
  
  //tag for the contraction file
  char pm_tag[2][2]={"+","-"};
  
  //diquark term and all-open dirac indices term
  diquark *diq=(diquark*)malloc(sizeof(diquark)*loc_vol);
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
	int iprop_like1=get_iprop(im_like,rlike,sme_sink,doub1);
	int iprop_like2=get_iprop(im_like,rlike,sme_sink,doub2);
	local_diquark(diq,S[iprop_like1],S[iprop_like2]);
	
	//now close with the third propagator leaving all dirac indices open, and makes global reduction
	for(int im_dislike=0;im_dislike<nmass;im_dislike++)
	  for(int rdislike=0;rdislike<2;rdislike++)
	    {
	      int iprop_dislike=get_iprop(im_dislike,rdislike,sme_sink,doub3);
	      if(rank==0) fprintf(output,"  # Two points for m_like=%g rlike=%d, m_dislike=%g rdislike=%d\n",
				  mass[im_like],rlike,mass[im_dislike],rdislike);
	      
	      close_diquark(prot6,diq,S[iprop_dislike]);
	      
	      //now we still have a structure with 6 dirac indices
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
			    complex_summ_the_prod(prot4[t][al][be][ga][ga1],o3.entr[al1],prot6[t][al][be][ga][al1][be1][ga1]);
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
		  
		  if(rank==0) //header
		    {
		      fprintf(output,"   # %s%s-P5S0\n",gtag[list_2pt_op1[icontr]],gtag[list_2pt_op2[icontr]]);
		      
		      for(int t=0;t<glb_size[0];t++)
			for(int ind_p=0;ind_p<4;ind_p++)
			  for(int ico=0;ico<2;ico++)
			    {
			      int ip=Proj_couples[ind_p][ico];
			      
			      int ga=Proj_ind1[ip];
			      int ga1=Proj_ind2[ip];
			      
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
			  if(rank==0)
			    {
			      if(nns<2) fprintf(output,"# Contraction with (1%sg4)/2\n",pm_tag[nns]);
			      else      fprintf(output,"# Contraction with (1+g4)_00(spinorial)/2\n");
			      
			      for(int tt=0;tt<glb_size[0];tt++)
				{
				  int t=(tt+source_pos[0])%glb_size[0];
				  
				  complex contr={0,0};
				  for(int ind_p=0;ind_p<4;ind_p++)
				    {
				      if(Proj_entr[nns][ind_p]==+1) complex_summ(contr,contr,prot[t][ind_p]);
				      if(Proj_entr[nns][ind_p]==-1) complex_subt(contr,contr,prot[t][ind_p]);
				    }
				  fprintf(output," %+016.16g\t%+016.16g\n",contr[0]/2,contr[1]/2);
				}
			      fprintf(output,"\n");
			    }
			}
		    }
		}
	    }
      }
  
  tcontr_2pt+=take_time();
  
  if(rank==0)
    {
      printf("contractions finished\n");
      fclose(output);
    }
  
  free(diq);
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
  
  //loop over configurations
  for(int iconf=0;iconf<nconf;iconf++)
    {
      read_conf_and_put_antiperiodic(conf,conf_path[iconf],source_pos[0]);
      
      calculate_S();
      calculate_SS();

      for(int doub1=0;doub1<2;doub1++)
	for(int doub2=0;doub2<2;doub2++)
	  for(int doub3=0;doub3<2;doub3++)
	    {
	      //prepare the SL output path
	      char path_SL[1024],path_SS[1024];
	      sprintf(path_SL,"%s_SL_%d%d%d",out_path[iconf],doub1,doub2,doub3);
	      sprintf(path_SS,"%s_SS_%d%d%d",out_path[iconf],doub1,doub2,doub3);
	      
	      //now print both SL and SS
	      if(compute_also_SL_2pts) calculate_all_2pts(path_SL,0,doub1,doub2,doub3);
	      calculate_all_2pts(path_SS,1,doub1,doub2,doub3);
	    }
    }

  tot_time+=take_time();

  ///////////////////////////////////////////

  if(rank==0)
    {
      printf("Total time: %g s\n",tot_time);
      printf("-inversion time: %g%s avg: %d s\n",tinv/tot_time*100,"%",(int)(tinv/nconf/12));
      printf("-contraction time for 2pts: %g%s\n",tcontr_2pt/tot_time*100,"%");
    }

  close_appretto();

  return 0;
}
