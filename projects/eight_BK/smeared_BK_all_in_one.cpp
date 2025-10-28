// This program create the two sources, invert them and perform all the needed contractions

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Three-point disconnected part (Mezzotos)
//               _                                          _       
//   Sum_x{ Tr [ S(md;x,tL) gsource S(ms,x,tL) OP_L  ] Tr [ S(md;x,tR) gsource S(ms,x,tR) OP_R  ]  }      
//
// Three-point connected part (Otto)
//
//   Sum_x{ Tr [ S(md;x,tL) gsource S(ms,x,tL) OP_L  S(md;x,tR) gsource S(ms,x,tR) OP_R  ]  }    
//        _
// where  S is the revert propagator
//   _     +
//   S=g5 S g5
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This is the reference scheme:

/*                             
                S(ms)                           S(ms) 
                .....          			.....
              ..     ..        		      ..     ..
             .         .         (x)	     .	       .	
   gsource  X           X  OP_L        OP_R  X           X gsource
     (tL)    .         .               (tR)  .         .
              ..     ..                       ..      ..
                .....          		  	.....
                S(md)                            S(md)  



                S(ms)        S(ms)
                .....       ,.....
              ..     ..    ,.     ..   
             .         .  ,         .      
   gsource  X           X            X gsource
     (tL)    .         , .          .   (tR)
              ..     .,    ..     ..
                ....,        .....
                S(md)         S(md)

                                    
source |------>---->----->---->| sink

*/

#include <stdlib.h>

#include "nissa.hpp"

using namespace nissa;

//use or not use static limit
int include_static;
double hyp_alpha0,hyp_alpha1,hyp_alpha2;

//gauge info
int igauge_conf=0,ngauge_conf,nanalized_conf;
char conf_path[1024],outfolder[1024];
quad_su3 *conf,*sme_conf,*hyp_conf;
double kappa;

//masses and theta
int ndyn_mass,nmass;
double *mass;
double put_theta[4],old_theta[4]={0,0,0,0};

//source data
int seed,noise_type;
int nsepa,*tsepa;
int nwall,*twall;
char **wall_name;
spincolor *source;
colorspinspin *original_source;

//smearing parameters
double gaussian_kappa,ape_alpha;
int **so_gnit,*so_gnlv; //one per wall
int *si_gnit,si_gnlv;
int ape_niter;

//vectors for the spinor data
int nprop;
colorspinspin **S;
spincolor **cgm_solution,*temp_vec[2];

//cgm inverter parameters
double *stopping_residues;
int niter_max=1000000;

//ottos contractions
complex *contr_otto,*contr_mezzotto;

//two points contractions
int ncontr_2pts;
complex *contr_2pts,*loc_2pts;
int *op1_2pts,*op2_2pts;

//timings
int wall_time,ntot_inv=0;
int ntot_contr_2pts=0;
int ntot_contr_3pts=0;
double tot_prog_time=0,tot_inv_time=0;
double tot_contr_2pts_time=0;
double tot_contr_3pts_time=0;

//number of spectator masses
int nspec;

//return the index of the prop
int iS(int iwall,int sm_lv,int imass,int r)
{
  int i=0;
  for(int iwall_pass=0;iwall_pass<iwall;iwall_pass++) i+=so_gnlv[iwall_pass];
  return r+2*(imass+nmass*(sm_lv+i));
}

//generate the source 
void generate_source(int iwall)
{
  enum rnd_t type[5]={RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_Z2,RND_Z2,RND_Z4};
  generate_spindiluted_source(original_source,type[noise_type],twall[iwall]);
}

//Parse all the input file
void initialize_Bk(int narg,char **arg)
{

  // 0) base initializations
  
  //Check arguments
  if(narg<2) CRASH("No input file specified!\n");
  //Take init time
  tot_prog_time-=take_time();
  //Open input
  open_input(arg[1]);

  // 1) Read information about the gauge conf
  
  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Read the walltime
  read_str_int("Walltime",&wall_time);
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the source and masses
  
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  //Read the position of additional walls
  read_list_of_ints("NSeparations",&nsepa,&tsepa);
  read_list_of_chars("WallNames",&nwall,&wall_name,10);
  //Allocate twall space
  if(nwall!=nsepa+1) CRASH("nwall=%d != nsepa+1=%d",nwall,nsepa+1);
  twall=nissa_malloc("twal",nwall,int);
  so_gnlv=nissa_malloc("so_gnlv",nwall,int);
  so_gnit=nissa_malloc("so_gnit",nwall,int*);

  //Read the noise type
  read_str_int("NoiseType",&noise_type);
  //Read list of masses
  read_list_of_double_pairs("MassResidues",&nmass,&mass,&stopping_residues);
  ndyn_mass=nmass;
  //Read if to include static limit
  read_str_int("IncludeStatic",&include_static);
  if(include_static)
    {
      read_str_double("HypAlpha0",&hyp_alpha0);
      read_str_double("HypAlpha1",&hyp_alpha1);
      read_str_double("HypAlpha2",&hyp_alpha2);
    }
  
  //If needed add the static point
  if(include_static)
    {
      nmass++;
      
      mass=(double*)realloc((void*)mass,(nmass)*sizeof(double));
      stopping_residues=(double*)realloc((void*)stopping_residues,(nmass)*sizeof(double));
      
      //Add the static point
      mass[ndyn_mass]=1000000;
      stopping_residues[ndyn_mass]=1.e-300;
    }
    
  // 3) Smearing parameters
  
  //Smearing parameters
  read_str_double("ApeAlpha",&ape_alpha);
  read_str_int("ApeNiter",&ape_niter);
  read_str_double("GaussianKappa",&gaussian_kappa);
  for(int iwall=0;iwall<nwall;iwall++)
    {
      char tag_source[100];
      sprintf(tag_source,"Source%01dGaussianNiters",iwall);
      read_list_of_ints(tag_source,&(so_gnlv[iwall]),&(so_gnit[iwall]));
      for(int glv=1;glv<so_gnlv[iwall];glv++)
	if(so_gnit[iwall][glv]<so_gnit[iwall][glv-1])
	  CRASH("gaussian levels have to be sorted in ascending order, but this is not the case for %02 wall!",iwall);
    }
  read_list_of_ints("SinkGaussianNiters",&si_gnlv,&si_gnit);
  for(int glv=1;glv<si_gnlv;glv++)
    if(si_gnit[glv]<si_gnit[glv-1])
      CRASH("gaussian levels have to be sorted in ascending order, but this is not the case for sink!");
  
  // 4) contraction list for eight

  contr_otto=nissa_malloc("contr_otto",16*glb_size[0],complex);
  contr_mezzotto=nissa_malloc("contr_mezzotto",16*glb_size[0],complex);
  
  read_str_int("NSpec",&nspec);
  if(nspec>nmass) CRASH("Nspec>nmass!!!");
  read_str_int("NContrTwoPoints",&ncontr_2pts);
  contr_2pts=nissa_malloc("contr_2pts",ncontr_2pts*glb_size[0],complex);
  loc_2pts=nissa_malloc("loc_2pts",ncontr_2pts*glb_size[0],complex);
  op1_2pts=nissa_malloc("op1",ncontr_2pts,int);
  op2_2pts=nissa_malloc("op2",ncontr_2pts,int);
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1_2pts[icontr]));
      read_int(&(op2_2pts[icontr]));
    }
  
  read_str_int("NGaugeConf",&ngauge_conf);
  
  ////////////////////////////////////// end of input reading/////////////////////////////////

  //Allocate gauge conf
  conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  sme_conf=nissa_malloc("sme_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  if(include_static) hyp_conf=nissa_malloc("hyp_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //Allocate all the propagators colorspinspin vectors
  nprop=0;
  for(int iwall=0;iwall<nwall;iwall++) nprop+=2*so_gnlv[iwall]*nmass;
  MASTER_PRINTF("Number of propagator to be allocated: %d\n",nprop);
  S=nissa_malloc("S",nprop,colorspinspin*);
  for(int iprop=0;iprop<nprop;iprop++) S[iprop]=nissa_malloc("S[i]",loc_vol+bord_vol,colorspinspin);
  
  //Allocate nmass spincolors, for the cgm solutions
  cgm_solution=nissa_malloc("cgm_solution",ndyn_mass,spincolor*);
  for(int imass=0;imass<ndyn_mass;imass++) cgm_solution[imass]=nissa_malloc("cgm_solution",loc_vol+bord_vol,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol,spincolor);

  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  
  //Allocate the colorspinspin of the source
  original_source=nissa_malloc("original_source",loc_vol,colorspinspin);
}

//find a new conf
int read_conf_parameters()
{
  int ok_conf;

  do
    {
      //Gauge path
      read_str(conf_path,1024);
      
      //Source coord
      read_int(twall);
      
      //Out folder
      read_str(outfolder,1024);
      
      //Check if the conf exist
      MASTER_PRINTF("Considering configuration \"%s\" with output path \"%s\".\n",outfolder,conf_path);
      ok_conf=!(dir_exists(outfolder));
      if(ok_conf)
        {
          int ris=create_dir(outfolder);
          if(ris==0) MASTER_PRINTF(" Output path \"%s\" not present: configuration \"%s\" not yet analyzed, starting.\n",outfolder,conf_path);
          else
            {
              ok_conf=0;
              MASTER_PRINTF(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
            }
        }
      else
        MASTER_PRINTF(" Output path \"%s\" already present: configuration \"%s\" already analyzed, skipping.\n",outfolder,conf_path);
      igauge_conf++;
    }
  while(!ok_conf && igauge_conf<ngauge_conf);
  
  MASTER_PRINTF("\n");
  
  return ok_conf;
}

//check if there is enough residual time for another conf
int check_residual_time()
{
  int enough_time;
  
  if(igauge_conf<ngauge_conf)
    {
      double spent_time=take_time()+tot_prog_time;
      double remaining_time=wall_time-spent_time;
      double ave_time=(nanalized_conf>0) ? (spent_time/nanalized_conf) : 0;
      double pess_time=ave_time*1.1;
      
      enough_time=remaining_time>pess_time;
      
      MASTER_PRINTF("\n");
      MASTER_PRINTF("-average running time: %lg secs per conf,\n",ave_time);
      MASTER_PRINTF("-pessimistical estimate: %lg secs per conf\n",pess_time);
      MASTER_PRINTF("-remaining time: %lg secs\n",remaining_time);
      if(!enough_time) MASTER_PRINTF("Not enough time for another conf, so exiting.\n");
    }
  else
    {
      MASTER_PRINTF("Finished all the confs, so exiting.\n");
      enough_time=0;
    }

  return enough_time;
}

//load the conf, smear it and put boundary cond
void load_gauge_conf()
{
  //load the gauge conf
  double time=-take_time();
  read_ildg_gauge_conf(conf,conf_path);
  time+=take_time();
  MASTER_PRINTF("\nTime needed to load conf %s: %g s.\n\n",conf_path,time);
  
  //compute plaquette
  MASTER_PRINTF("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  
  //prepare the smerded version
  ape_spatial_smear_conf(sme_conf,conf,ape_alpha,ape_niter);
  
  //calculate smerded plaquette
  MASTER_PRINTF("smerded plaq: %.18g\n",global_plaquette_lx_conf(sme_conf));
  
  //Put the anti-periodic condition on the temporal border
  old_theta[0]=0;
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,0);
  
  //prepare the hypped conf if needed
  if(include_static)
    {
      MASTER_PRINTF("Hypping the conf\n");
      hyp_smear_conf_dir(hyp_conf,conf,hyp_alpha0,hyp_alpha1,hyp_alpha2,0);
      MASTER_PRINTF("hypped plaq: %.18g\n",global_plaquette_lx_conf(hyp_conf));
    }
}

//calculate the propagators
void calculate_S(int iwall)
{
  MASTER_PRINTF("\n");
  
  for(int id=0;id<4;id++)
    {
      //loop over the source dirac index
      get_spincolor_from_colorspinspin(source,original_source,id);
      safe_dirac_prod_spincolor(source,base_gamma+5,source);
      
      //loop over smerding levels of the source
      for(int so_glv=0;so_glv<so_gnlv[iwall];so_glv++)
	{
	  MASTER_PRINTF("\n");
	  
	  int so_gnit_to_app=((so_glv==0) ? so_gnit[iwall][so_glv] : (so_gnit[iwall][so_glv]-so_gnit[iwall][so_glv-1]));
	  gaussian_smearing(source,source,sme_conf,gaussian_kappa,so_gnit_to_app);
	  
	  double part_time=-take_time();
	  MASTER_PRINTF("\n");
	  inv_tmQ2_cgm(cgm_solution,conf,kappa,mass,ndyn_mass,niter_max,stopping_residues,source);
	  part_time+=take_time();ntot_inv++;tot_inv_time+=part_time;
	  MASTER_PRINTF("\nFinished the wall %d inversion, dirac index %d, sm lev %d in %g sec\n\n",
			     iwall,id,so_glv,part_time);
	  
	  for(int imass=0;imass<ndyn_mass;imass++)
	    { //reconstruct the doublet
	      reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass[imass],cgm_solution[imass]);
	      
	      MASTER_PRINTF("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	      for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
		{
		  int iprop=iS(iwall,so_glv,imass,r);
		  NISSA_LOC_VOL_LOOP(ivol) put_spincolor_into_colorspinspin(S[iprop][ivol],temp_vec[r][ivol],id);
		}
	    }
	}
    }
  
  //rotate dynamical quarks to physical basis
  MASTER_PRINTF("\nRotating propagators\n");
  for(int so_glv=0;so_glv<so_gnlv[iwall];so_glv++)
    for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
      for(int imass=0;imass<ndyn_mass;imass++) //put the (1+ig5)/sqrt(2) factor
	{
	  int iprop=iS(iwall,so_glv,imass,r);
	  rotate_vol_colorspinspin_to_physical_basis(S[iprop],!r,!r);
	}
  
  //compute static prop, if needed
  if(include_static)
    {
      color *stat_source=nissa_malloc("StatSource",loc_vol+bord_vol,color);
      
      MASTER_PRINTF("\nComputing static prop\n");
      get_color_from_colorspinspin(stat_source,original_source,0,0);

      for(int so_glv=0;so_glv<so_gnlv[iwall];so_glv++)
	{
	  MASTER_PRINTF("\n");
	  
	  int so_gnit_to_app=((so_glv==0) ? so_gnit[iwall][so_glv] : (so_gnit[iwall][so_glv]-so_gnit[iwall][so_glv-1]));
	  gaussian_smearing(stat_source,stat_source,sme_conf,gaussian_kappa,so_gnit_to_app);
	  
	  //compute for r=0 and copy in r=1
	  compute_Wstat_stoch_prop(S[iS(iwall,so_glv,ndyn_mass,0)],hyp_conf,0,twall[iwall],stat_source);
	  vector_copy(S[iS(iwall,so_glv,ndyn_mass,1)],S[iS(iwall,so_glv,ndyn_mass,0)]);
	}
      
      nissa_free(stat_source);
    }
}

//Compute the connected part of bag correlator
void Bk_eights(colorspinspin *SL1,colorspinspin *SL2,colorspinspin *SR1,colorspinspin *SR2)
{
  //Temporary vectors for the internal gamma
  dirac_matr tsink[16];//tsource[16]
  for(int igamma=0;igamma<16;igamma++)
    {
      //Put the two gamma5 needed for the revert of the d spinor
      //tsource[igamma]=base_gamma[0]; //g5*g5, so commented and never used
      dirac_prod(&(tsink[igamma]),&(base_gamma[5]),&(base_gamma[igamma]));
    }
  
  //Call the routine which does the real contraction for the Mezzotto and the Otto
  trace_id_css_dag_g_css_id_css_dag_g_css(contr_otto,SL1,tsink,SL2,SR1,tsink,SR2,16);
  trace_id_css_dag_g_css_times_trace_id_css_dag_g_css(contr_mezzotto,SL1,tsink,SL2,SR1,tsink,SR2,16);
}

//Compute mesonic two points
void meson_two_points(colorspinspin *s1,colorspinspin *s2)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr_2pts],t2[ncontr_2pts];
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Put the two gamma5 needed for the revert of the first spinor
      dirac_prod(&(t1[icontr]),&(base_gamma[op1_2pts[icontr]]),&(base_gamma[5]));
      dirac_prod(&(t2[icontr]),&(base_gamma[5]),&(base_gamma[op2_2pts[icontr]]));
    }
  //Call for the routine which perform actually real contraction
  trace_g_css_dag_g_css(contr_2pts,loc_2pts,t1,s1,t2,s2,ncontr_2pts);
}

//print all the passed contractions to the file
void print_ottos_contractions_to_file(FILE *fout)
{
  double norm=1.0;
  
  if(rank==0)
    {
      for(int icontr=0;icontr<16;icontr++)
	{
	  fprintf(fout,"\n");
	  print_contraction_to_file(fout,icontr,5,contr_mezzotto+icontr*glb_size[0],twall[0],"DISCONNECTED ",norm);
	  fprintf(fout,"\n");
	  print_contraction_to_file(fout,icontr,5,contr_otto+icontr*glb_size[0],twall[0],"CONNECTED ",norm);
	}
      fprintf(fout,"\n");
    }
}

//print all the passed contractions to the file
void print_two_points_contractions_to_file(FILE *fout)
{
  double norm=1.0;
  
  if(rank==0)
    {
      for(int icontr=0;icontr<ncontr_2pts;icontr++)
	{
	  fprintf(fout,"\n");
	  print_contraction_to_file(fout,op2_2pts[icontr],op1_2pts[icontr],contr_2pts+icontr*glb_size[0],twall[0],"",norm);
	}
      fprintf(fout,"\n");
    }
}

//Calculate and print to file all the contractions
void calculate_all_contractions()
{
  MASTER_PRINTF("Computing all contractions\n");
  tot_contr_3pts_time-=take_time();
  
  //loop over left-right wall combo
  for(int iwL=0;iwL<nwall;iwL++)
    for(int iwR=iwL+1;iwR<nwall;iwR++)
      //loop over smearing of the left and right wall
      for(int sm_lv_L=0;sm_lv_L<so_gnlv[iwL];sm_lv_L++)
	for(int sm_lv_R=0;sm_lv_R<so_gnlv[iwR];sm_lv_R++)
	  {
	    char path_bag[1024];
	    sprintf(path_bag,"%s/otto_w%s_%s_%02d_%02d",outfolder,wall_name[iwL],wall_name[iwR],
		    so_gnit[iwL][sm_lv_L],so_gnit[iwR][sm_lv_R]);
	    FILE *fout_bag=open_text_file_for_output(path_bag);
	    
	    int tsepar=(twall[iwR]+glb_size[0]-twall[iwL])%glb_size[0];
	    
	    master_fprintf(fout_bag," # LEFT_WALL_t=%d , RIGHT_WALL_t=%d , tseparat=%d\n\n",twall[iwL],twall[iwR],tsepar);
	      
	    //Ottos
	    for(int im2=0;im2<nmass;im2++)
	      for(int r2=0;r2<2;r2++)
		for(int im1=0;im1<nspec;im1++)
		  for(int r1=0;r1<2;r1++)
		    for(int r3=0;r3<2;r3++)
		      {
			int r4=1-(r1+r2+r3)%2;
			
			master_fprintf(fout_bag," # m1=%f r1=%d , m2=%f r2=%d , m3=%f r3=%d , m4=%f r4=%d",
				       mass[im1],r1,mass[im2],r2,mass[im1],r3,mass[im2],r4);
			master_fprintf(fout_bag," . Tseparation=%d\n",tsepar);
			
			int ip1=iS(iwL,sm_lv_L,im1,r1),ip2=iS(iwL,sm_lv_L,im2,r2);
			int ip3=iS(iwR,sm_lv_R,im1,r3),ip4=iS(iwR,sm_lv_R,im2,r4);
			
			Bk_eights(S[ip1],S[ip2],S[ip3],S[ip4]);
			
			ntot_contr_3pts+=32;
			
			print_ottos_contractions_to_file(fout_bag);
		      }
	    
	    close_file(fout_bag);
	  }
  
  tot_contr_3pts_time+=take_time();
  
  ///////////////two points///////////////
  
  tot_contr_2pts_time-=take_time();
  
  //loop over the wall
  for(int iwall=0;iwall<nwall;iwall++)
    //loop over smearing levels of the sink
    for(int si_glv=0;si_glv<si_gnlv;si_glv++)
      {
	//smear all the propagators of iwall on the sink
	int si_gnit_to_app=((si_glv==0) ? si_gnit[si_glv] : (si_gnit[si_glv]-si_gnit[si_glv-1])); 
	if(si_gnit_to_app!=0)
	  for(int so_glv=0;so_glv<so_gnlv[iwall];so_glv++)
	    for(int im=0;im<nmass;im++)
	      for(int r=0;r<2;r++)
		{
		  int iprop=iS(iwall,so_glv,im,r);
		  gaussian_smearing(S[iprop],S[iprop],sme_conf,gaussian_kappa,si_gnit_to_app);
		}
	
	//loop over all source smearing level
	for(int so_glv=0;so_glv<so_gnlv[iwall];so_glv++)
	  {
	    char path_2pts[1024];
	    sprintf(path_2pts,"%s/two_points_w%s_%02d_%02d",outfolder,wall_name[iwall],so_gnit[iwall][so_glv],si_gnit[si_glv]);
	    FILE *fout_2pts=open_text_file_for_output(path_2pts);
	    MASTER_PRINTF("opening file %s\n",path_2pts);
	    //loop over all the combos
	    for(int im2=0;im2<nmass;im2++)
	      for(int r2=0;r2<2;r2++)
		for(int im1=0;im1<nspec;im1++)
		  for(int r1=0;r1<2;r1++)
		    {
		      master_fprintf(fout_2pts," # m1=%f r1=%d , m2=%f r2=%d , wall=%d , ",mass[im1],r1,mass[im2],r2,iwall);
		      master_fprintf(fout_2pts," sme_source=%d sme_sink=%d\n",so_gnit[iwall][so_glv],si_gnit[si_glv]);
		      
		      int iprop1=iS(iwall,so_glv,im1,r1);
		      int iprop2=iS(iwall,so_glv,im2,r2);
		      
		      meson_two_points(S[iprop1],S[iprop2]);
		      
		      print_two_points_contractions_to_file(fout_2pts);
		      VERBOSITY_LV3_MASTER_PRINTF("printing contractions between %d and %d\n",iprop1,iprop2);
		      ntot_contr_2pts+=ncontr_2pts;
		    }
	    close_file(fout_2pts);
	  }
      }
  
  tot_contr_2pts_time+=take_time();
}

//analize a single configuration
void analize_conf()
{
  //Determine the position of all the wall starting from the distance
  for(int iwall=1;iwall<nwall;iwall++) twall[iwall]=(twall[0]+tsepa[iwall-1])%glb_size[0];
  
  //Invert propagators
  MASTER_PRINTF("Going to invert: %d walls\n",nwall);
  for(int iwall=0;iwall<nwall;iwall++)
    {
      generate_source(iwall);
      calculate_S(iwall);
    }
  
  //Perform all the contractions
  calculate_all_contractions();
  
  nanalized_conf++;
}

//Finalization
void close_Bk()
{
  //take final time
  tot_prog_time+=take_time();

  MASTER_PRINTF("\n");
  MASTER_PRINTF("Total time: %g secs to analize %d configurations (%f secs avg), of which:\n",
		tot_prog_time,nanalized_conf,tot_prog_time/nanalized_conf);
  MASTER_PRINTF(" - %02.2f%s to perform %d inversions (%f secs avg)\n",
		tot_inv_time/tot_prog_time*100,"%",ntot_inv,tot_inv_time/ntot_inv);
  MASTER_PRINTF(" - %02.2f%s to perform %d 3pts contr. (%f secs avg)\n",
		tot_contr_3pts_time/tot_prog_time*100,"%",ntot_contr_3pts,tot_contr_3pts_time/ntot_contr_3pts);
  MASTER_PRINTF(" - %02.2f%s to perform %d 2pts contr. (%f secs avg)\n",
		tot_contr_2pts_time/tot_prog_time*100,"%",ntot_contr_2pts,tot_contr_2pts_time/ntot_contr_2pts);
  
  nissa_free(twall);
  nissa_free(op1_2pts);nissa_free(op2_2pts);
  nissa_free(contr_otto);nissa_free(contr_mezzotto);
  nissa_free(conf);nissa_free(sme_conf);nissa_free(contr_2pts);nissa_free(loc_2pts);
  if(include_static) nissa_free(hyp_conf);
  for(int imass=0;imass<ndyn_mass;imass++) nissa_free(cgm_solution[imass]);
  nissa_free(cgm_solution);
  nissa_free(source);nissa_free(original_source);
  for(int iprop=0;iprop<nprop;iprop++) nissa_free(S[iprop]);
  nissa_free(S);
  nissa_free(temp_vec[1]);nissa_free(temp_vec[0]);
  nissa_free(so_gnlv);
  nissa_free(so_gnit);
}

void in_main(int narg,char **arg)
{
  //Inner initialization
  initialize_Bk(narg,arg);

  //Loop over configurations

  //Find if there is another conf to analize and time to analize it
  while(check_residual_time() && read_conf_parameters())
    {
      //Load the gauge conf
      load_gauge_conf();      
      
      //Analize it
      analize_conf();
    }
  
  //Finalization
  close_input();
  close_Bk(); 
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();

  return 0;
}
