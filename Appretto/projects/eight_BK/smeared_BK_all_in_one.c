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

#include <mpi.h>
#include <lemon.h>

#include "appretto.h"
#include "addendum.c"

//gauge info
int ngauge_conf;
char conf_path[1024];
quad_su3 *conf,*sme_conf;
double kappa;

int nmass;
double *mass;
double put_theta[4],old_theta[4]={0,0,0,0};

//source data
int seed,noise_type;
int nsepa,*tsepa;
int nwall,*twall;
double rad2=1.414213562373095048801688724209;
spincolor *source;
colorspinspin *original_source;

//smearing parameters
double jacobi_kappa,ape_alpha;
int *jacobi_niter,ape_niter;
int nsm_lv;

//vectors for the spinor data
int nprop;
colorspinspin **S;
spincolor **cgmms_solution,*reco_solution[2];

//cgmms inverter parameters
double stopping_residue;
double minimal_residue;
int stopping_criterion;
int niter_max;

//ottos contractions
complex *contr_otto,*contr_mezzotto;
char basepath_bag[1024];

//two points contractions
int ncontr_2pts;
complex *contr_2pts;
int *op1_2pts,*op2_2pts;
char basepath_2pts[1024];

//timings
int ninv_tot=0,ncontr_tot=0;
double tot_time=0,inv_time=0,contr_time=0;

//number of spectator masses
int nspec;

//return the index of the prop
int iS(int iwall,int sm_lv,int imass,int r) {return r+2*(imass+nmass*(sm_lv+nsm_lv*iwall));}

//generate the source for the dirac index
void generate_source(int iwall)
{ //reset
  memset(original_source,0,sizeof(colorspinspin)*loc_vol);
  
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    if(glb_coord_of_loclx[loc_site][0]==twall[iwall])
      for(int ic=0;ic<3;ic++)
	{ //real part
	  if(noise_type>=2) original_source[loc_site][ic][0][0][0]=pm_one(loc_site)/rad2;
	  else original_source[loc_site][ic][0][0][0]=noise_type;
	  //imaginary part
	  if(noise_type==4) original_source[loc_site][ic][0][0][1]=pm_one(loc_site)/rad2;
	  
	  for(int d=1;d<4;d++) //copy the other three dirac indexes
	    memcpy(original_source[loc_site][ic][d][d],original_source[loc_site][ic][0][0],sizeof(complex));
	}
}

//Parse all the input file
void initialize_Bk(char *input_path)
{
  open_input(input_path);

  // 1) Read information about the gauge conf
  
  //Read the volume
  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  //Init the MPI grid 
  init_grid(); 
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the source and masses
  
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  init_random(seed);
  //Read the position of additional walls
  read_list_of_ints("NSeparations",&nsepa,&tsepa);
  //Allocate twall space
  nwall=nsepa+1;
  twall=(int*)malloc(sizeof(int)*(nwall));
  //Read the noise type
  read_str_int("NoiseType",&noise_type);
  //Read list of masses
  read_list_of_doubles("NMass",&nmass,&mass);
  
  // 3) Smearing parameters
  
  //Smearing parameters
  read_str_double("ApeAlpha",&ape_alpha);
  read_str_int("ApeNiter",&ape_niter);
  read_str_double("JacobiKappa",&jacobi_kappa);
  read_list_of_ints("JacobiNiters",&nsm_lv,&jacobi_niter);
  if(jacobi_niter[0]!=0)
    {
      if(rank==0) fprintf(stderr,"Error, jacobi level of smearing 0 have to be null, obtained: %d!\n",jacobi_niter[0]);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  for(int iter=1;iter<nsm_lv;iter++)
    if(jacobi_niter[iter]<jacobi_niter[iter-1])
      {
        if(rank==0) fprintf(stderr,"Error, jacobi level %d minor than %d (%d, %d)!\n",iter,iter-1,jacobi_niter[iter],jacobi_niter[iter-1]);
        MPI_Abort(MPI_COMM_WORLD,1);
      }
  
  // 4) Info about the inverter

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
  
  // 5) contraction list for eight

  contr_otto=(complex*)malloc(sizeof(complex)*16*glb_size[0]); 
  contr_mezzotto=(complex*)malloc(sizeof(complex)*16*glb_size[0]); 
  
  read_str_int("NSpec",&nspec);
  read_str_int("NContrTwoPoints",&ncontr_2pts);
  contr_2pts=(complex*)malloc(sizeof(complex)*ncontr_2pts*glb_size[0]); 
  op1_2pts=(int*)malloc(sizeof(int)*ncontr_2pts);
  op2_2pts=(int*)malloc(sizeof(int)*ncontr_2pts);
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1_2pts[icontr]));
      read_int(&(op2_2pts[icontr]));

      if(rank==0 && debug) printf(" contr.%d %d %d\n",icontr,op1_2pts[icontr],op2_2pts[icontr]);
    }
  
  read_str_int("NGaugeConf",&ngauge_conf);
  
  ////////////////////////////////////// end of input reading/////////////////////////////////

  //allocate gauge conf
  conf=allocate_quad_su3(loc_vol+loc_bord,"conf");
  sme_conf=allocate_quad_su3(loc_vol+loc_bord,"sme_conf");
  
  //Allocate all the propagators colorspinspin vectors
  nprop=nwall*nsm_lv*nmass*2;
  if(rank==0) printf("nprop: %d\n",nprop);
  S=(colorspinspin**)malloc(sizeof(colorspinspin*)*nprop);
  for(int iprop=0;iprop<nprop;iprop++) S[iprop]=allocate_colorspinspin(loc_vol,"S");
  
  //Allocate nmass spincolors, for the cgmms solutions
  cgmms_solution=(spincolor**)malloc(sizeof(spincolor*)*nmass);
  for(int imass=0;imass<nmass;imass++) cgmms_solution[imass]=allocate_spincolor(loc_vol+loc_bord,"cgmms_solution");
  reco_solution[0]=allocate_spincolor(loc_vol,"reco_solution[0]");
  reco_solution[1]=allocate_spincolor(loc_vol,"reco_solution[1]");
  //Allocate one spincolor for the source
  source=allocate_spincolor(loc_vol+loc_bord,"source");

  //Allocate the colorspinspin of the source
  original_source=allocate_colorspinspin(loc_vol,"original_source");
}

//load the conf, smear it and put boundary cond
void load_gauge_conf()
{
  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_local_gauge_conf(conf,conf_path);
  //prepared the smerded version and  calculate plaquette
  ape_smearing(sme_conf,conf,ape_alpha,ape_niter);

  communicate_gauge_borders(conf);
  communicate_gauge_borders(sme_conf);
  
  double gplaq=global_plaquette(conf);
  if(rank==0) printf("plaq: %.18g\n",gplaq);
  gplaq=global_plaquette(sme_conf);
  if(rank==0) printf("smerded plaq: %.18g\n",gplaq);

  //Put the anti-periodic condition on the temporal border
  old_theta[0]=0;
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,0);
}

//Finalization
void close_Bk()
{
  if(rank==0)
    {
      printf("\n");
      printf("Total time: %g, of which:\n",tot_time);
      printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
      printf(" - %02.2f%s to perform %d contr. (%2.2gs avg)\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
    }
  
  check_free(conf);check_free(sme_conf);check_free(mass);
  for(int iprop=0;iprop<nprop;iprop++) check_free(S[iprop]);
  check_free(reco_solution[1]);check_free(reco_solution[0]);
  check_free(S);
  check_free(contr_otto);check_free(contr_mezzotto);
  close_appretto();
}

//calculate the propagators
void calculate_S(int iwall)
{
  if(rank==0) printf("\n");
  
  for(int id=0;id<4;id++)
    { //loop over the source dirac index
      for(int ivol=0;ivol<loc_vol;ivol++)
	{
	  get_spincolor_from_colorspinspin(source[ivol],original_source[ivol],id);
	  //put the g5
	  for(int id1=2;id1<4;id1++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic][ri]*=-1;
	}
      
      //loop over smerding levels of the source
      for(int sm_lv=0;sm_lv<nsm_lv;sm_lv++)
	{
	  if(rank==0) printf("\n");
	  
	  int sm_to_app=((sm_lv==0) ? jacobi_niter[sm_lv] : (jacobi_niter[sm_lv]-jacobi_niter[sm_lv-1]));
	  if(rank==0 && debug) printf("Source ");
	  jacobi_smearing(source,source,sme_conf,jacobi_kappa,sm_to_app);
	  
	  double part_time=-take_time();
	  communicate_lx_spincolor_borders(source);
	  if(rank==0) printf("\n");
	  inv_Q2_cgmms(cgmms_solution,source,NULL,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
	  part_time+=take_time();ninv_tot++;inv_time+=part_time;
	  if(rank==0) printf("Finished the wall %d inversion, dirac index %d, sm lev %d in %g sec\n",iwall,id,sm_lv,part_time);
	  
	  for(int imass=0;imass<nmass;imass++)
	    { //reconstruct the doublet
	      reconstruct_doublet(reco_solution[0],reco_solution[1],cgmms_solution[imass],conf,kappa,mass[imass]);
	      if(rank==0) printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	      for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
		{
		  int iprop=iS(iwall,sm_lv,imass,r);
		  for(int i=0;i<loc_vol;i++) put_spincolor_into_colorspinspin(S[iprop][i],reco_solution[r][i],id);
		}
	    }
	}
    }
  
  for(int sm_lv=0;sm_lv<nsm_lv;sm_lv++)
    for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
      for(int imass=0;imass<nmass;imass++) //put the (1+ig5)/sqrt(2) factor
	{
	  int iprop=iS(iwall,sm_lv,imass,r);
	  rotate_vol_colorspinspin_to_physical_basis(S[iprop],!r,!r);
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
      //tsource[igamma]=base_gamma[0]; //g5*g5
      dirac_prod(&(tsink[igamma]),&(base_gamma[5]),&(base_gamma[igamma]));
    }
  
  //Call the routine which does the real contraction for the Mezzotto and the Otto
  trace_id_sdag_g_s_id_sdag_g_s(contr_otto,SL1,tsink,SL2,SR1,tsink,SR2,16);
  sum_trace_id_sdag_g_s_times_trace_id_sdag_g_s(contr_mezzotto,SL1,tsink,SL2,SR1,tsink,SR2,16);
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
  trace_g_sdag_g_s(&contr_2pts,t1,s1,t2,s2,ncontr_2pts);
}

//print all the passed contractions to the file
void print_ottos_contractions_to_file(FILE *fout)
{
  double norm=glb_size[1]*glb_size[2]*glb_size[3];
  
  if(rank==0)
    for(int icontr=0;icontr<16;icontr++)
      {
        fprintf(fout,"\n");
	print_contraction_to_file(fout,icontr,5,contr_mezzotto+icontr*glb_size[0],twall[0],"DISCONNECTED ",norm);
        fprintf(fout,"\n");
	print_contraction_to_file(fout,icontr,5,contr_otto+icontr*glb_size[0],twall[0],"CONNECTED ",norm);
      }
}

//print all the passed contractions to the file
void print_two_points_contractions_to_file(FILE *fout)
{
  double norm=glb_size[1]*glb_size[2]*glb_size[3];
  
  if(rank==0)
    for(int icontr=0;icontr<ncontr_2pts;icontr++)
      {
        fprintf(fout,"\n");
	print_contraction_to_file(fout,op2_2pts[icontr],op1_2pts[icontr],contr_2pts+icontr*glb_size[0],twall[0],"",norm);
      }
}

//Calculate and print to file all the contractions
void calculate_all_contractions()
{
  contr_time-=take_time();
  
  //loop over smearing of the left and right wall
  for(int sm_lv_L=0;sm_lv_L<nsm_lv;sm_lv_L++)
    for(int sm_lv_R=0;sm_lv_R<nsm_lv;sm_lv_R++)
      {
	char path_bag[1024];
	sprintf(path_bag,"%s_%d%d",basepath_bag,sm_lv_L,sm_lv_R);
	FILE *fout_bag=open_text_file_for_output(path_bag);
	
	//loop over left-right wall combo
	for(int iwL=0;iwL<nwall;iwL++)
	  for(int iwR=iwL+1;iwR<nwall;iwR++)
	    {
	      int tsepar=abs(twall[iwL]-twall[iwR]);
	      
	      MPI_Barrier(MPI_COMM_WORLD);
	      if(rank==0) fprintf(fout_bag," # LEFT_WALL_t=%d , RIGHT_WALL_t=%d , tseparat=%d\n\n",
		      twall[iwL],twall[iwR],tsepar);
	      
	      //Ottos
	      for(int im2=0;im2<nmass;im2++)
		for(int r2=0;r2<2;r2++)
		  for(int im1=0;im1<nspec;im1++)
		    for(int r1=0;r1<2;r1++)
		      for(int r3=0;r3<2;r3++)
			{
			  int r4=1-(r1+r2+r3)%2;
			  
			  if(rank==0)
			    {
			      fprintf(fout_bag," # m1=%f r1=%d , m2=%f r2=%d , m3=%f r3=%d , m4=%f r4=%d",
				      mass[im1],r1,mass[im2],r2,mass[im1],r3,mass[im2],r4);
			      fprintf(fout_bag," . Tseparation=%d\n",tsepar);
			    }
			  
			  int ip1=iS(iwL,sm_lv_L,im1,r1),ip2=iS(iwL,sm_lv_L,im2,r2);
			  int ip3=iS(iwR,sm_lv_R,im1,r3),ip4=iS(iwR,sm_lv_R,im2,r4);
			  
			  Bk_eights(S[ip1],S[ip2],S[ip3],S[ip4]);
			  
			  ncontr_tot+=32;
			  
			  if(rank==0)
			    {
			      print_ottos_contractions_to_file(fout_bag);
			      fprintf(fout_bag,"\n");
			    }
			}
	    }
	
	if(rank==0) fclose(fout_bag);
      }

  //loop over smearing of the source and sink
  for(int sm_lv_sour=0;sm_lv_sour<nsm_lv;sm_lv_sour++)
    for(int sm_lv_sink=0;sm_lv_sink<nsm_lv;sm_lv_sink++)
      {
	char path_2pts[1024];
	sprintf(path_2pts,"%s_%d%d",basepath_2pts,sm_lv_sour,sm_lv_sink);
	FILE *fout_2pts=open_text_file_for_output(path_2pts);

	//smear all the sink
	for(int iwall=0;iwall<nwall;iwall++)
	  for(int r=0;r<2;r++)
	    for(int imass=0;imass<nmass;imass++)
	      for(int id=0;id<4;id++)
		{
		  int sm_to_app=((sm_lv_sink==0) ? 
				 jacobi_niter[sm_lv_sink] : (jacobi_niter[sm_lv_sink]-jacobi_niter[sm_lv_sink-1])); 
		  
		  //smerd each source with the level of smerding of the sink
		  int iprop=iS(iwall,sm_lv_sour,imass,r);
		  
		  for(int ivol=0;ivol<loc_vol;ivol++)
		    get_spincolor_from_colorspinspin(reco_solution[0][ivol],S[iprop][ivol],id);
		  if(rank==0 && debug) printf("Iwall %d, r=%d, imass=%d, id=%d ",iwall,r,imass,id);
		  jacobi_smearing(reco_solution[1],reco_solution[0],sme_conf,jacobi_kappa,sm_to_app);
		  for(int ivol=0;ivol<loc_vol;ivol++)
		    put_spincolor_into_colorspinspin(S[iprop][ivol],reco_solution[1][ivol],id);
              }
		
	//two points
	for(int im2=0;im2<nmass;im2++)
	  for(int r2=0;r2<2;r2++)
	    for(int im1=0;im1<nspec;im1++)
	      for(int r1=0;r1<2;r1++)
		for(int iwall=0;iwall<nwall;iwall++)
		  {
		    if(rank==0)
		      {
		       fprintf(fout_2pts," # m1=%f r1=%d , m2=%f r2=%d , wall=%d , ",mass[im1],r1,mass[im2],r2,iwall);
		       fprintf(fout_2pts," sme_source=%d sme_sink=%d\n",jacobi_niter[sm_lv_sour],jacobi_niter[sm_lv_sink]);
		      }
		    
		    int iprop1=iS(iwall,sm_lv_sour,im1,r1);
		    int iprop2=iS(iwall,sm_lv_sour,im2,r2);
		    
		    meson_two_points(S[iprop1],S[iprop2]);
		    
		    if(rank==0)
		      {
			print_two_points_contractions_to_file(fout_2pts);
			fprintf(fout_2pts,"\n");
		      }
		    ncontr_tot+=ncontr_2pts;
		  }
	if(rank==0) fclose(fout_2pts);
      }
  
  contr_time+=take_time();
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
      {
	fprintf(stderr,"Use: %s input_file\n",arg[0]);
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD,1);
      }

  //Inner initialization
  tot_time-=take_time();
  initialize_Bk(arg[1]);

  //Loop over configurations
  for(int iconf=0;iconf<ngauge_conf;iconf++)
  {
      //Read the conf
      read_str(conf_path,1024);
      load_gauge_conf();

      //Read position of first wall, reconstruct the others
      read_int(twall);
      for(int iwall=1;iwall<nwall;iwall++) twall[iwall]=(twall[0]+tsepa[iwall-1])%glb_size[0];

      //Read outpaths
      read_str(basepath_bag,1024);
      read_str(basepath_2pts,1024);

      //Invert propagators
      if(rank==0) printf("Going to invert: %d walls\n",nwall);
      for(int iwall=0;iwall<nwall;iwall++)
	{
	  generate_source(iwall);
	  calculate_S(iwall);
	}
      
      //Perform all the contractions
      calculate_all_contractions();
  }
  
  tot_time+=take_time();

  //Finalization
  close_input();
  close_Bk();
  
  return 0;
}
