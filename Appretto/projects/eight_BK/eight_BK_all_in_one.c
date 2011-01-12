
//This program calculates a list of three-point functions for the BK 
//contracting all the propagators in the four lists only if 
//if one of the r is different. 
// The first and third list are uploaded in blocks as big as possible 
// propagators  the second and fourth are uploaded one by one

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

void BK_eight_disconnected(complex **mezzo, int *list_opQ,colorspinspin *sdL,colorspinspin *ssL, colorspinspin *sdR,colorspinspin *ssR,int ncontr,int fdL,int rdL,int fsL,int rsL,int fdR,int rdR,int fsR,int rsR )
{
  
  //Temporary vectors for the internal gamma
  dirac_matr gsource_times_g5_L[ncontr],g5_times_gQ_L[ncontr];
  dirac_matr gsource_times_g5_R[ncontr],g5_times_gQ_R[ncontr];
  int gsource=5; //Kaon-Kaon
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Put the two gamma5 needed for the revert of the d spinor
      dirac_prod(&(gsource_times_g5_L[icontr]),&(base_gamma[gsource]),&(base_gamma[5]));
      dirac_prod(&(g5_times_gQ_L[icontr]),&(base_gamma[5]),&(base_gamma[list_opQ[icontr]]));
      dirac_prod(&(gsource_times_g5_R[icontr]),&(base_gamma[gsource]),&(base_gamma[5]));
      dirac_prod(&(g5_times_gQ_R[icontr]),&(base_gamma[5]),&(base_gamma[list_opQ[icontr]]));
    }
  
  //Call for the routine which does the real contraction for the Mezzotto and the otto:
  sum_trace_g_sdag_g_s_times_trace_g_sdag_g_s(mezzo,gsource_times_g5_L,sdL,g5_times_gQ_L,ssL,gsource_times_g5_R,sdR,g5_times_gQ_R,ssR,ncontr);
}

void BK_eight_connected(complex **otto, int *list_opQ,colorspinspin *sdL,colorspinspin *ssL, colorspinspin *sdR,colorspinspin *ssR,int ncontr,int fdL,int rdL,int fsL,int rsL,int fdR,int rdR,int fsR,int rsR)
{
  //Temporary vectors for the internal gamma
  dirac_matr gsource_times_g5_L[ncontr],g5_times_gQ_L[ncontr];
  dirac_matr gsource_times_g5_R[ncontr],g5_times_gQ_R[ncontr];
  int gsource=5; //Kaon-Kaon
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Put the two gamma5 needed for the revert of the d spinor
      dirac_prod(&(gsource_times_g5_L[icontr]),&(base_gamma[gsource]),&(base_gamma[5]));
      dirac_prod(&(g5_times_gQ_L[icontr]),&(base_gamma[5]),&(base_gamma[list_opQ[icontr]]));
      dirac_prod(&(gsource_times_g5_R[icontr]),&(base_gamma[gsource]),&(base_gamma[5]));
      dirac_prod(&(g5_times_gQ_R[icontr]),&(base_gamma[5]),&(base_gamma[list_opQ[icontr]]));
    }
  
  //Call for the routine which does the real contraction for the Mezzotto and the otto:
  trace_g_sdag_g_s_g_sdag_g_s(otto,gsource_times_g5_L,sdL,g5_times_gQ_L,ssL,gsource_times_g5_R,sdR,g5_times_gQ_R,ssR,ncontr);
}

#include "appretto.h"

char conf_path[1024];
quad_su3 *conf;
double kappa;

int nmass;
double *mass;
double put_theta[4],old_theta[4]={0,0,0,0};

//source data
int seed,noise_type,twall[2];
double rad2=1.414213562373095048801688724209;
spincolor *source;
colorspinspin *original_source;

//vectors for the spinor data
colorspinspin **S[2][2];
spincolor **cgmms_solution,*reco_solution[2];

//cgmms inverter parameters
double stopping_residue;
int stopping_criterion;
int niter_max;

//two points contractions
int ncontr_2pts;
complex **contr_2pts;
int *op1_2pts,*op2_2pts;
char outfile_2pts[1024];

//timings
int ninv_tot=0,ncontr_tot=0;
double tot_time=0,inv_time=0,contr_time=0;

//This function takes care to make the revert on the FIRST spinor, putting the needed gamma5
void meson_two_points(complex **corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Put the two gamma5 needed for the revert of the first spinor
      dirac_prod(&(t1[icontr]), &(base_gamma[list_op1[icontr]]),&(base_gamma[5]));
      dirac_prod(&(t2[icontr]), &(base_gamma[5]),&(base_gamma[list_op2[icontr]]));
    }
  
  //Call the routine which does the real contraction
  trace_g_sdag_g_s(corr,t1,s1,t2,s2,ncontr);
}

//Generate the source for the dirac index
void generate_source(int LR)
{ //reset
  memset(original_source,0,sizeof(colorspinspin)*loc_vol);
  
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    if(glb_coord_of_loclx[loc_site][0]==twall[LR])
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
  //Gauge path
  read_str_str("GaugeConfPath",conf_path,1024);
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the source
  
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  init_random(seed);
  //Read the location of the left wall
  read_str_int("TLeftWall",&(twall[0]));
  read_str_int("TRightWall",&(twall[1]));
  //Read the noise type
  read_str_int("NoiseType",&noise_type);

  // 3) Read list of masses

  read_list_of_doubles("NMass",&nmass,&mass);

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
  
  //Number of iterations
  read_str_int("NiterMax",&niter_max);
  
  // 5) contraction list for two points

  read_str_int("NContrTwoPoints",&ncontr_2pts);
  contr_2pts=(complex**)malloc(sizeof(complex*)*ncontr_2pts);
  contr_2pts[0]=(complex*)malloc(sizeof(complex)*ncontr_2pts*glb_size[0]); 
  op1_2pts=(int*)malloc(sizeof(int)*ncontr_2pts);
  op2_2pts=(int*)malloc(sizeof(int)*ncontr_2pts);
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      contr_2pts[icontr]=contr_2pts[0]+icontr*glb_size[0];

      //Read the operator pairs
      read_int(&(op1_2pts[icontr]));
      read_int(&(op2_2pts[icontr]));

      if(rank==0 && debug) printf(" contr.%d %d %d\n",icontr,op1_2pts[icontr],op2_2pts[icontr]);
    }
  
  read_str_str("OutfileTwoPoints",outfile_2pts,1024);

  close_input();

  ////////////////////////////////////// end of input reading/////////////////////////////////

  //allocate gauge conf, Pmunu and all the needed spincolor and colorspinspin
  conf=allocate_quad_su3(loc_vol+loc_bord,"conf");

  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_local_gauge_conf(conf,conf_path);
  communicate_gauge_borders(conf);
  
  double gplaq=global_plaquette(conf);
  if(rank==0) printf("plaq: %.18g\n",gplaq);
  
  //Put the anti-periodic condition on the temporal border
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,0);

  //Allocate all the propagators colorspinspin vectors
  for(int LR=0;LR<2;LR++)
    for(int UD=0;UD<2;UD++)
      {
	S[LR][UD]=(colorspinspin**)malloc(sizeof(colorspinspin*)*nmass);
	for(int iprop=0;iprop<nmass;iprop++) S[LR][UD][iprop]=allocate_colorspinspin(loc_vol,"S");
      }
  
  //Allocate nmass spincolors, for the cgmms solutions
  cgmms_solution=(spincolor**)malloc(sizeof(spincolor*)*nmass);
  for(int imass=0;imass<nmass;imass++) cgmms_solution[imass]=allocate_spincolor(loc_vol+loc_bord,"cgmms_solution");
  reco_solution[0]=allocate_spincolor(loc_vol,"reco_solution[0]");
  reco_solution[1]=allocate_spincolor(loc_vol,"reco_solution[1]");
  
  //Allocate one spincolor for the source
  source=allocate_spincolor(loc_vol+loc_bord,"source");

  //Allocate the colorspinspin of the LR source
  original_source=allocate_colorspinspin(loc_vol,"original_source");
}

//Finalization
void close_semileptonic()
{
  if(rank==0)
    {
      printf("\n");
      printf("Total time: %g, of which:\n",tot_time);
      printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
      printf(" - %02.2f%s to perform %d contr. (%2.2gs avg)\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
    }

  free(conf);free(mass);
  for(int iprop=0;iprop<nmass;iprop++) for(int LR=0;LR<2;LR++) for(int UD=0;UD<2;UD++) free(S[LR][UD][iprop]);
  free(reco_solution[0]);free(reco_solution[1]);
  //free(op1_2pts);free(op2_2pts);free(ch_op1_2pts);free(ch_op2_2pts);
  //free(contr_2pts[0]);free(contr_2pts);free(ch_contr_2pts[0]);free(ch_contr_2pts);
  //free(op1_3pts);free(op2_3pts);free(ch_op1_3pts);free(ch_op2_3pts);
  //free(contr_3pts[0]);free(contr_3pts);free(ch_contr_3pts[0]);free(ch_contr_3pts);
  close_appretto();
}

//calculate the standard propagators
void calculate_S(int LR)
{
  for(int id=0;id<4;id++)
    { //loop over the source dirac index
      for(int ivol=0;ivol<loc_vol;ivol++)
	{
	  get_spincolor_from_colorspinspin(source[ivol],original_source[ivol],id);
	  //put the g5
	  for(int id1=2;id1<4;id1++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic][ri]*=-1;
	}
      
      communicate_lx_spincolor_borders(source);
      
      double part_time=-take_time();
      inv_Q2_cgmms(cgmms_solution,source,NULL,conf,kappa,mass,nmass,niter_max,stopping_residue,stopping_criterion);
      part_time+=take_time();ninv_tot++;inv_time+=part_time;
      char tag[2][10]={"left","right"};
      if(rank==0) printf("Finished the %s source inversion, dirac index %d in %g sec\n",tag[LR],id,part_time);
      
      for(int imass=0;imass<nmass;imass++)
	{ //reconstruct the doublet
	  reconstruct_doublet(reco_solution[0],reco_solution[1],cgmms_solution[imass],conf,kappa,mass[imass]);
	  if(rank==0) printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	  for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
	    for(int i=0;i<loc_vol;i++)
	      put_spincolor_into_colorspinspin(S[LR][r][imass][i],reco_solution[r][i],id);
	}
    }
  
  for(int LR=0;LR<2;LR++)
    for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
      for(int ipropS=0;ipropS<nmass;ipropS++) //put the (1+ig5)/sqrt(2) factor
	rotate_vol_colorspinspin_to_physical_basis(S[LR][r][ipropS],!r,!r);
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

  tot_time-=take_time();
  initialize_Bk(arg[1]);
  
  for(int LR=0;LR<2;LR++)
    {
      generate_source(LR);
      calculate_S(LR);
    }
  
  tot_time+=take_time();
  close_semileptonic();
  
  return 0;
}
