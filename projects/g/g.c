#include "nissa.h"

/*
  This code is a simplified and improved version of the "semileptonic_smeared" program.
  It computes only source-smeared standing two points for light-light, light-charm and charm-charm mesons,
  and source-sink smeared moving three points for light-charm and charm-light mesons.
  Only one theta is taken as input.
  This makes the program suitable for the computation of g(D*->Dp) and g(D*D\gamma)) form factors.
  
  This in particular are the three point:
  
  |          / \	   |          / \	    |   
  |  c,th=0 /   \ c,th!=0  |  l,th=0 /   \ l,th!=0  |
  |        /     \	   |        /     \	    |   
  |       /_l,th0_\ g5     |       /_c,th0_\ g5     |

 */

//fix the r of the spectator
int r_spec;
int sign_r[2]={-1,+1};

//gauge info
int ngauge_conf,nanalized_conf=0;
char conf_path[1024];
quad_su3 *conf,*sme_conf;
double kappa;
char outfolder[1024];
double put_theta[4],old_theta[4]={0,0,0,0};

//source data
int seed,noise_type,twall;
spincolor *source;
colorspinspin *original_source;

//masses and theta
double lmass,cmass;
double theta;

//smearing parameters
double jacobi_kappa,ape_alpha;
int ape_niter;
int jacobi_niter;

//vectors for the spinor data
colorspinspin *S0L,*S0C;
colorspinspin *S1LC,*S1CL;

//inverter parameters
double stopping_residue_l,stopping_residue_c;
int niter_max;

//two points contractions
int ncontr_2pts;
complex *contr_2pts;
int *op1_2pts,*op2_2pts;

//sequential props info
int tsep;
colorspinspin *sequential_source;

//three points contractions
int ncontr_3pts;
complex *contr_3pts;
int *op1_3pts,*op2_3pts;

//timings
int ninv_tot=0,ncontr_tot=0;
int wall_time;
double tot_time=0,inv_time=0;
double load_time=0,contr_time;
double contr_save_time=0;
double contr_3pts_time=0;
double contr_2pts_time=0;

//smear a colorspinspin
void smear_colorspinspin(colorspinspin *out,colorspinspin *in)
{
  spincolor *temp=nissa_malloc("temp",loc_vol+loc_bord,spincolor);
  
  for(int id=0;id<4;id++)
    {
      for(int ivol=0;ivol<loc_vol;ivol++)
	get_spincolor_from_colorspinspin(temp[ivol],in[ivol],id);
      
      jacobi_smearing(temp,temp,sme_conf,jacobi_kappa,jacobi_niter);
      
      for(int ivol=0;ivol<loc_vol;ivol++)
	put_spincolor_into_colorspinspin(out[ivol],temp[ivol],id);
    }
  
  nissa_free(temp);
}

//This function takes care to make the revert on the FIRST spinor, putting the needed gamma5
void meson_two_points(complex *corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr)
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

//This function contract a source with a sequential spinor putting the passed list of operators
void contract_with_source(complex *corr,colorspinspin *S1,int *list_op,colorspinspin *source)
{
  //Temporary vector for the internal matrices
  dirac_matr t1[ncontr_2pts],t2[ncontr_2pts];

  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    { //Init the second 
      t1[icontr]=base_gamma[0];
      t2[icontr]=base_gamma[list_op[icontr]];
    }

  //Call the routine which does the real contraction
  trace_g_sdag_g_s(corr,t1,S1,t2,source,ncontr_2pts);
}

//generate the source
void generate_source()
{
  enum rnd_type type[5]={RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_Z2,RND_Z2,RND_Z4};
  generate_spindiluted_source(original_source,type[noise_type],twall);
  smear_colorspinspin(original_source,original_source);
}

//Generate a sequential source for S1
void generate_sequential_source(colorspinspin *S0)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      //put to zero everything but the slice
      if(glb_coord_of_loclx[ivol][0]!=(twall+tsep)%glb_size[0])
	memset(sequential_source[ivol],0,sizeof(colorspinspin));
      else
	{
	  //take ths S0
	  memcpy(sequential_source[ivol],S0[ivol],sizeof(colorspinspin));
	  for(int c=0;c<3;c++)
	    {
	      //rotate as r_spec because it's D^-1
	      rotate_spinspin_to_physical_basis(sequential_source[ivol][c],r_spec,r_spec);
	      //put g5
	      for(int id_sink=2;id_sink<4;id_sink++)
		for(int id_sour=0;id_sour<4;id_sour++)
		  for(int ri=0;ri<2;ri++)
		    sequential_source[ivol][c][id_sink][id_sour][ri]*=-1;
	    }
	}
    }
  //smear the seq source (twice!)
  smear_colorspinspin(sequential_source,sequential_source);
  smear_colorspinspin(sequential_source,sequential_source);
  
  master_printf("Sequential source created!\n\n");
}  

//Parse all the input file
void initialize_semileptonic(char *input_path)
{
  tot_time-=take_time();
  
  open_input(input_path);

  // 1) Read information about the gauge conf
  
  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Wall time
  read_str_int("WallTime",&wall_time);
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the source
  
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  //Read the noise type
  read_str_int("NoiseType",&noise_type);

  // 3) Smearing parameters
  
  //Smearing parameters
  read_str_double("ApeAlpha",&ape_alpha);
  read_str_int("ApeNiter",&ape_niter);
  read_str_double("JacobiKappa",&jacobi_kappa);
  read_str_int("JacobiNiter",&jacobi_niter);
  
  // 4) Read list of masses and of thetas

  read_str_double("lMass",&lmass);
  read_str_double("cMass",&cmass);
  read_str_double("Theta",&theta);

  // 5) Info about the inverter
  
  //Stopping Residue
  read_str_double("StoppingResidueL",&stopping_residue_l);
  read_str_double("StoppingResidueC",&stopping_residue_c);
  
  //Number of iterations
  read_str_int("NiterMax",&niter_max);
  
  // 6) contraction list for two points
  
  read_str_int("NContrTwoPoints",&ncontr_2pts);
  contr_2pts=nissa_malloc("contr_2pts",ncontr_2pts*glb_size[0],complex);
  op1_2pts=nissa_malloc("op1_2pts",ncontr_2pts,int);
  op2_2pts=nissa_malloc("op2_2pts",ncontr_2pts,int);
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1_2pts[icontr]));
      read_int(&(op2_2pts[icontr]));
      
      master_printf(" contr.%d %d %d\n",icontr,op1_2pts[icontr],op2_2pts[icontr]);
    }
  
  // 7) three points functions
  
  sequential_source=nissa_malloc("Sequential source",loc_vol,colorspinspin);
  read_str_int("TSep",&tsep);
  read_str_int("NContrThreePoints",&ncontr_3pts);
  contr_3pts=nissa_malloc("contr_3pts",ncontr_3pts*glb_size[0],complex); 
  op1_3pts=nissa_malloc("op1_3pts",ncontr_3pts,int);
  op2_3pts=nissa_malloc("op2_3pts",ncontr_3pts,int);
  for(int icontr=0;icontr<ncontr_3pts;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1_3pts[icontr]));
      read_int(&(op2_3pts[icontr]));
      
      master_printf(" contr.%d %d %d\n",icontr,op1_3pts[icontr],op2_3pts[icontr]);
    }

  read_str_int("NGaugeConf",&ngauge_conf);
  
  master_printf("\n");
  
  ////////////////////////////////////// end of input reading/////////////////////////////////

  //allocate gauge conf, Pmunu and all the needed spincolor and colorspinspin
  conf=nissa_malloc("or_conf",loc_vol+loc_bord+loc_edge,quad_su3);
  sme_conf=nissa_malloc("sm_conf",loc_vol+loc_bord,quad_su3);

  //Allocate all the S0 colorspinspin vectors
  S0L=nissa_malloc("S0L",loc_vol,colorspinspin);
  S0C=nissa_malloc("S0C",loc_vol,colorspinspin);
  
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+loc_bord,spincolor);
  original_source=nissa_malloc("original_source",loc_vol,colorspinspin);

  //Allocate all the S1 colorspinspin vectors
  S1LC=nissa_malloc("S1LC",loc_vol,colorspinspin);
  S1CL=nissa_malloc("S1CL",loc_vol,colorspinspin);
}

//find a new conf
int read_conf_parameters(int *iconf)
{
  int ok_conf;

  do
    {
      //Gauge path
      read_str(conf_path,1024);
      
      //Twall
      read_int(&(twall));
      
      //Out folder
      read_str(outfolder,1024);
      
      //Check if the conf exist
      master_printf("Considering configuration %s\n",conf_path);
      ok_conf=!(dir_exists(outfolder));
      if(ok_conf)
	{
	  int ris=create_dir(outfolder);
	  if(ris==0) master_printf("Configuration %s not already analized, starting.\n",conf_path);
	  else
	    {
	      ok_conf=0;
	      master_printf("Configuration %s taken by someone else.\n",conf_path);
	    }
	}
      else
	master_printf("Configuration %s already analized, skipping.\n",conf_path);
      (*iconf)++;
    }
  while(!ok_conf && (*iconf)<ngauge_conf);
  
  master_printf("\n");
  
  return ok_conf;
}

//read the conf and setup it
void setup_conf()
{
  //load the gauge conf, propagate borders
  read_ildg_gauge_conf(conf,conf_path);
  communicate_lx_gauge_borders(conf);
  
  //prepare the smerded version
  ape_smearing(sme_conf,conf,ape_alpha,ape_niter);
  communicate_lx_gauge_borders(sme_conf);
  
  //compute plaquette
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  master_printf("smerded plaq: %.18g\n",global_plaquette_lx_conf(sme_conf));
  
  //put the anti-periodic condition on the temporal border
  memset(old_theta,0,4*sizeof(double));
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,1);
}

//Finalization
void close_semileptonic()
{
  tot_time+=take_time();

  contr_time=contr_2pts_time+contr_3pts_time+contr_save_time;
  
  master_printf("\n");
  master_printf("Total time: %g, of which:\n",tot_time);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d contr. (%2.2gs avg) of which:\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
  master_printf("   * %02.2f%s to compute two points\n",contr_2pts_time*100.0/contr_time,"%");
  master_printf("   * %02.2f%s to compute three points\n",contr_3pts_time*100.0/contr_time,"%");
  master_printf("   * %02.2f%s to save correlations\n",contr_save_time*100.0/contr_time,"%");
  
  nissa_free(conf);nissa_free(sme_conf);
  nissa_free(S0L);nissa_free(S0C);
  nissa_free(S1LC);nissa_free(S1CL);
  nissa_free(sequential_source);
  nissa_free(contr_2pts);nissa_free(contr_3pts);
  nissa_free(op1_2pts);nissa_free(op2_2pts);
  nissa_free(op1_3pts);nissa_free(op2_3pts);
  nissa_free(source);nissa_free(original_source);
  close_nissa();
}

//calculate the standard propagators
void calculate_S0(colorspinspin *S0,double mass,double stopping_residue)
{
  //temporary vector for sol
  spincolor *temp_sol=nissa_malloc("Temp_sol",loc_vol,spincolor);
  
  //loop over the source dirac index
  for(int id=0;id<4;id++)
    { 
      for(int ivol=0;ivol<loc_vol;ivol++)
	get_spincolor_from_colorspinspin(source[ivol],original_source[ivol],id);
      
      //inverting the light quark
      double part_time=-take_time();
      inv_tmD_cg_eoprec_eos(temp_sol,source,NULL,conf,kappa,sign_r[r_spec]*mass,niter_max,stopping_residue);
      part_time+=take_time();ninv_tot++;inv_time+=part_time;
      master_printf("Finished the inversion of S0 mass=%lg, dirac index %d in %g sec\n",mass,id,part_time);

      //convert the id-th spincolor into the colorspinspin
      for(int i=0;i<loc_vol;i++)
	  put_spincolor_into_colorspinspin(S0[i],temp_sol[i],id);
    }
  
  //rotate to physical basis; remember that D^-1 rotate opposite than D!
  rotate_vol_colorspinspin_to_physical_basis(S0,!r_spec,!r_spec);
  
  master_printf("Propagators rotated\n");
  
  master_printf("\n");
  
  nissa_free(temp_sol);
}

//calculate the sequential propagators
void calculate_S1(colorspinspin *S1,double mass,colorspinspin *S0,double stopping_residue)
{
  generate_sequential_source(S0);
  
  //temporary solution
  spincolor *temp_sol=nissa_malloc("Temp_sol",loc_vol,spincolor);
  
  //loop over dirac index of the source
  for(int id=0;id<4;id++)
    { 
      for(int ivol=0;ivol<loc_vol;ivol++)
	get_spincolor_from_colorspinspin(source[ivol],sequential_source[ivol],id);
      communicate_lx_spincolor_borders(source);

      double part_time=-take_time();
      inv_tmD_cg_eoprec_eos(temp_sol,source,NULL,conf,kappa,-sign_r[r_spec]*mass,niter_max,stopping_residue);
      part_time+=take_time();ninv_tot++;inv_time+=part_time;
      master_printf("Finished the inversion of S1 mass=%lg dirac index %d in %g sec\n",mass,id,part_time);
      
      //put in the solution
      for(int i=0;i<loc_vol;i++)
	put_spincolor_into_colorspinspin(S1[i],temp_sol[i],id);
    }
  
  //put the (1+-ig5)/sqrt(2) factor. On the source rotate as r_spec, on the sink as !r_spec
  //but, being D^-1, everything is swapped
  rotate_vol_colorspinspin_to_physical_basis(S1,!r_spec,r_spec);
  master_printf("Propagators rotated\n\n");

  nissa_free(temp_sol);
}

//Calculate and print to file the 2pts
void calculate_all_2pts(colorspinspin *SA,colorspinspin *SB,double SA_mass,double SB_mass,const char *tag)
{
  contr_2pts_time-=take_time();

  meson_two_points(contr_2pts,op1_2pts,SA,op2_2pts,SB,ncontr_2pts);
  ncontr_tot+=ncontr_2pts;
  
  char path[1024];
  sprintf(path,"%s/2pts_%s",outfolder,tag);
  FILE *fout=open_text_file_for_output(path);

  if(rank==0) fprintf(fout," # mA=%f, mB=%f\n",SA_mass,SB_mass);
      
  contr_save_time-=take_time();
  print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,twall,"",1.0);
  contr_save_time+=take_time();
  
  if(rank==0)
    {
      fprintf(fout,"\n");
      fclose(fout);
    }
  
  contr_2pts_time+=take_time();
}

//Calculate and print to file the 3pts
void calculate_all_3pts(colorspinspin *S0,colorspinspin *S1,double S0_mass,double S1_mass,const char *tag)
{
  contr_3pts_time-=take_time();
  
  meson_two_points(contr_3pts,op1_3pts,S0,op2_3pts,S1,ncontr_3pts);
  ncontr_tot+=ncontr_3pts;
  
  char path[1024];
  sprintf(path,"%s/3pts_%s",outfolder,tag);
  FILE *fout=open_text_file_for_output(path);
  
  if(rank==0) fprintf(fout," # S0_mass=%f , S1_mass=%f\n",S0_mass,S1_mass);
  
  contr_save_time-=take_time();
  print_contractions_to_file(fout,ncontr_3pts,op1_3pts,op2_3pts,contr_3pts,twall,"",1.0);
  contr_save_time+=take_time();
  
  if(rank==0)
    {
      fprintf(fout,"\n");
      fclose(fout);
    }
  
  contr_3pts_time+=take_time();
}

//check if the time is enough
int check_remaining_time()
{
  int enough_time;

  //check remaining time
  double temp_time=take_time()+tot_time;
  double ave_time=temp_time/nanalized_conf;
  double left_time=wall_time-temp_time;
  enough_time=left_time>(ave_time*1.1);

  master_printf("Remaining time: %lg sec\n",left_time);
  master_printf("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
  if(enough_time) master_printf("Continuing with next conf!\n");
  else master_printf("Not enough time, exiting!\n");
  
  return enough_time;
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();

  if(narg<3) crash("Use: %s 0,1 input_file",arg[0]);
  r_spec=atoi(arg[1]);
  master_printf("R_spec: %d\n",r_spec);
  
  initialize_semileptonic(arg[2]);
  
  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf && enough_time && read_conf_parameters(&iconf))
    {
      setup_conf();
      generate_source();
      
      ////////////////// standing two points /////////////////
      
      calculate_S0(S0L,lmass,stopping_residue_l);
      calculate_S0(S0C,cmass,stopping_residue_c);
      
      calculate_all_2pts(S0L,S0L,lmass,lmass,"ll_sl");
      calculate_all_2pts(S0L,S0C,lmass,cmass,"lc_sl");
      calculate_all_2pts(S0C,S0C,cmass,cmass,"cc_sl");
      
      ////////////////// moving three points /////////////////
      
      put_theta[1]=put_theta[2]=put_theta[3]=theta;
      adapt_theta(conf,old_theta,put_theta,1,1);
      
      calculate_S1(S1LC,cmass,S0L,stopping_residue_c);
      calculate_S1(S1CL,lmass,S0C,stopping_residue_l);
      
      calculate_all_3pts(S0L,S1CL,lmass,lmass,"charm_spec");
      calculate_all_3pts(S0C,S1LC,cmass,cmass,"light_spec");
      
      /////////////// standing smeared two points ////////////
      
      smear_colorspinspin(S0L,S0L);
      smear_colorspinspin(S0C,S0C);
      
      calculate_all_2pts(S0L,S0L,lmass,lmass,"ll_ss");
      calculate_all_2pts(S0L,S0C,lmass,cmass,"lc_ss");
      calculate_all_2pts(S0C,S0C,cmass,cmass,"cc_ss");
      
      ////////////////////////////////////////////////////////
      
      nanalized_conf++;      
      enough_time=check_remaining_time();
    }
  
  if(iconf==ngauge_conf) master_printf("Finished all the conf!\n");
  
  close_semileptonic();
  
  return 0;
}
