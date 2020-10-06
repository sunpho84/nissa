#include <stdlib.h>

#include "nissa.hpp"
using namespace nissa;

/*
  This code is a simplified and improved version of the "semileptonic_smeared" program.
  Only one theta per heavy and compute automatically all the +-theta pairs
  It needs to compute both l(0)h(th) and l(th)h(0) to eliminate Z in three pts
  
  This in particular are the three points:
  
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
spincolor *source,*temp_sol;
colorspinspin *original_source;

//masses and theta
int nmass,nheavy,ntheta;
double *mass;
double *theta;

//smearing parameters
double gaussian_kappa,ape_alpha;
int ape_niter;
int gaussian_niter;

//vectors for the spinor data
colorspinspin **S0L,**S0H;
colorspinspin **S1LH,**S1HL;

//inverter parameters
double *stopping_residue;
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
int ninv_tot=0,ncontr_tot=0,nsmear_tot=0;
int wall_time;
double tot_prog_time=0,inv_time=0;
double load_time=0,contr_time;
double smear_time=0;
double contr_save_time=0;
double contr_3pts_time=0;
double contr_2pts_time=0;

//output files
FILE *fout_2pts_sl;
FILE *fout_2pts_ss;
FILE *fout_3pts;

int iS0L(int mp,int itheta)
{
  switch(mp)
    {
    case  0: return 0;break;
    case -1: return 1+0+2*itheta;break;
    case +1: return 1+1+2*itheta;break;
    default: crash("unknown mp: %d",mp);return 0;break;
    }
}

int iS0H(int m0p,int itheta)
{
  m0p++;
  return m0p+3*itheta;
}

int iS1(int m0p,int itheta)
{
  m0p++;
  return m0p+3*itheta;
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
  complex *loc_corr=nissa_malloc("loc_corr",ncontr*glb_size[0],complex);
  trace_g_css_dag_g_css(corr,loc_corr,t1,s1,t2,s2,ncontr);
  nissa_free(loc_corr);
}

//generate the source
void generate_source()
{
  enum rnd_t type[5]={RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_Z2,RND_Z2,RND_Z4};
  generate_spindiluted_source(original_source,type[noise_type],twall);
  gaussian_smearing(original_source,original_source,sme_conf,gaussian_kappa,gaussian_niter);
}

//Generate a sequential source for S1
void generate_sequential_source(colorspinspin *S0)
{
  NISSA_LOC_VOL_LOOP(ivol)
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
  set_borders_invalid(sequential_source);
  
  //smear the seq source (twice!)
  gaussian_smearing(sequential_source,sequential_source,sme_conf,gaussian_kappa,gaussian_niter*2);
  
  master_printf("Sequential source created!\n\n");
}  

//Parse all the input file
void initialize_semileptonic(char *input_path)
{
  tot_prog_time-=take_time();
  
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
  //RSpec
  read_str_int("RSpec",&r_spec);

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
  read_str_double("GaussianKappa",&gaussian_kappa);
  read_str_int("GaussianNiter",&gaussian_niter);
  
  // 4) Read list of masses and of thetas, as well as niter_max

  read_list_of_double_triples("MassThetaResidue",&nmass,&mass,&theta,&stopping_residue);
  if(theta[0]!=0) crash("theta for light cannot be %lg",theta[0]);
  read_str_int("NiterMax",&niter_max);
  nheavy=nmass-1;
  
  //suppress 0 theta (AAAAARGHHHH!!!!)
  theta++;
  ntheta=nheavy;
  
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
  
  sequential_source=nissa_malloc("Sequential source",loc_vol+bord_vol,colorspinspin);
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
  conf=nissa_malloc("or_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  sme_conf=nissa_malloc("sm_conf",loc_vol+bord_vol,quad_su3);

  //Allocate all the S0 colorspinspin vectors
  S0L=nissa_malloc("S0L*",2*nheavy+1,colorspinspin*);
  S0H=nissa_malloc("S0H*",3*nheavy,colorspinspin*);
  for(int imass=0;imass<2*nheavy+1;imass++) S0L[imass]=nissa_malloc("S0L",loc_vol,colorspinspin);
  for(int imass=0;imass<3*nheavy;imass++) S0H[imass]=nissa_malloc("S0H",loc_vol,colorspinspin);
  
  S1LH=nissa_malloc("S1LH*",3*nheavy,colorspinspin*);
  S1HL=nissa_malloc("S1HL*",3*nheavy,colorspinspin*);
  for(int imass=0;imass<3*nheavy;imass++)
    {
      S1LH[imass]=nissa_malloc("S1LH",loc_vol,colorspinspin);
      S1HL[imass]=nissa_malloc("S1HL",loc_vol,colorspinspin);
    }
  
  //Allocate one spincolor for source, sequential source and temporary sol
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  temp_sol=nissa_malloc("temp_sol",loc_vol,spincolor);
  original_source=nissa_malloc("original_source",loc_vol+bord_vol,colorspinspin);
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
  
  //prepare the smerded version
  ape_spatial_smear_conf(sme_conf,conf,ape_alpha,ape_niter);
  
  //compute plaquette
  master_printf("plaq: %.16g\n",global_plaquette_lx_conf(conf));
  master_printf("smerded plaq: %.16g\n",global_plaquette_lx_conf(sme_conf));
  
  //put the anti-periodic condition on the temporal border
  memset(old_theta,0,4*sizeof(double));
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,1);
  
  //open output files
  char path[1024];
  sprintf(path,"%s/2pts_sl_corr",outfolder);
  fout_2pts_sl=open_text_file_for_output(path);
  sprintf(path,"%s/2pts_ss_corr",outfolder);
  fout_2pts_ss=open_text_file_for_output(path);
  sprintf(path,"%s/3pts_corr",outfolder);
  fout_3pts=open_text_file_for_output(path);
}

//calculate the standard propagators
void calculate_S0(colorspinspin *S0,double mass,double theta,double stopping_residue)
{
  //adapt with passed theta
  put_theta[1]=put_theta[2]=put_theta[3]=theta;
  adapt_theta(conf,old_theta,put_theta,1,1);

  //loop over the source dirac index
  for(int id=0;id<4;id++)
    { 
      get_spincolor_from_colorspinspin(source,original_source,id);
      
      //inverting the light quark
      double part_time=-take_time();
      inv_tmD_cg_eoprec_eos(temp_sol,NULL,conf,kappa,sign_r[r_spec]*mass,niter_max,stopping_residue,source);
      part_time+=take_time();ninv_tot++;inv_time+=part_time;
      master_printf("Finished the inversion of S0 mass=%lg, dirac index %d in %g sec\n",mass,id,part_time);

      //convert the id-th spincolor into the colorspinspin
      put_spincolor_into_colorspinspin(S0,temp_sol,id);
    }
  
  //rotate to physical basis; remember that D^-1 rotate opposite than D!
  rotate_vol_colorspinspin_to_physical_basis(S0,!r_spec,!r_spec);
  
  master_printf("Propagators rotated\n");
  
  master_printf("\n");
}

//calculate the sequential propagators
void calculate_S1(colorspinspin *S1,double mass,double theta,colorspinspin *S0,double stopping_residue)
{
  double S1_time=-take_time();
  
  double ini_time=-take_time();
  
  //deal with sequential source
  static colorspinspin *old_S0=NULL;
  if(old_S0==NULL||S0!=old_S0) generate_sequential_source(S0);
  else master_printf("no need to generate again sequential source");
  old_S0=S0;
  
  //adapt with passed theta
  put_theta[1]=put_theta[2]=put_theta[3]=theta;
  adapt_theta(conf,old_theta,put_theta,1,1);
  
  ini_time+=take_time();
  master_printf("Initialization time: %lg s\n",ini_time);
  
  //loop over dirac index of the source
  for(int id=0;id<4;id++)
    {
      get_spincolor_from_colorspinspin(source,sequential_source,id);

      double part_time=-take_time();
      inv_tmD_cg_eoprec_eos(temp_sol,NULL,conf,kappa,-sign_r[r_spec]*mass,niter_max,stopping_residue,source);
      part_time+=take_time();ninv_tot++;inv_time+=part_time;
      master_printf("Finished the inversion of dirac index %d in %lg sec\n",id,part_time);
      
      //put in the solution
      put_spincolor_into_colorspinspin(S1,temp_sol,id);
    }
  
  //put the (1+-ig5)/sqrt(2) factor. On the source rotate as r_spec, on the sink as !r_spec
  //but, being D^-1, everything is swapped
  rotate_vol_colorspinspin_to_physical_basis(S1,!r_spec,r_spec);
  master_printf("Propagators rotated\n\n");
  
  S1_time+=take_time();
  master_printf("Finished the inversion of S1 mass=%lg theta=%lg in %lg sec\n",mass,theta,S1_time);
}

//compute all the S0 prop
void calculate_all_S0()
{
  //compute the standing S0 light prop
  calculate_S0(S0L[0],mass[0],0.0,stopping_residue[0]);
  
  //compute the -theta_i and +theta_i S0 light prop for each heavy-related theta
  for(int mp=-1;mp<=1;mp+=2)
    for(int itheta=0;itheta<ntheta;itheta++)
      calculate_S0(S0L[iS0L(mp,itheta)],mass[0],mp*theta[itheta],stopping_residue[0]);
  
  //compute the standing and moving S0 heavy prop
  for(int m0p=-1;m0p<=1;m0p++)
    for(int imass=0;imass<nheavy;imass++)
      calculate_S0(S0H[iS0H(m0p,imass)],mass[1+imass],m0p*theta[imass],stopping_residue[1+imass]);
}

//compute all the sequential propagator
void calculate_all_S1()
{
  for(int mp=-1;mp<=1;mp++)
    for(int imass=0;imass<nheavy;imass++)
      {
	calculate_S1(S1LH[iS1(mp,imass)],mass[1+imass],mp*theta[imass],S0L[iS0L(0,0)],stopping_residue[1+imass]);
	calculate_S1(S1HL[iS1(mp,imass)],mass[0],mp*theta[imass],S0H[iS0H(0,imass)],stopping_residue[0]);
      }
}

//smear all the sinks of S0
void smear_all_S0()
{
  smear_time-=time(0);
  
  master_printf("Smearing all S0 propagators\n\n");
  for(int i=0;i<2*nheavy+1;i++) gaussian_smearing(S0L[i],S0L[i],sme_conf,gaussian_kappa,gaussian_niter);
  for(int i=0;i<3*nheavy;i++) gaussian_smearing(S0H[i],S0H[i],sme_conf,gaussian_kappa,gaussian_niter);

  nsmear_tot+=2*nheavy+1+3*nheavy;
  smear_time+=time(0);
}

//Calculate and print to file the 2pts
void calculate_2pts_prop_combo(FILE *fout,colorspinspin *SA,colorspinspin *SB,double SA_mass,double SA_th,double SB_mass,double SB_th)
{
  meson_two_points(contr_2pts,op1_2pts,SA,op2_2pts,SB,ncontr_2pts);
  ncontr_tot+=ncontr_2pts;
  
  master_fprintf(fout," # mA=%lg thA=%lg, mB=%lg thB=%lg\n",SA_mass,SA_th,SB_mass,SB_th);
  
  contr_save_time-=take_time();
  print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,twall,"",1.0);
  contr_save_time+=take_time();
  
  master_fprintf(fout,"\n");
}

//calculate all the 2pts
void calculate_all_2pts(FILE *fout)
{
  contr_2pts_time-=take_time();
  
  master_printf("Computing 2pts contractions\n\n");
  
  //pion
  calculate_2pts_prop_combo(fout,S0L[iS0L(0,0)],S0L[iS0L(0,0)],mass[0],0,mass[0],0);
  for(int m0p=-1;m0p<=1;m0p++)
    for(int imass=0;imass<nheavy;imass++)
      {
	//D with moving (or at rest) c
	calculate_2pts_prop_combo(fout,S0L[iS0L(0,0)],S0H[iS0H(m0p,imass)],mass[0],0,mass[1+imass],m0p*theta[imass]);
	//D with moving (or at rest) l - certain redundant present
	calculate_2pts_prop_combo(fout,S0L[iS0L(m0p,0)],S0H[iS0H(0,imass)],mass[0],m0p*theta[imass],mass[1+imass],0);
      }
  
  contr_2pts_time+=take_time();
}

//Calculate and print to file the 3pts
void calculate_all_3pts_prop_combo(colorspinspin *S0,colorspinspin *S1,double h_mass,double theta,const char *tag)
{
  contr_3pts_time-=take_time();
  
  master_fprintf(fout_3pts," # h_mass=%lg theta=%lg %s\n",h_mass,theta,tag);
  fflush(fout_3pts);
  
  vector_reset(contr_3pts);
  
  meson_two_points(contr_3pts,op1_3pts,S0,op2_3pts,S1,ncontr_3pts);
  ncontr_tot+=ncontr_3pts;
  
  contr_save_time-=take_time();
  print_contractions_to_file(fout_3pts,ncontr_3pts,op1_3pts,op2_3pts,contr_3pts,twall,"",1.0);
  contr_save_time+=take_time();
  
  master_fprintf(fout_3pts,"\n");
  
  contr_3pts_time+=take_time();
}

//compute all the combinations of three point correlators
void calculate_all_3pts()
{
  for(int mp=-1;mp<=1;mp++)
    for(int imass=0;imass<nheavy;imass++)
      {
	master_printf("computing 3pts contractions %d %d\n",mp,imass);
	calculate_all_3pts_prop_combo(S0H[iS0H(0,imass)],S1LH[iS1(mp,imass)],mass[1+imass],mp*theta[imass],"light_spec");
	calculate_all_3pts_prop_combo(S0L[iS0L(0,0)],S1HL[iS1(mp,imass)],mass[1+imass],mp*theta[imass],"charm_spec");
      }
}

//check if the time is enough
int check_remaining_time()
{
  int enough_time;

  //check remaining time
  double temp_time=take_time()+tot_prog_time;
  double ave_time=temp_time/nanalized_conf;
  double left_time=wall_time-temp_time;
  enough_time=left_time>(ave_time*1.1);

  master_printf("Remaining time: %lg sec\n",left_time);
  master_printf("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
  if(enough_time) master_printf("Continuing with next conf!\n");
  else master_printf("Not enough time, exiting!\n");
  
  return enough_time;
}

//Finalization
void close_semileptonic()
{
  tot_prog_time+=take_time();

  contr_time=contr_2pts_time+contr_3pts_time+contr_save_time;
  
  master_printf("\n");
  master_printf("Total time: %g, of which:\n",tot_prog_time);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf(" - %02.2f%s to smear %d propagators (%2.2gs avg)\n",smear_time/tot_prog_time*100,"%",nsmear_tot,smear_time/nsmear_tot);
  master_printf(" - %02.2f%s to perform %d contr. (%2.2gs avg) of which:\n",contr_time/tot_prog_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
  master_printf("   * %02.2f%s to compute two points\n",contr_2pts_time*100.0/contr_time,"%");
  master_printf("   * %02.2f%s to compute three points\n",contr_3pts_time*100.0/contr_time,"%");
  master_printf("   * %02.2f%s to save correlations\n",contr_save_time*100.0/contr_time,"%");
  
  nissa_free(conf);nissa_free(sme_conf);
  for(int imass=0;imass<3*nheavy;imass++)
    {
      nissa_free(S1LH[imass]);
      nissa_free(S1HL[imass]);
    }
  for(int imass=0;imass<2*nheavy+1;imass++) nissa_free(S0L[imass]);
  for(int imass=0;imass<3*nheavy;imass++) nissa_free(S0H[imass]);

  nissa_free(S1LH);
  nissa_free(S1HL);
  nissa_free(S0L);
  nissa_free(S0H);
  nissa_free(sequential_source);
  nissa_free(contr_2pts);nissa_free(contr_3pts);
  nissa_free(op1_2pts);nissa_free(op2_2pts);
  nissa_free(op1_3pts);nissa_free(op2_3pts);
  nissa_free(temp_sol);
  nissa_free(source);nissa_free(original_source);
}

void in_main(int narg,char **arg)
{
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  initialize_semileptonic(arg[1]);

  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf && enough_time && read_conf_parameters(&iconf))
    {
      setup_conf();
      generate_source();

      calculate_all_S0();
      calculate_all_2pts(fout_2pts_sl);
      calculate_all_S1();
      calculate_all_3pts();
      
      /////////////// standing smeared two points ////////////
      
      //smear all S0 propagators to compute ss corr
      smear_all_S0();
      
      //compute smeared smeared corelators
      calculate_all_2pts(fout_2pts_ss);
      
      ////////////////////////////////////////////////////////
      
      close_file(fout_2pts_sl);
      close_file(fout_2pts_ss);
      close_file(fout_3pts);
      
      nanalized_conf++;      
      enough_time=check_remaining_time();
    }
  
  if(iconf==ngauge_conf) master_printf("Finished all the confs!\n");
  
  close_semileptonic();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
    
  return 0;
}
