#include <nissa.hpp>

#include "nissa.hpp"
#include "driver_corr.hpp"

char outfolder[100],conf_path[100];
int ngauge_conf,nanalyzed_conf;
int T,L,seed;
int nsources,twall,tsep;
int ninv_tot=0;
double wall_time,conf_smear_time=0,tot_prog_time=0,inv_time=0,smear_time=0,corr_time=0;
double theta_needed,mass,residue,kappa;
ape_pars_t ape_smearing_pars;
spincolor *source,*cgm_solution[1],*temp_vec[2];
quad_su3 *conf,*sme_conf;
colorspinspin *loc_source,*sme_source,*seq_source,*C_loc_source,*C_sme_source,*S1;
colorspinspin *S0_std_lcsi_smso,*S0_mov_lcsi_smso;
colorspinspin *S0_std_smsi_smso,*S0_mov_smsi_smso;
colorspinspin *S0_std_C_lcsi_smso,*S0_std_C_smsi_smso;
colorspinspin *S0_std_C_lcsi_lcso,*S0_std_C_smsi_lcso;
colorspinspin *S0_std_lcsi_C_lcso,*S0_std_smsi_C_lcso;
colorspinspin *S0_std_lcsi_C_smso,*S0_std_smsi_C_smso;
colorspinspin *S0_std_lcsi_lcso,*S0_std_smsi_lcso;

double *loc_corr,*glb_corr;
int ncorr_tot=0;
momentum_t old_theta,put_theta;
int rspec=0;

int gaussian_niter;
double gaussian_kappa;

const int i2pts_VsVs_st=0,i2pts_VsVs_mv=1;
const int i2pts_ClCl=2,i2pts_PsCl=3,i2pts_CsCl=4;
const int i2pts_ClPs=5,i2pts_PsPs=6,i2pts_CsPs=7;
const int i2pts_ClCs=8,i2pts_PsCs=9,i2pts_CsCs=10;
const int i3pts_base=11;
const int nop_tot=i3pts_base+3;

//apply the C operator
void apply_C_operator(colorspinspin *out,colorspinspin *in)
{
  colorspinspin *temp_add=nissa_malloc("temp_add",loc_vol+bord_vol,colorspinspin);
  colorspinspin *temp_res=nissa_malloc("temp_res",loc_vol+bord_vol,colorspinspin);
  vector_reset(temp_res);
  
  //add the three contributions
  for(int mu=0;mu<3;mu++)
    {
      int idir=1+mu,iga=16+mu;
      apply_nabla_i(temp_add,in,sme_conf,idir);
      NISSA_LOC_VOL_LOOP(ivol)
	for(int ic=0;ic<3;ic++)
	  {
	    spinspin t;
	    unsafe_dirac_prod_spinspin(t,base_gamma+iga,temp_add[ivol][ic]);
	    spinspin_summassign(temp_res[ivol][ic],t);
	  }
    }
  
  vector_copy(out,temp_res);
  
  nissa_free(temp_add);
  nissa_free(temp_res);

  set_borders_invalid(out);
}

void initialize_semileptonic(char *input_path)
{
  open_input(input_path);

  // 1) Read information about the gauge conf
  
  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Wall time
  read_str_double("WallTime",&wall_time);
  //Kappa is read really only for tm
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the source
  
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  read_str_int("NSources",&nsources);
  
  // 3) Read smearing info
  
  read_str_double("GaussianKappa",&gaussian_kappa);
  read_str_int("GaussianNiter",&gaussian_niter);
  read_ape_pars(ape_smearing_pars);
  
  // 4) Read masses, residue and list of confs

  read_str_double("Mass",&mass);
  read_str_double("Residue",&residue);
  read_str_double("ThetaNeeded",&theta_needed);
  read_str_int("TSep",&tsep);
  read_str_int("NGaugeConf",&ngauge_conf);
  
  ////////////////////////////////// allocate //////////////////////
  
  loc_corr=nissa_malloc("loc_corr",nop_tot*glb_size[0],double);
  glb_corr=nissa_malloc("glb_corr",nop_tot*glb_size[0],double);
  
  //temp vecs
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol+bord_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol+bord_vol,spincolor);
  cgm_solution[0]=nissa_malloc("cgm_solution",loc_vol+bord_vol,spincolor);

  //allocate gauge conf, all the needed spincolor and propagators
  conf=nissa_malloc("or_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  sme_conf=nissa_malloc("sm_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //Allocate S0
  S0_std_lcsi_lcso=nissa_malloc("S0_std_lcsi_lcso",loc_vol+bord_vol,colorspinspin);
  S0_std_smsi_lcso=nissa_malloc("S0_std_smsi_lcso",loc_vol+bord_vol,colorspinspin);
  S0_std_lcsi_smso=nissa_malloc("S0_std_lcsi_smso",loc_vol+bord_vol,colorspinspin);
  S0_std_smsi_smso=nissa_malloc("S0_std_smsi_smso",loc_vol+bord_vol,colorspinspin);
  S0_std_C_lcsi_smso=nissa_malloc("S0_std_C_lcsi_smso",loc_vol+bord_vol,colorspinspin);
  S0_std_C_smsi_smso=nissa_malloc("S0_std_C_smsi_smso",loc_vol+bord_vol,colorspinspin);
  S0_std_C_lcsi_lcso=nissa_malloc("S0_std_C_lcsi_lcso",loc_vol+bord_vol,colorspinspin);
  S0_std_C_smsi_lcso=nissa_malloc("S0_std_C_smsi_lcso",loc_vol+bord_vol,colorspinspin);
  S0_std_lcsi_C_lcso=nissa_malloc("S0_std_lcsi_C_lcso",loc_vol+bord_vol,colorspinspin);
  S0_std_smsi_C_lcso=nissa_malloc("S0_std_smsi_C_lcso",loc_vol+bord_vol,colorspinspin);
  S0_std_lcsi_C_smso=nissa_malloc("S0_std_lcsi_C_smso",loc_vol+bord_vol,colorspinspin);
  S0_std_smsi_C_smso=nissa_malloc("S0_std_smsi_C_smso",loc_vol+bord_vol,colorspinspin);
  S0_mov_lcsi_smso=nissa_malloc("S0_mov_lcsi_smso",loc_vol+bord_vol,colorspinspin);
  S0_mov_smsi_smso=nissa_malloc("S0_mov_smsi_smso",loc_vol+bord_vol,colorspinspin);
  S1=nissa_malloc("S1",loc_vol+bord_vol,colorspinspin);
  
  //Allocate source
  loc_source=nissa_malloc("loc_source",loc_vol+bord_vol,colorspinspin);
  sme_source=nissa_malloc("sme_source",loc_vol+bord_vol,colorspinspin);
  seq_source=nissa_malloc("seq_source",loc_vol+bord_vol,colorspinspin);
  C_loc_source=nissa_malloc("C_loc_source",loc_vol+bord_vol,colorspinspin);
  C_sme_source=nissa_malloc("C_sme_source",loc_vol+bord_vol,colorspinspin);
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
}

//find a new conf
int read_conf_parameters(int *iconf)
{
  int ok_conf=0;

  do
    {
      //Gauge path
      read_str(conf_path,1024);
      
      //Out folder
      read_str(outfolder,1024);
      
      //Check if the conf has been finished or is already running
      MASTER_PRINTF("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);

      //if not chosen
      if(!dir_exists(outfolder))
	{
          MASTER_PRINTF(" Configuration \"%s\" not yet analyzed, starting\n",conf_path);
	  if(create_dir(outfolder)) CRASH(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
	  ok_conf=1;
	}
      (*iconf)++;
    }
  while(!ok_conf && (*iconf)<ngauge_conf);
  
  MASTER_PRINTF("\n");
  
  return ok_conf;
}

//read the conf and setup it
void setup_conf()
{
  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_ildg_gauge_conf(conf,conf_path);
  MASTER_PRINTF("plaq: %.16g\n",global_plaquette_lx_conf(conf));

  conf_smear_time-=take_time();
  ape_spatial_smear_conf(sme_conf,conf,ape_smearing_pars.alpha,ape_smearing_pars.nlevels);
  conf_smear_time+=take_time();
  
  MASTER_PRINTF("smerded plaq: %.16g\n",global_plaquette_lx_conf(sme_conf));
  
  //put the anti-periodic condition on the temporal border
  old_theta[0]=old_theta[1]=old_theta[2]=old_theta[3]=0;
  put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
  adapt_theta(conf,old_theta,put_theta,1,1);
  
  //reset the corr
  vector_reset(loc_corr);
}

//Finalization
void close_gevp_3pts()
{
  MASTER_PRINTF("\n");
  MASTER_PRINTF("Inverted %d configurations.\n",nanalyzed_conf);
  MASTER_PRINTF("Total time: %g s, of which:\n",tot_prog_time);
  MASTER_PRINTF(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",
		ninv_tot,inv_time/ninv_tot);
  MASTER_PRINTF("  of which  %02.2f%s for %d cgm inversion overhead (%2.2gs avg)\n",cgm_inv_over_time/inv_time*100,"%",
                ninv_tot,cgm_inv_over_time/ninv_tot);
  MASTER_PRINTF(" - %02.2f%s to smear configuration\n",conf_smear_time*100.0/tot_prog_time,"%");
  MASTER_PRINTF(" - %02.2f%s to sink-smear propagators\n",smear_time*100.0/tot_prog_time,"%");
  MASTER_PRINTF(" - %02.2f%s to compute %d correlations (%2.2gs avg)\n",corr_time/tot_prog_time*100,"%",
		ncorr_tot,corr_time/ncorr_tot);
  
  nissa_free(conf);nissa_free(sme_conf);
  nissa_free(S0_std_lcsi_lcso);
  nissa_free(S0_std_lcsi_smso);
  nissa_free(S0_std_smsi_lcso);
  nissa_free(S0_std_smsi_smso);
  nissa_free(S0_std_C_lcsi_smso);
  nissa_free(S0_std_C_smsi_smso);
  nissa_free(S0_std_C_lcsi_lcso);
  nissa_free(S0_std_C_smsi_lcso);
  nissa_free(S0_std_lcsi_C_lcso);
  nissa_free(S0_std_smsi_C_lcso);
  nissa_free(S0_std_lcsi_C_smso);
  nissa_free(S0_std_smsi_C_smso);
  nissa_free(S0_mov_lcsi_smso);
  nissa_free(S0_mov_smsi_smso);
  nissa_free(S1);
  nissa_free(loc_source);
  nissa_free(sme_source);
  nissa_free(C_loc_source);
  nissa_free(C_sme_source);
  nissa_free(seq_source);
  nissa_free(source);
  nissa_free(glb_corr);
  nissa_free(loc_corr);
  nissa_free(cgm_solution[0]);
  nissa_free(temp_vec[0]);
  nissa_free(temp_vec[1]);

  close_nissa();
}

//generate the source
void generate_source(int isource)
{
  twall=(int)rnd_get_unif(&glb_rnd_gen,0,glb_size[0]-1);
  MASTER_PRINTF("Source %d on t=%d\n",isource,twall);
  generate_spindiluted_source(loc_source,rnd_type_map[4],twall);

  //apply C
  apply_C_operator(C_loc_source,loc_source);

  //smear the source
  MASTER_PRINTF("Source Smearing\n");
  smear_time-=take_time();
  gaussian_smearing(sme_source,loc_source,sme_conf,gaussian_kappa,gaussian_niter);
  gaussian_smearing(C_sme_source,C_loc_source,sme_conf,gaussian_kappa,gaussian_niter);
  smear_time+=take_time();
}

//calculate S0_lcsi propagator
void calculate_S0_lcsi(colorspinspin *out,double theta,colorspinspin *in)
{
  put_theta[1]=put_theta[2]=put_theta[3]=theta; //if multiple sources
  adapt_theta(conf,old_theta,put_theta,1,1);
  
  //loop over source dirac index
  for(int id=0;id<4;id++)
    { 
      //add gamma5
      get_spincolor_from_colorspinspin(source,in,id);
      safe_dirac_prod_spincolor(source,base_gamma+5,source);
      
      //invert
      const int niter_max=10000;
      double part_time=-take_time();
      inv_tmQ2_cgm(cgm_solution,conf,kappa,&mass,1,niter_max,&residue,source);
      part_time+=take_time();
      inv_time+=part_time;
      ninv_tot++;
      
      MASTER_PRINTF("Finished the inversion of dirac index %d in %g sec\n",id,part_time);
      
      //reconstruct the doublet
      reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass,cgm_solution[0]);
      put_spincolor_into_colorspinspin(out,temp_vec[rspec],id);
    }

  //rotate to physical basis, remember that D^-1 rotate opposite than D!
  rotate_vol_colorspinspin_to_physical_basis(out,!rspec,!rspec);
  MASTER_PRINTF("Propagators rotated\n");
}

//smear the sink
void calculate_S0_smsi()
{
  MASTER_PRINTF("S0 smearing\n");
  smear_time-=take_time();
  gaussian_smearing(S0_std_smsi_lcso,S0_std_lcsi_lcso,sme_conf,gaussian_kappa,gaussian_niter);
  gaussian_smearing(S0_std_smsi_smso,S0_std_lcsi_smso,sme_conf,gaussian_kappa,gaussian_niter);
  gaussian_smearing(S0_std_smsi_C_lcso,S0_std_lcsi_C_lcso,sme_conf,gaussian_kappa,gaussian_niter);
  gaussian_smearing(S0_mov_smsi_smso,S0_mov_lcsi_smso,sme_conf,gaussian_kappa,gaussian_niter);
  gaussian_smearing(S0_std_smsi_C_smso,S0_std_lcsi_C_smso,sme_conf,gaussian_kappa,gaussian_niter);
  smear_time+=take_time();
}

//generate the sequential source
void generate_sequential_source(int iseq)
{
  MASTER_PRINTF("\nCreating the sequential source for sequential operator %d\n",iseq);
  
  //seq source smearing
  switch(iseq)
    {
    case 0:
      apply_C_operator(seq_source,S0_std_lcsi_smso);
      break;
    case 1:
      safe_dirac_prod_colorspinspin(seq_source,base_gamma+5,S0_std_smsi_smso);
      break;
    case 2:
      apply_C_operator(seq_source,S0_std_smsi_smso);
      break;
    default:
      CRASH("unknown operator");
    }

  //seq source smearing
  if(iseq==1||iseq==2)
    {
      MASTER_PRINTF("Seq source smearing\n");
      smear_time-=take_time();
      gaussian_smearing(seq_source,seq_source,sme_conf,gaussian_kappa,gaussian_niter);
      smear_time+=take_time();
    }
  
  //rotate as r because it's D^-1
  rotate_vol_colorspinspin_to_physical_basis(seq_source,rspec,rspec);
  select_propagator_timeslice(seq_source,seq_source,(tsep+twall)%glb_size[0]);
  
  MASTER_PRINTF("Sequential source created\n\n");
}
  
//calculate S1 propagator
void calculate_S1(int iseq)
{
  //put the theta
  put_theta[1]=put_theta[2]=put_theta[3]=0;
  adapt_theta(conf,old_theta,put_theta,1,1);
  
  generate_sequential_source(iseq);
  
  //loop over source dirac index
  for(int id=0;id<4;id++)
    { 
      //add gamma5
      get_spincolor_from_colorspinspin(source,seq_source,id);
      safe_dirac_prod_spincolor(source,base_gamma+5,source);
      
      //invert
      const int niter_max=10000;
      double part_time=-take_time();
      inv_tmQ2_cgm(cgm_solution,conf,kappa,&mass,1,niter_max,&residue,source);
      part_time+=take_time();
      inv_time+=part_time;
      ninv_tot++;
      
      MASTER_PRINTF("Finished the inversion of dirac index %d in %g sec\n",id,part_time);
      
      //reconstruct the doublet: r(S1)=!r(spec), so we have to multiply by Q+ if r(spec)==1 and Q- if 0
      double reco_mass=-mass;
      apply_tmQ(temp_vec[0],conf,kappa,reco_mass,cgm_solution[0]);
      put_spincolor_into_colorspinspin(S1,temp_vec[0],id);
    }
  
  //put the (1+-ig5)/sqrt(2) factor if tm. On the source rotate as r_spec, on the sink as !r_spec
  rotate_vol_colorspinspin_to_physical_basis(S1,!rspec,rspec);
  
  MASTER_PRINTF("S1 propagators %d computed\n",iseq);
}

//calculate all 3pts
void calculate_3pts(int iop_seq)
{
  double *loc=loc_corr+glb_size[0]*(i3pts_base+iop_seq);
  
  //loop on the three polarizations
  if(1)
  for(int ipol=0;ipol<3;ipol++)
    {
      //find the two possible insertions
      int ins1=(ipol+1)%3,ins2=(ins1+1)%3;
      int ig_pol=5+1+ipol,ig1=5+1+ins1,ig2=5+1+ins2;
      
      //loop over all the volume
      NISSA_LOC_VOL_LOOP(ivol)
	{
	  int glb_t=glb_coord_of_loclx[ivol][0],dst_t=(glb_size[0]+glb_t-twall)%glb_size[0];
	  
	  //take the two contractions
	  complex ctemp1,ctemp2;
	  trace_g_css_dag_g_css(ctemp1,base_gamma+ig_pol,S0_mov_lcsi_smso[ivol],base_gamma+ig1,S1[ivol]);
	  trace_g_css_dag_g_css(ctemp2,base_gamma+ig_pol,S0_mov_lcsi_smso[ivol],base_gamma+ig2,S1[ivol]);
	  loc[dst_t]+=-(ctemp1[IM]-ctemp2[IM])/3;
	}
    }
  
  //checks
  if(0)
    {
      vector_reset(loc_corr);
      
      NISSA_LOC_VOL_LOOP(ivol)
      {
	int glb_t=glb_coord_of_loclx[ivol][0],dst_t=(glb_size[0]+glb_t-twall)%glb_size[0];
	
	//take the two contractions
	complex ctemp={0,0};
	if(1) trace_g_css_dag_g_css(ctemp,base_gamma+6,S0_std_lcsi_smso[ivol],base_gamma+7,S1[ivol]); //V1V2
	loc[dst_t]+=ctemp[IM];
	if(0) trace_g_css_dag_g_css(ctemp,base_gamma,S0_std_lcsi_smso[ivol],base_gamma+9,S1[ivol]); //V0P5
	if(0) trace_g_css_dag_g_css(ctemp,base_gamma,S0_std_lcsi_smso[ivol],base_gamma,
				    S0_std_lcsi_smso[ivol]); //2pts
	if(0) trace_g_css_dag_g_css(ctemp,base_gamma+5,sme_source[ivol],base_gamma+0,S1[ivol]); //2pts check
	//loc[dst_t]+=ctemp[RE];
      }
    }
}

//calculate all 2pts
void calculate_2pts()
{
  //C operator on the sink must be applied after smearing
  apply_C_operator(S0_std_C_lcsi_smso,S0_std_lcsi_smso);
  apply_C_operator(S0_std_C_smsi_smso,S0_std_smsi_smso);
  apply_C_operator(S0_std_C_lcsi_lcso,S0_std_lcsi_lcso);
  apply_C_operator(S0_std_C_smsi_lcso,S0_std_smsi_lcso);
  
  //loop on the three polarizations
  NISSA_LOC_VOL_LOOP(ivol)
  {
    int glb_t=glb_coord_of_loclx[ivol][0],dst_t=(glb_size[0]+glb_t-twall)%glb_size[0];
    
    //take the three vector polarization
    complex temp;
    for(int mu=0;mu<3;mu++)
      {
	//VV(st)
	trace_g_css_dag_g_css(temp,base_gamma+6+mu,S0_std_smsi_smso[ivol],base_gamma+6+mu,S0_std_smsi_smso[ivol]);
	loc_corr[i2pts_VsVs_st*glb_size[0]+dst_t]+=-temp[RE]/3;
	//VV(mv)
	trace_g_css_dag_g_css(temp,base_gamma+6+mu,S0_std_smsi_smso[ivol],base_gamma+6+mu,S0_mov_smsi_smso[ivol]);
	loc_corr[i2pts_VsVs_mv*glb_size[0]+dst_t]+=-temp[RE]/3;
      }
    
    //ClCl
    trace_g_css_dag_g_css(temp,base_gamma+5,S0_std_lcsi_C_lcso[ivol],base_gamma+5,S0_std_C_lcsi_lcso[ivol]);
    loc_corr[i2pts_ClCl*glb_size[0]+dst_t]+=temp[RE];
    //PsCl
    trace_g_css_dag_g_css(temp,base_gamma+5,S0_std_smsi_C_lcso[ivol],base_gamma+0,S0_std_smsi_lcso[ivol]);
    loc_corr[i2pts_PsCl*glb_size[0]+dst_t]+=temp[RE];
    //CsCl
    trace_g_css_dag_g_css(temp,base_gamma+5,S0_std_smsi_C_lcso[ivol],base_gamma+5,S0_std_C_smsi_lcso[ivol]);
    loc_corr[i2pts_CsCl*glb_size[0]+dst_t]+=temp[RE];
    
    //ClPs
    trace_g_css_dag_g_css(temp,base_gamma+0,S0_std_lcsi_smso[ivol],base_gamma+5,S0_std_C_lcsi_smso[ivol]);
    loc_corr[i2pts_ClPs*glb_size[0]+dst_t]+=temp[RE];
    //PsPs
    trace_g_css_dag_g_css(temp,base_gamma+0,S0_std_smsi_smso[ivol],base_gamma+0,S0_std_smsi_smso[ivol]);
    loc_corr[i2pts_PsPs*glb_size[0]+dst_t]+=temp[RE];
    //CsPs
    trace_g_css_dag_g_css(temp,base_gamma+0,S0_std_smsi_smso[ivol],base_gamma+5,S0_std_C_smsi_smso[ivol]);
    loc_corr[i2pts_CsPs*glb_size[0]+dst_t]+=temp[RE];

    //ClCs
    trace_g_css_dag_g_css(temp,base_gamma+5,S0_std_lcsi_C_smso[ivol],base_gamma+5,S0_std_C_lcsi_smso[ivol]);
    loc_corr[i2pts_ClCs*glb_size[0]+dst_t]+=temp[RE];
    //PsCs
    trace_g_css_dag_g_css(temp,base_gamma+5,S0_std_smsi_C_smso[ivol],base_gamma+0,S0_std_smsi_smso[ivol]);
    loc_corr[i2pts_PsCs*glb_size[0]+dst_t]+=temp[RE];
    //CsCs
    trace_g_css_dag_g_css(temp,base_gamma+5,S0_std_smsi_C_smso[ivol],base_gamma+5,S0_std_C_smsi_smso[ivol]);
    loc_corr[i2pts_CsCs*glb_size[0]+dst_t]+=temp[RE];
  }
}

//check if the time is enough
int check_remaining_time()
{
  int enough_time;

  //check remaining time
  double temp_time=take_time()+tot_prog_time;
  double ave_time=temp_time/nanalyzed_conf;
  double left_time=wall_time-temp_time;
  enough_time=left_time>(ave_time*1.1);

  MASTER_PRINTF("Remaining time: %lg sec\n",left_time);
  MASTER_PRINTF("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
  if(enough_time) MASTER_PRINTF("Continuing with next conf!\n");
  else MASTER_PRINTF("Not enough time, exiting!\n");
  
  return enough_time;
}

//reduce and write
void write_corr()
{
  //open output file
  char path[1024];
  sprintf(path,"%s/corr",outfolder);
  FILE *fout=open_text_file_for_output(path);
  
  //reduce
  MPI_Reduce(loc_corr,glb_corr,glb_size[0]*nop_tot,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  for(int iop=0;iop<nop_tot;iop++)
    {
      const char tag[nop_tot][20]={"Vs_Vs(std)","Vs_Vs(theta)",
				   "Cl_Cl","Ps_Cl","Cs_Cl",
				   "Cl_Ps","Ps_Ps","Cs_Ps",
				   "Cl_Cs","Ps_Cs","Cs_Cs",
				   "Cl_VJ_VK","Ps_VJ_VK","Cs_VJ_VK"};
      master_fprintf(fout," # OPERATOR %s\n",tag[iop]);
      for(int t=0;t<glb_size[0];t++) master_fprintf(fout,"%+16.16lg\n",glb_corr[t+iop*glb_size[0]]/nsources);
      
      master_fprintf(fout,"\n");
    }
  
  close_file(fout);
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //initialize the program
  if(narg<2) CRASH("Use: %s input_file",arg[0]);
  initialize_semileptonic(arg[1]);
  
  //loop over the configs
  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf && enough_time && read_conf_parameters(&iconf))
    {
      //smear the conf and generate the source
      setup_conf();
      
      //average all the required sources
      for(int isource=0;isource<nsources;isource++)
	{
	  generate_source(isource);
	  
	  //compute S0_lcsi propagator
	  calculate_S0_lcsi(S0_std_lcsi_lcso,0,loc_source);
	  calculate_S0_lcsi(S0_std_lcsi_smso,0,sme_source);
	  calculate_S0_lcsi(S0_std_lcsi_C_lcso,0,C_loc_source);
	  calculate_S0_lcsi(S0_std_lcsi_C_smso,0,C_sme_source);
	  calculate_S0_lcsi(S0_mov_lcsi_smso,theta_needed,sme_source);
	  calculate_S0_smsi();
	  calculate_2pts();
	  
	  //compute all operators
	  for(int iop_seq=0;iop_seq<3;iop_seq++)
	    {
	      calculate_S1(iop_seq);
	      calculate_3pts(iop_seq);
	    }
	}
      
      //reduce and write three points      
      write_corr();
      
      nanalyzed_conf++;
      enough_time=check_remaining_time();
    }
  
  if(iconf==ngauge_conf) MASTER_PRINTF("Finished all the conf!\n");
  
  tot_prog_time+=take_time();
  close_gevp_3pts();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  
  return 0;
}
