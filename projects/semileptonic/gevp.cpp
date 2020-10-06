#include <nissa.hpp>

#include "nissa.hpp"
#include "driver_corr.hpp"

char outfolder[100],conf_path[100];
int ngauge_conf,nanalyzed_conf;
int nsm_lev,*gaussian_niter;
int T,L,seed;
int ncorr_computed_tot,nsources;
int ndoubles,twall;
int ninv_tot=0;
double wall_time,conf_smear_time=0,tot_prog_time=0,inv_time=0,smear_time=0,corr_time=0;
double *glb_corr,*loc_corr;
double mass,residue,kappa,gaussian_kappa;
two_pts_comp_t two_pts_comp[4][4][4][4];
int two_pts_comp_off[4][4][4][4];
ape_pars_t ape_smearing_pars;
spincolor *source,*temp_vec[2],*cgm_solution[1];
quad_su3 *conf,*sme_conf;
colorspinspin *ori_source,*temp_transp,*S0[2][4][4];
int ncorr_tot=0;

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
  
  // 3) Smearing parameters
  
  //Smearing parameters
  read_str_double("GaussianKappa",&gaussian_kappa);
  read_list_of_ints("GaussianNiter",&nsm_lev,&gaussian_niter);
  for(int iter=1;iter<nsm_lev;iter++)
    if(gaussian_niter[iter]<gaussian_niter[iter-1])
      crash("Error, gaussian lev %d minor than %d (%d, %d)!\n",iter,iter-1,gaussian_niter[iter],gaussian_niter[iter-1]);
  read_ape_pars(ape_smearing_pars);
  
  // 4) Correlations
  
  int nop[4][4],*op_list[4][4];
  for(int mu=0;mu<4;mu++)
    for(int nu=mu;nu<4;nu++)
      {
	char tag[100];
	sprintf(tag,"OPList%d%d",mu,nu);
	read_list_of_ints(tag,&(nop[mu][nu]),&(op_list[mu][nu]));
      }
  
  // 5) Read masses, residue and list of confs

  read_str_double("Mass",&mass);
  read_str_double("Residue",&residue);
  read_str_int("NGaugeConf",&ngauge_conf);
  
  ////////////////////////////////////////////////////////////////////////////
  
  //insert all corr possibility
  int ncontr_tot=0;
  for(int mu_so=0;mu_so<4;mu_so++)
    for(int nu_so=0;nu_so<4;nu_so++)
      for(int mu_si=0;mu_si<4;mu_si++)
	for(int nu_si=0;nu_si<4;nu_si++)
	  {
	    //access to minmax list
	    int so1=std::min(mu_so,nu_so),so2=std::max(mu_so,nu_so);
	    int si1=std::min(mu_si,nu_si),si2=std::max(mu_si,nu_si);
	    
	    //fill it
	    uint16_t icorr=0;
	    for(int iop_so=0;iop_so<nop[so1][so2];iop_so++)
	      for(int iop_si=0;iop_si<nop[si1][si2];iop_si++)
		{
		  int ga_so=op_list[so1][so2][iop_so];
		  int ga_si=op_list[si1][si2][iop_si];
		  two_pts_comp[mu_so][nu_so][mu_si][nu_si].corr_name[icorr]+=gtag[ga_si]+(std::string)gtag[ga_so];
		  two_pts_comp[mu_so][nu_so][mu_si][nu_si].add_sink_source_corr(icorr++,1,0,ga_si,ga_so);
		}		  

	    //take the offset and increment it
	    two_pts_comp[mu_so][nu_so][mu_si][nu_si].ncorr=icorr;
	    two_pts_comp_off[mu_so][nu_so][mu_si][nu_si]=ncorr_tot;
	    ncorr_tot+=icorr;
	    ncontr_tot+=two_pts_comp[mu_so][nu_so][mu_si][nu_si].size();
	  }
  verbosity_lv2_master_printf("Read %d corrs, corresponding to %d contr\n",ncorr_tot,ncontr_tot);
  ndoubles=nsm_lev*nsm_lev*glb_size[0]*ncorr_tot;
  
  //correlations
  glb_corr=nissa_malloc("glb_corr",ndoubles,double);
  loc_corr=nissa_malloc("loc_corr",ndoubles,double);

  //temp vecs
  temp_transp=nissa_malloc("temp_transp",loc_vol+bord_vol,colorspinspin);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol+bord_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol+bord_vol,spincolor);
  cgm_solution[0]=nissa_malloc("cgm_solution",loc_vol+bord_vol,spincolor);

  //allocate gauge conf, all the needed spincolor and propagators
  conf=nissa_malloc("or_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  sme_conf=nissa_malloc("sm_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //Allocate all the S0 PROP_TYPE vectors
  for(int r=0;r<2;r++)
    for(int imu=0;imu<4;imu++)
      for(int inu=0;inu<4;inu++)
	S0[r][imu][inu]=nissa_malloc("S0[r]",loc_vol+bord_vol,colorspinspin);
  
  //Allocate one spincolor for the source
  ori_source=nissa_malloc("ori_source",loc_vol+bord_vol,colorspinspin);
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
      master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);

      //if not chosen
      if(!dir_exists(outfolder))
	{
          master_printf(" Configuration \"%s\" not yet analyzed, starting\n",conf_path);
	  if(create_dir(outfolder)) crash(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
	  ok_conf=1;
	}
      (*iconf)++;
    }
  while(!ok_conf && (*iconf)<ngauge_conf);
  
  master_printf("\n");
  
  return ok_conf;
}

//read the conf and setup it
void setup_conf()
{
  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_ildg_gauge_conf(conf,conf_path);
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));

  conf_smear_time-=take_time();
  ape_spatial_smear_conf(sme_conf,conf,ape_smearing_pars.alpha,ape_smearing_pars.nlevels);
  conf_smear_time+=take_time();
  
  master_printf("smeared plaq: %.16g\n",global_plaquette_lx_conf(sme_conf));
  
  //put the anti-periodic condition on the temporal border
  momentum_t old_theta,put_theta;
  old_theta[0]=0;
  put_theta[0]=1;
  for(int mu=1;mu<NDIM;mu++) old_theta[mu]=put_theta[mu]=0;
  adapt_theta(conf,old_theta,put_theta,1,1);
  
  //reset the corr
  vector_reset(loc_corr);
}

//Finalization
void close_semileptonic()
{
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g s, of which:\n",tot_prog_time);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",
		ninv_tot,inv_time/ninv_tot);
  master_printf("  of which  %02.2f%s for %d cgm inversion overhead (%2.2gs avg)\n",cgm_inv_over_time/inv_time*100,"%",
                ninv_tot,cgm_inv_over_time/ninv_tot);
  master_printf(" - %02.2f%s to smear configuration\n",conf_smear_time*100.0/tot_prog_time,"%");
  master_printf(" - %02.2f%s to sink-smear propagators\n",smear_time*100.0/tot_prog_time,"%");
  master_printf(" - %02.2f%s to compute %d correlations (%2.2gs avg)\n",corr_time/tot_prog_time*100,"%",
		ncorr_computed_tot,corr_time/ncorr_computed_tot);
  
  nissa_free(conf);nissa_free(sme_conf);
  for(int r=0;r<2;r++) for(int imu=0;imu<4;imu++) for(int inu=0;inu<4;inu++)nissa_free(S0[r][imu][inu]);
  nissa_free(source);
  nissa_free(ori_source);
  nissa_free(loc_corr);
  nissa_free(glb_corr);
  nissa_free(cgm_solution[0]);
  nissa_free(temp_vec[0]);
  nissa_free(temp_vec[1]);
  nissa_free(temp_transp);
}

//calculate the standard propagators
void calculate_S0(int ism_lev_so)
{
  //smear additively the source
  master_printf("\nSource Smearing level: %d\n",ism_lev_so);
  int nsme=gaussian_niter[ism_lev_so];
  if(ism_lev_so>0) nsme-=gaussian_niter[ism_lev_so-1];
  gaussian_smearing(ori_source,ori_source,sme_conf,gaussian_kappa,nsme);
  master_printf("\n");
  
  //loop over derivative of the source
  for(int id=0;id<4;id++)
    for(int imu=0;imu<4;imu++)
      { 
	//add gamma5
	get_spincolor_from_colorspinspin(source,ori_source,id);
	safe_dirac_prod_spincolor(source,base_gamma+5,source);
      
	if(imu>0) apply_nabla_i(source,source,sme_conf,imu);
	
	//invert
	const int niter_max=10000;
	double part_time=-take_time();
	inv_tmQ2_cgm(cgm_solution,conf,kappa,&mass,1,niter_max,&residue,source);
	part_time+=take_time();
	inv_time+=part_time;
	ninv_tot++;
	
	master_printf("Finished the inversion of ");
	master_printf("source displaced in %d ",imu);
	master_printf("dirac index %d in %g sec\n",id,part_time);
	
	//reconstruct the doublet
	reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass,cgm_solution[0]);
	for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
	  put_spincolor_into_colorspinspin(S0[r][imu][0],temp_vec[r],id);
      }

  //rotate to physical basis
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int imu=0;imu<4;imu++)
      rotate_vol_colorspinspin_to_physical_basis(S0[r][imu][0],!r,!r);
  master_printf("Propagators rotated\n");
  
  //create sink-derived propagators
  for(int r=0;r<2;r++)
    for(int imu=0;imu<4;imu++)
      for(int inu=1;inu<4;inu++)
	apply_nabla_i(S0[r][imu][inu],S0[r][imu][0],sme_conf,inu);
  
  master_printf("\n");
}

//transpose dirac indices with spatial ones, so that destination
//is arranged as: t,id_si,id_so,ri,spat,color
THREADABLE_FUNCTION_2ARG(prepare_prop_for_new_contraction, colorspinspin*,prop, colorspinspin*,aux)
{
  GET_THREAD_ID();
  
  int spat_color_vol=3*loc_vol/loc_size[0];
  
  //fill aux
  for(int t=0;t<loc_size[0];t++)
    NISSA_PARALLEL_LOOP(ispat_color,0,spat_color_vol)
      for(int idirac_ri=0;idirac_ri<32;idirac_ri++)
        ((double*)aux)[ispat_color+spat_color_vol*(idirac_ri+32*t)]=
          ((double*)prop)[idirac_ri+32*(ispat_color+spat_color_vol*t)];
  THREAD_BARRIER();
  
  //copy aux to prop
  parallel_memcpy(prop,aux,sizeof(colorspinspin)*loc_vol);
}}
//do the opposite
THREADABLE_FUNCTION_2ARG(revert_prop_from_new_contraction, colorspinspin*,prop, colorspinspin*,aux)
{
  GET_THREAD_ID();
  
  int spat_color_vol=3*loc_vol/loc_size[0];
  
  //fill aux
  for(int t=0;t<loc_size[0];t++)
    for(int idirac_ri=0;idirac_ri<32;idirac_ri++)
      NISSA_PARALLEL_LOOP(ispat_color,0,spat_color_vol)
        ((double*)aux)[idirac_ri+32*(ispat_color+spat_color_vol*t)]=
        ((double*)prop)[ispat_color+spat_color_vol*(idirac_ri+32*t)];
  THREAD_BARRIER();
  
  //copy aux to prop
  parallel_memcpy(prop,aux,sizeof(colorspinspin)*loc_vol);
}}

//compute 2 pts
void calculate_all_2pts(int ism_lev_so,int ism_lev_si)
{
  //take nsme
  master_printf("\nSink Smearing level: %d\n",ism_lev_si);
  int nsme=gaussian_niter[ism_lev_si];
  if(ism_lev_si>0) nsme-=gaussian_niter[ism_lev_si-1];

  //smear and put in the special layout
  smear_time-=take_time();
  for(int r=0;r<2;r++)
    for(int imu=0;imu<4;imu++)
      for(int inu=0;inu<4;inu++)
	{
	  gaussian_smearing(S0[r][imu][inu],S0[r][imu][inu],sme_conf,gaussian_kappa,nsme);
	  prepare_prop_for_new_contraction(S0[r][imu][inu],temp_transp);
	}
  master_printf("\n");
  
  //take intermediate timing
  double temp_time=take_time();
  smear_time+=temp_time;
  corr_time-=temp_time;
  
  //compute correlations
  for(int imu1_so=0;imu1_so<4;imu1_so++)
    for(int imu2_so=0;imu2_so<4;imu2_so++)
      for(int imu1_si=0;imu1_si<4;imu1_si++)
	for(int imu2_si=0;imu2_si<4;imu2_si++)
	  if(two_pts_comp[imu1_so][imu2_so][imu1_si][imu2_si].ncorr)
	    {
	      //master_printf("%d%d%d%d\n",imu1_so,imu2_so,imu1_si,imu2_si);
	      //if(rank==0)
	      //{
	      //for(int icorr=0;icorr<two_pts_comp[imu1_so][imu2_so][imu1_si][imu2_si].ncorr;icorr++)
	      //printf("%d %s\n",icorr,two_pts_comp[imu1_so][imu2_so][imu1_si][imu2_si].corr_name[icorr].c_str());
	      //two_pts_comp[imu1_so][imu2_so][imu1_si][imu2_si].print();
	      //}
	      double *out=loc_corr+glb_size[0]*((ism_lev_so*nsm_lev+ism_lev_si)*
			   ncorr_tot+two_pts_comp_off[imu1_so][imu2_so][imu1_si][imu2_si]);
	      int slice_vol=3*loc_vol/loc_size[0];
	      
	      //contract
	      for(int r=0;r<2;r++)
		{
		  double *S1=(double*)(S0[r][imu1_so][imu1_si]);
		  double *S2=(double*)(S0[r][imu2_so][imu2_si]);
		  two_pts_comp[imu1_so][imu2_so][imu1_si][imu2_si].
		    summ_the_loc_forw_back_contractions(out,S1,S2,slice_vol,twall);
		  ncorr_computed_tot+=two_pts_comp[imu1_so][imu2_so][imu1_si][imu2_si].ncorr;
		}
	    }
  
  //come back to normal layout
  for(int r=0;r<2;r++)
    for(int imu=0;imu<4;imu++)
      for(int inu=0;inu<4;inu++)
	revert_prop_from_new_contraction(S0[r][imu][inu],temp_transp);
  
  corr_time+=take_time();
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

  master_printf("Remaining time: %lg sec\n",left_time);
  master_printf("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
  if(enough_time) master_printf("Continuing with next conf!\n");
  else master_printf("Not enough time, exiting!\n");
  
  return enough_time;
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //initialize the program
  if(narg<2) crash("Use: %s input_file",arg[0]);
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
	  twall=(int)rnd_get_unif(&glb_rnd_gen,0,glb_size[0]-1);
	  master_printf("Source %d on t=%d\n",isource,twall);
	  generate_spindiluted_source(ori_source,rnd_type_map[4],twall);
	  
	  //loop on smearing of the source
	  for(int ism_lev_so=0;ism_lev_so<nsm_lev;ism_lev_so++)
	    {
	      //compute S0 propagator smearing the config the appropriate number of time the source
	      calculate_S0(ism_lev_so);
	      
	      //loop on the smearing of the sink
	      for(int ism_lev_si=0;ism_lev_si<nsm_lev;ism_lev_si++) calculate_all_2pts(ism_lev_so,ism_lev_si);
	    }
	}
      
      //open output file
      char path[1024];
      sprintf(path,"%s/2pts_corr",outfolder);
      FILE *fout=open_text_file_for_output(path);
      
      //reduce and write
      MPI_Reduce(loc_corr,glb_corr,ndoubles,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      double_vector_prod_double(glb_corr,glb_corr,0.5/nsources,ndoubles);

      if(rank==0)
	{
	  int bin_write=0;
	  if(bin_write)
	    {
	      int nw=fwrite(glb_corr,sizeof(double),ndoubles,fout);
	      if(nw!=ndoubles) crash("wanted to write %d, obtained %d",ndoubles,nw);
	    }
	  else
	    {
	      for(int ism_lev_so=0;ism_lev_so<nsm_lev;ism_lev_so++)
		for(int imu1_so=0;imu1_so<4;imu1_so++)
		  for(int imu2_so=0;imu2_so<4;imu2_so++)
		    for(int ism_lev_si=0;ism_lev_si<nsm_lev;ism_lev_si++)
		      for(int imu1_si=0;imu1_si<4;imu1_si++)
			for(int imu2_si=0;imu2_si<4;imu2_si++)
			  if(two_pts_comp[imu1_so][imu2_so][imu1_si][imu2_si].ncorr)
			    {
			      double *out=glb_corr+glb_size[0]*
				((ism_lev_so*nsm_lev+ism_lev_si)*ncorr_tot+
				 two_pts_comp_off[imu1_so][imu2_so][imu1_si][imu2_si]);
			      
			      fprintf(fout," # sm_so %d mu1_so %d mu2_so %d "  ,ism_lev_so,imu1_so,imu2_so);
			      fprintf(fout,   "sm_si %d mu1_si %d mu2_si %d \n",ism_lev_si,imu1_si,imu2_si);
			      fprintf(fout,"\n");
			      
			      for(int icorr=0;icorr<two_pts_comp[imu1_so][imu2_so][imu1_si][imu2_si].ncorr;icorr++)
				{
				  fprintf(fout," # %s \n",two_pts_comp[imu1_so][imu2_so][imu1_si][imu2_si]
					  .corr_name[icorr].c_str());
				  for(int t=0;t<glb_size[0];t++) fprintf(fout,"%+016.16lg\n",out[icorr*glb_size[0]+t]);
				  fprintf(fout,"\n");
				}
			    }
	    }
	  fclose(fout);
	}
      
      nanalyzed_conf++;
      enough_time=check_remaining_time();
    }
  
  if(iconf==ngauge_conf) master_printf("Finished all the conf!\n");
  
  tot_prog_time+=take_time();
  close_semileptonic();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
