#include "nissa.hpp"

/*
            	 
           .....           
         ..     ..         
        .         .        
   op  X           .  S
        .         .        
         ..     ..         
           .....           
                      
*/

int ncontr,nch_contr;
complex *contr,*ch_contr;
int *op,*ch_op;

int wall_time;
quad_su3 *conf,*sme_conf;
as2t_su3 *Pmunu;

colorspinspin *source,*ch_prop,***S,*edm_prop[3];
spincolor *inv_source,*temp_vec[2],**QQ;

int starting_source,ending_source,noise_type;
int compute_edm;
int source_pos[4];

int nmass;
double kappa,*mass;

double *stopping_residues;
int niter_max=100000000;

int ngauge_conf,nanalyzed_conf=0;
char conf_path[1024],outfolder[1024];

double tot_time=0,inv_time=0,contr_time=0;
int ncontr_tot,ninv_tot;

//smearing parameters
double jacobi_kappa,ape_alpha;
int ape_niter;
int jacobi_niter;

const double rad2=1.414213562373095048801688724209;

//contract a source with a propagator putting the passed list of operators
void contract_with_source(complex *corr,colorspinspin *prop,int *list_op,colorspinspin *source,int ncontr)
{
  //Temporary vector for the internal matrices
  dirac_matr t1[ncontr];
  dirac_matr t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Init the second 
      t1[icontr]=base_gamma[list_op[icontr]];
      t2[icontr]=base_gamma[0];
    }

  //Call for the routine which does the real contraction
  trace_g_sdag_g_s(corr,t2,source,t1,prop,ncontr);
}

//This perform the contraction using the so-callled (by Silvano) Chris Micheal trick
//otherwise called "twisted trick"
void contract_with_chris_michael(complex *corr,int *list_op,colorspinspin *s1,colorspinspin *s2,int ncontr,double mass)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {      
      //copy the dirac matrices
      t1[icontr]=base_gamma[list_op[icontr]];
      t2[icontr]=base_gamma[0];
    }
  
  //Call for the routine which does the real contraction
  trace_g_sdag_g_s(corr,t2,s2,t1,s1,ncontr);
  
  //multiply by 2*mass
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int t=0;t<glb_size[0];t++)
      for(int ri=0;ri<2;ri++)
	(corr+icontr*glb_size[0])[t][ri]*=2*mass;  
}
	  
//print all the passed contractions to the file
void print_bubbles_to_file(FILE *fout,int ncontr,int *op,complex *contr,const char *tag)
{
  master_fprintf(fout,"\n");
  
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      master_fprintf(fout," # %s%s\n\n",tag,gtag[op[icontr]]);
      for(int t=0;t<glb_size[0];t++)
	master_fprintf(fout,"%+016.16g\t%+016.16g\n",(contr+icontr*glb_size[0])[t][0],(contr+icontr*glb_size[0])[t][1]);
      master_fprintf(fout,"\n");
    }
  master_fprintf(fout,"\n");
}

void initialize_bubbles(char *input_path)
{
  //take initial time
  tot_time=-take_time();
  
  open_input(input_path);
  
  //read lattice size
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  //init the grid 
  init_grid(T,L);
  
  //read the walltime
  read_str_int("WallTime",&wall_time);
  //read kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read information about the source
  
  //read the seed
  int seed;
  read_str_int("Seed",&seed);
  //read the number of the starting sources
  read_str_int("StartingSource",&starting_source);
  //read the number of the ending sources
  read_str_int("EndingSource",&ending_source);
  //read the noise type
  read_str_int("NoiseType",&noise_type);

  // 3) Smearing parameters
  
  read_str_double("ApeAlpha",&ape_alpha);
  read_str_int("ApeNiter",&ape_niter);
  read_str_double("JacobiKappa",&jacobi_kappa);
  read_str_int("JacobiNiter",&jacobi_niter);

  // 4) Masses
  read_list_of_double_pairs("MassResidues",&nmass,&mass,&stopping_residues);

  // 5) Contractions
  
  //read the list of correlations and the list of operators
  read_list_of_ints("NContr",&ncontr,&op);
  read_list_of_ints("NChromoContr",&nch_contr,&ch_op);
  //read whether we want to compute EDM
  read_str_int("ComputeEDM",&compute_edm);
  
  read_str_int("NGaugeConf",&ngauge_conf);
  
  ////////////////////////////// allocation & initialization //////////////////////////
  
  //start loc rnd gen
  start_loc_rnd_gen(seed);
  
  //allocate contr
  contr=nissa_malloc("contr",ncontr*glb_size[0],complex);
  if(nch_contr!=0) ch_contr=nissa_malloc("ch_contr",nch_contr*glb_size[0],complex);
  
  //allocate gauge conf and Pmunu
  conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  sme_conf=nissa_malloc("sme_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  if(compute_edm) Pmunu=nissa_malloc("clover",loc_vol,as2t_su3);
  
  //allocate the source and the chromo-prop
  source=nissa_malloc("source",loc_vol+bord_vol,colorspinspin);
  inv_source=nissa_malloc("inv_source",loc_vol+bord_vol,spincolor);
  if(nch_contr!=0) ch_prop=nissa_malloc("chromo_prop",loc_vol,colorspinspin);
  if(compute_edm) for(int idir=0;idir<3;idir++) edm_prop[idir]=nissa_malloc("edm_prop",loc_vol,colorspinspin);
  
  //temporary vectors
  for(int r=0;r<2;r++) temp_vec[r]=nissa_malloc("temp_vec",loc_vol,spincolor);
  
  //allocate spinors for the cgm
  S=nissa_malloc("S",nmass,colorspinspin**);
  QQ=nissa_malloc("QQ",nmass,spincolor*);
  for(int imass=0;imass<nmass;imass++)
    {
      S[imass]=nissa_malloc("S[imass]",2,colorspinspin*);
      for(int r=0;r<2;r++) S[imass][r]=nissa_malloc("S",loc_vol,colorspinspin);
      QQ[imass]=nissa_malloc("QQ[i]",loc_vol+bord_vol,spincolor);
    }
}

//wrapper
void generate_volume_source(int isource)
{
  enum rnd_type type[5]={RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_Z2,RND_Z2,RND_Z4};
  generate_spindiluted_source(source,type[noise_type],-1);
  
  for(int idso=0;idso<4;idso++)
    {
      NISSA_LOC_VOL_LOOP(ivol)
	get_spincolor_from_colorspinspin(inv_source[ivol],source[ivol],idso);
      set_borders_invalid(inv_source);
      
      jacobi_smearing(inv_source,inv_source,sme_conf,jacobi_kappa,jacobi_niter);
      NISSA_LOC_VOL_LOOP(ivol)
	put_spincolor_into_colorspinspin(source[ivol],inv_source[ivol],idso);
    }
}

//compute propagators
void calculate_S()
{
  inv_time-=take_time();
  
  //Invert
  for(int idso=0;idso<4;idso++)
    {
      //prepare the source
      NISSA_LOC_VOL_LOOP(ivol)
	{
	  get_spincolor_from_colorspinspin(inv_source[ivol],source[ivol],idso);
	  
	  //put the g5
	  for(int idsi=2;idsi<4;idsi++)
	    for(int icol=0;icol<3;icol++)
	      for(int ri=0;ri<2;ri++)
		inv_source[ivol][idsi][icol][ri]=-inv_source[ivol][idsi][icol][ri];
	}
      set_borders_invalid(inv_source);

      inv_tmQ2_cgm(QQ,conf,kappa,mass,nmass,niter_max,stopping_residues,inv_source);
      ninv_tot++;
      //put the solution inside the S vector
      for(int imass=0;imass<nmass;imass++) 
	{
	  reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass[imass],QQ[imass]);
	  for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
	    NISSA_LOC_VOL_LOOP(ivol)
	      put_spincolor_into_colorspinspin(S[imass][r][ivol],temp_vec[r][ivol],idso);
	}
    }
  
  //rotate to physical basis
  for(int imass=0;imass<nmass;imass++) //put the (1+ig5)/sqrt(2) factor
    for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
      rotate_vol_colorspinspin_to_physical_basis(S[imass][r],!r,!r);
  
  inv_time+=take_time();
}

//Apply the dipole operator on a su3spinspin
void apply_dipole_operator(colorspinspin *S_out,colorspinspin *S_in,int dir)
{
  NISSA_LOC_VOL_LOOP(ivol)
    {
      int coor=(glb_coord_of_loclx[ivol][dir]-source_pos[dir]+glb_size[dir])%glb_size[dir];
      
      if(coor> glb_size[dir]/2) coor-=glb_size[dir]; //minor distance
      if(coor==glb_size[dir]/2) coor=0; //take away the border (Simpson)
      
      for(int icsi=0;icsi<3;icsi++)
	for(int idso=0;idso<4;idso++)
	  for(int idsi=0;idsi<4;idsi++)
	    for(int ri=0;ri<2;ri++)
	      S_out[ivol][icsi][idsi][idso][ri]=S_in[ivol][icsi][idsi][idso][ri]*coor;
    }
}

//contract
void calculate_all_contractions(int isource)
{
  contr_time-=take_time();
  
  char outp[1024];
  sprintf(outp,"%s/bubble_%04d",outfolder,isource);
  FILE *fout=open_text_file_for_output(outp);
  
  for(int imass=0;imass<nmass;imass++)
    for(int r=0;r<2;r++)
      {
	//header
	master_fprintf(fout," # mass=%g r=%d\n",mass[imass],r);
	
	//simple contractions
	contract_with_source(contr,S[imass][r],op,source,ncontr);
	print_bubbles_to_file(fout,ncontr,op,contr,"");
	
	ncontr_tot+=ncontr;
	
	//apply the chromo magnetic operator to the second spinor
	if(nch_contr!=0)
	  {
	    unsafe_apply_chromo_operator_to_colorspinspin(ch_prop,Pmunu,S[imass][r]);
	    contract_with_source(ch_contr,ch_prop,ch_op,source,nch_contr);
	    print_bubbles_to_file(fout,nch_contr,ch_op,ch_contr,"CHROMO-");

	    ncontr_tot+=nch_contr;
	  }
	
	//edm contractions
	if(compute_edm)
	  {
	    for(int idir=0;idir<3;idir++) apply_dipole_operator(edm_prop[idir],S[imass][r],idir+1);
	    
	    for(int idir=0;idir<3;idir++)
	      {
		int ninsertions=1;
		int dipole_insertion[1]={4};
		contract_with_source(contr,edm_prop[idir],dipole_insertion,source,ninsertions);
		char tag[1024];
		sprintf(tag,"EDM_%d",idir+1);
		print_bubbles_to_file(fout,ninsertions,dipole_insertion,contr,tag);
	      }
	  }
	
	//do everything with the chris-michael trick
	if(r==0)
	  {
	    contract_with_chris_michael(contr,op,S[imass][r],S[imass][r],ncontr,mass[imass]);
	    print_bubbles_to_file(fout,ncontr,op,contr,"CHRIS-");
	    
	    ncontr_tot+=ncontr;
	    
	    if(nch_contr!=0)
	      {
		contract_with_chris_michael(ch_contr,ch_op,ch_prop,S[imass][r],nch_contr,mass[imass]);
		print_bubbles_to_file(fout,nch_contr,ch_op,ch_contr,"CHROMO-CHRIS-");
		
		ncontr_tot+=nch_contr;
	      }
	  }
      }

  contr_time+=take_time();

  if(rank==0) fclose(fout);
}

void close_bubbles()
{
  close_input();

  //take final time
  tot_time+=take_time();
  
  master_printf("\n");
  master_printf("Total time: %g, of which:\n",tot_time);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d contr. (%2.2gs avg)\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
  
  ///////////////////////////////////////////
  
  master_printf("\nEverything ok, exiting!\n");
  
  nissa_free(contr);
  nissa_free(conf);
  nissa_free(sme_conf);
  if(compute_edm) nissa_free(Pmunu);
  if(nch_contr!=0) nissa_free(ch_contr);
  
  nissa_free(source);
  for(int imass=0;imass<nmass;imass++)
    {
      nissa_free(QQ[imass]);
      for(int r=0;r<2;r++) nissa_free(S[imass][r]);
      nissa_free(S[imass]);
    }
  nissa_free(S);
  nissa_free(QQ);
  
  nissa_free(inv_source);
  if(nch_contr!=0) nissa_free(ch_prop);
  if(compute_edm) for(int idir=0;idir<3;idir++) nissa_free(edm_prop[idir]);
  for(int r=0;r<2;r++) nissa_free(temp_vec[r]);
  
  close_nissa();
}

int read_conf_parameters(int *iconf)
{
  int ok_conf;

  do
    {
      //read path of conf and output
      read_str(conf_path,1024);
      read_str(outfolder,1024);
      if(compute_edm)
	{
	  //read the position of the source (for the edm)
	  expect_str("SourcePosition");
	  for(int idir=0;idir<4;idir++)
	    {
	      read_int(&(source_pos[idir]));
	      master_printf("%d ",source_pos[idir]);
	    }
	}

      //check if the conf exist
      master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
      ok_conf=!(dir_exists(outfolder));
      if(ok_conf)
        {
          int ris=create_dir(outfolder);
          if(ris==0) master_printf(" Output path \"%s\" not present: configuration \"%s\" not yet analyzed, starting.\n",outfolder,conf_path);
          else
            {
              ok_conf=0;
              master_printf(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
            }
        }
      else
        master_printf(" Output path \"%s\" already present: configuration \"%s\" already analyzed, skipping.\n",outfolder,conf_path);
      (*iconf)++;
    }
  while(!ok_conf && (*iconf)<ngauge_conf);
  
  master_printf("\n");
  
  return ok_conf;
}

//read the conf and setup it
void setup_conf()
{
  //load conf, calculate plaquette and Pmunu
  read_ildg_gauge_conf(conf,conf_path);
  ape_spatial_smear_conf(sme_conf,conf,ape_alpha,ape_niter);
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  master_printf("smerded plaq: %.18g\n",global_plaquette_lx_conf(sme_conf));

  if(nch_contr!=0) Pmunu_term(Pmunu,conf);
      
  //put the anti-periodic condition on the temporal border
  double put_theta[4],old_theta[4]={0,0,0,0};
  put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
  adapt_theta(conf,old_theta,put_theta,1,1);
}

//check if the time is enough
int check_remaining_time()
{
  int enough_time;

  //check remaining time                                                                                                                                                                        
  double temp_time=take_time()+tot_time;
  double ave_time=temp_time/nanalyzed_conf;
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
  
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  initialize_bubbles(arg[1]);
  
  //loop over the configs
  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf && enough_time && read_conf_parameters(&iconf))
    {
      setup_conf();
      
      //loop over the sources
      for(int isource=starting_source;isource<ending_source;isource++)
	{
	  generate_volume_source(isource);
	  calculate_S();
	  calculate_all_contractions(isource);
	}
      
      //pass to the next conf if there is enough time
      nanalyzed_conf++;
      enough_time=check_remaining_time();
    }
  if(iconf==ngauge_conf) master_printf("Finished all the conf!\n");
  
  close_bubbles();
  
  return 0;
}
