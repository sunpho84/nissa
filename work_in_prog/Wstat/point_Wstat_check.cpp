#include <math.h>

#include "nissa.h"

//gauge info
char conf_path[1024],outfolder[1024];
quad_su3 *conf,*hyp_conf;
double kappa;
double put_theta[4],old_theta[4]={0,0,0,0};

//list of masses and theta
int nmass;
double *mass;

//source data
coords source_coord;
spincolor *source;
su3spinspin *original_source;

int nprop;
su3spinspin **S[2];
int ncgm_solution;
spincolor **cgm_solution,*temp_vec[2];

//cgm inverter parameters
double *stopping_residues;
int niter_max=100000;

//two points contractions
int ncontr_2pts;
complex *contr_2pts;
int *op1_2pts,*op2_2pts;


void make_point_test(quad_su3 *conf)
{
  su3spinspin *prop=nissa_malloc("prop",loc_vol+bord_vol,su3spinspin);
  spincolor *in=nissa_malloc("in",loc_vol+bord_vol,spincolor);
  spincolor *out=nissa_malloc("out",loc_vol+bord_vol,spincolor);
  
  ////////////////////////////////////////////////
  
  //choose a site
  int ivol=glblx_of_coord_list(1,1,1,1);
  
  //compute the prop
  compute_Wstat_prop_point(prop,conf,0,source_coord);
  
  //take a dirac-color index
  int id=0,ic=0;
  get_spincolor_from_su3spinspin(in,prop,id,ic);
  master_printf("Prop at site %d:\n",ivol);
  spincolor_print(in[ivol]);
  
  //verify that gamma4*prop=prop
  safe_dirac_prod_spincolor(out,base_gamma[4],in);  
  master_printf("Gamma4*Prop at site %d:\n",ivol);
  spincolor_print(out[ivol]);
  
  //verify that W*in=delta
  apply_Wstat(out,conf,in,0,0);
  nissa_loc_vol_loop(x) if(glb_coord_of_loclx[x][0]==0) out[x][id][ic][0]-=1;
  master_printf("W*Prop at site %d:\n",ivol);
  spincolor_print(out[ivol]);
  
  //verify the norm of the out
  nissa_loc_vol_loop(x)
    {
      double n=0;
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    n+=sqr(out[x][id][ic][ri]);
      master_printf("%d %d %lg\n",glb_coord_of_loclx[x][0],x,n);
      spincolor_print(out[x]);
    }
  master_printf("Error: %lg\n",glb_reduce_double(double_vector_loc_scalar_prod((double*)out,(double*)out,24*loc_vol)));
  
  ////////////////////////////////////////////////
  
  nissa_free(out);
  nissa_free(in);
  nissa_free(prop);
}

void init(char *input_path)
{
  open_input(input_path);

  // 1) Read information about the gauge conf
  
  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Kappa
  read_str_double("Kappa",&kappa);
  
  // 2) Read list of masses and of thetas

  read_list_of_double_pairs("MassResidues",&nmass,&mass,&stopping_residues);
  mass=(double*)realloc((void*)mass,(nmass+1)*sizeof(double));
  stopping_residues=(double*)realloc((void*)stopping_residues,(nmass+1)*sizeof(double));

  //Add the static point                                                                                                                                                                    
  mass[nmass]=1000000;
  stopping_residues[nmass]=1.e-300;

  // 3) contraction list for two points
  
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
  
  master_printf("\n");
  
  ////////////////////////////////////// end of input reading/////////////////////////////////
  
  //allocate gauge conf, hypped conf and all the needed spincolor and propagators
  conf=nissa_malloc("or_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  hyp_conf=nissa_malloc("hyp_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //Allocate all the S su3spinspin vectors
  nprop=nmass+1;
  S[0]=nissa_malloc("S[0]",nprop,su3spinspin*);
  S[1]=nissa_malloc("S[1]",nprop,su3spinspin*);
  for(int iprop=0;iprop<nprop;iprop++)
    {
      S[0][iprop]=nissa_malloc("S[0]",loc_vol,su3spinspin);
      S[1][iprop]=nissa_malloc("S[1]",loc_vol,su3spinspin);
    }
  
  //Allocate nmass spincolors, for the cgm solutions
  ncgm_solution=nmass;
  cgm_solution=nissa_malloc("cgm_solution",ncgm_solution,spincolor*);
  for(int imass=0;imass<ncgm_solution;imass++) cgm_solution[imass]=nissa_malloc("cgm_solution",loc_vol+bord_vol,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol,spincolor);
  
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  original_source=nissa_malloc("original_source",loc_vol,su3spinspin);
}

//generate the source
void generate_source()
{generate_delta_source(original_source,source_coord);}

//calculate the standard propagators
void calculate_S()
{
  //loop over the source dirac index
  for(int ic=0;ic<3;ic++)
    for(int id=0;id<4;id++)
      { 
	//put the g5
	nissa_loc_vol_loop(ivol)
	  {
	    get_spincolor_from_su3spinspin(source[ivol],original_source[ivol],id,ic);
	    safe_dirac_prod_spincolor(source[ivol],&(base_gamma[5]),source[ivol]);
	  }
	set_borders_invalid(source);
	
	inv_tmQ2_cgm(cgm_solution,conf,kappa,mass,nmass,niter_max,stopping_residues,source);
	master_printf("Finished the inversion of S, color index %d dirac index %d\n",id,ic);
	
	//reconstruct the doublet
	for(int imass=0;imass<nmass;imass++)
	  {
	    reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass[imass],cgm_solution[imass]);
	    master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	    for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
	      put_spincolor_into_su3spinspin(S[r][imass],temp_vec[r],id,ic);
	  }
      }

  //rotate to physical basis
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int imass=0;imass<nmass;imass++) //put the (1+ig5)/sqrt(2) factor
      rotate_vol_su3spinspin_to_physical_basis(S[r][imass],!r,!r);
  
  master_printf("Propagators rotated\n");

  compute_Wstat_prop_point(S[0][nmass],conf,0,source_coord);
  vector_copy(S[1][nmass],S[0][nmass]);
  
  complex cc={0,0};
  for(int ic=0;ic<3;ic++) complex_summassign(cc,S[0][0][0][ic][ic][0][0]);
  printf("%16.16lg %16.16lg\n",cc[0],cc[1]);
  master_printf("\n");
}

//read the conf and setup it
void setup_conf()
{
  //Gauge path
  read_str(conf_path,1024);
      
  //Source coord
  read_int(&(source_coord[0]));
  read_int(&(source_coord[1]));
  read_int(&(source_coord[2]));
  read_int(&(source_coord[3]));
      
  //Out folder
  read_str(outfolder,1024);
  create_dir(outfolder);
  
  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_ildg_gauge_conf(conf,conf_path);
  //landau_gauge_fix(conf,conf,1.e-28);
  
  //put the anti-periodic condition on the temporal border
  old_theta[0]=old_theta[1]=old_theta[2]=old_theta[3]=0;
  put_theta[0]=1;put_theta[1]=put_theta[2]=put_theta[3]=0;
  //adapt_theta(conf,old_theta,put_theta,1,1);
  
  //prepare the hypped version
  double hyp_alpha0=1,hyp_alpha1=1,hyp_alpha2=0.5;
  hyp_smear_conf_dir(hyp_conf,conf,hyp_alpha0,hyp_alpha1,hyp_alpha2,0);
  
  //compute plaquette
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  master_printf("hypped plaq: %.18g\n",global_plaquette_lx_conf(hyp_conf));
}

//Finalization
void close()
{
  nissa_free(conf);nissa_free(hyp_conf);
  for(int iprop=0;iprop<nprop;iprop++){nissa_free(S[0][iprop]);nissa_free(S[1][iprop]);}
  nissa_free(S[0]);nissa_free(S[1]);
  nissa_free(temp_vec[0]);nissa_free(temp_vec[1]);
  nissa_free(contr_2pts);
  nissa_free(op1_2pts);nissa_free(op2_2pts);
  for(int imass=0;imass<ncgm_solution;imass++) nissa_free(cgm_solution[imass]);
  nissa_free(cgm_solution);
  nissa_free(source);nissa_free(original_source);
  close_nissa();
}

//12 index prop                                                                                                                                                                                 
void site_trace_g_ccss_dag_g_ccss_mod(complex c,dirac_matr *g1,su3spinspin s1,dirac_matr *g2,su3spinspin s2)
{
  c[0]=c[1]=0;
  //Color loop                                                                                                                                                                                  
  for(int ic1=0;ic1<3;ic1++)
    for(int ic2=0;ic2<3;ic2++)
      {
        spinspin t1,t2;

        unsafe_dirac_prod_spinspin_dag(t1,g1,s1[ic2][ic1]);
        unsafe_dirac_prod_spinspin_dag(t2,g2,s2[ic2][ic1]);

	summ_the_trace_prod_spinspins(c,t1,t2);
      }
}

//Trace the product of gamma1 * su3spinspin1^dag * gamma2 * su3spinspin2,                                                                                                                       
void trace_g_ccss_dag_g_ccss_mod(complex *glb_c,dirac_matr *g1,su3spinspin *s1,dirac_matr *g2,su3spinspin *s2,const int ncontr)
{
  //Allocate a contiguous memory area where to store local node results                                                                                                                         
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      verbosity_lv3_master_printf("Contraction %d/%d\n",icontr+1,ncontr);

      //Local loop                                                                                                                                                                              
      nissa_loc_vol_loop(ivol)
      {
	int glb_t=glb_coord_of_loclx[ivol][0];

	complex ctemp;
	site_trace_g_ccss_dag_g_ccss_mod(ctemp,&(g1[icontr]),s1[ivol],&(g2[icontr]),s2[ivol]);
	complex_summassign(loc_c[icontr*glb_size[0]+glb_t],ctemp);
      }
    }

  //Final reduction                                                                                                                                                                             
  verbosity_lv3_master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  verbosity_lv3_master_printf("Reduction done\n");

  nissa_free(loc_c);
}


//This function takes care to make the revert on the FIRST spinor, putting the needed gamma5                                                                                                    
void meson_two_points_Wilson_prop_mod(complex *corr,int *list_op1,su3spinspin *s1,int *list_op2,su3spinspin *s2,int ncontr)
{
  //Temporary vectors for the internal gamma                                                                                                                                                    
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Put the two gamma5 needed for the revert of the first spinor                                                                                                                            
      dirac_prod(&(t1[icontr]), &(base_gamma[list_op1[icontr]]),&(base_gamma[5]));
      dirac_prod(&(t2[icontr]), &(base_gamma[5]),&(base_gamma[list_op2[icontr]]));
    }

  //Call the routine which perform the contraction                                                                                                                                              
  trace_g_ccss_dag_g_ccss_mod(corr,t1,s1,t2,s2,ncontr);
}


void calculate_all_2pts()
{
  char path[1024];
  sprintf(path,"%s/2pts",outfolder);
  FILE *fout=open_text_file_for_output(path);
  
  for(int im2=nmass;im2<nmass+1;im2++)
    for(int r2=0;r2<1;r2++)
      for(int im1=0;im1<nmass;im1++)
	for(int r1=0;r1<1;r1++)
	  {
	    //header
	    master_fprintf(fout," # m1=%lg res1=%lg r1=%d, m2=%lg res2=%lg r2=%d\n",
			   mass[im1],stopping_residues[im1],r1,
			   mass[im2],stopping_residues[im2],r2);
	    
	    //compute contractions
	    meson_two_points_Wilson_prop(contr_2pts,op1_2pts,S[r1][im1],op2_2pts,S[r2][im2],ncontr_2pts);
	    
	    //write 
	    print_contractions_to_file(fout,ncontr_2pts,op1_2pts,op2_2pts,contr_2pts,source_coord[0],"",1.0);
	  }
  
  if(rank==0) fclose(fout);
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();
  
  //initialize the program
  if(narg<2) crash("Use: %s input_file",arg[0]);
  init(arg[1]);
  
  //load the conf and generate the source
  setup_conf();
  generate_source();
  
  //compute S propagators
  calculate_S();
  
  calculate_all_2pts();
  
  close();

  return 0;
}
