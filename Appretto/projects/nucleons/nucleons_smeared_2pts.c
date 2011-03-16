#include <mpi.h>
#include <math.h>
#include <lemon.h>

#include "appretto.h"

//configuration
int nconf;
double put_theta[4]={0,0,0,0};
char **conf_path,**out_path;
quad_su3 *conf,*ori_conf;
double mass;
double kappa;

//source
int source_pos[4];
spincolor *source;
su3spinspin *original_source;

//smearing parameters
double jacobi_kappa;
int jacobi_niter;

//the propagators
su3spinspin *S0[2];

//inverter
int nitermax;
double residue;
spincolor *solDD,*sol_reco[2];

//insertion
int tseparation;
int tsink;

spinspin Proj[3]; //projectors over N and N*, and 00 compont of N (in the spinorial representation)
spinspin C5; //C*gamma5

//two points contractions
int nproton_2pt_contr=32;      //VA     AV       VV       AA       TT        BB        TB        BT
int list_2pt_op1[32]={5,0,5,0, 1,2,3,4, 6,7,8,9, 1,2,3,4, 6,7,8,9, 10,11,12, 13,14,15, 10,11,12, 13,14,15};
int list_2pt_op2[32]={0,5,5,0, 6,7,8,9, 1,2,3,4, 1,2,3,4, 6,7,8,9, 10,11,12, 13,14,15, 13,14,15, 10,11,12};

//                      e_00x   e_01x    e_02x     e_10x    e_11x   e_12x     e_20x   e_21x    e_22x
int epsilon[3][3][3]={{{0,0,0},{0,0,1},{0,-1,0}},{{0,0,-1},{0,0,0},{1,0,0}},{{0,1,0},{-1,0,0},{0,0,0}}};

//timings
double tinv=0,tcontr=0,tot_time=0;

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
  ori_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"ori_conf");

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
  //Put border condition and communicate

  //Kappa
  read_str_double("Kappa",&(kappa));
  
  // 2) Source position and smearing parameters
  expect_str("SourcePosition");
  if(rank==0) printf("Source position: ");
  for(int idir=0;idir<4;idir++) 
    {
      read_int(&(source_pos[idir]));
      if(rank==0) printf("%d ",source_pos[idir]);
    }
  if(rank==0) printf("\n");
  read_str_double("JacobiKappa",&jacobi_kappa);
  read_str_int("JacobiNiter",&jacobi_niter);

  // 3) mass list (to be added)
  read_str_double("Mass",&mass);

  // 4) inverter

  //Residue
  read_str_double("Residue",&residue);
  //Number of iterations
  read_str_int("NiterMax",&nitermax);
  
  close_input();

  ///////////////////// Allocate the various spinors ///////////////////////
  
  original_source=allocate_su3spinspin(loc_vol,"original_source");
  
  source=allocate_spincolor(loc_vol+loc_bord,"source");
  solDD=allocate_spincolor(loc_vol+loc_bord,"solDD");

  sol_reco[0]=allocate_spincolor(loc_vol,"solution_reco[0]");
  sol_reco[1]=allocate_spincolor(loc_vol,"solution_reco[1]");

  S0[0]=allocate_su3spinspin(loc_vol,"S0[0]");
  S0[1]=allocate_su3spinspin(loc_vol,"S0[1]");
}

//read a configuration and put anti-periodic condition at the slice tsource-1
void read_conf_and_put_antiperiodic(quad_su3 *conf,char *conf_path,int tsource)
{
  read_local_gauge_conf(ori_conf,conf_path);
  memcpy(conf,ori_conf,sizeof(quad_su3)*loc_vol);

  //commmunicate borders
  communicate_gauge_borders(conf);  
  communicate_gauge_edges(conf);
  
  //calculate plaquette
  double gplaq=global_plaquette(conf);
  if(rank==0) printf("plaq: %.18g\n",gplaq);

  //Put the anti-periodic condition on the temporal border
  put_theta[0]=1;
  put_boundaries_conditions(conf,put_theta,1,1);

  //re-communicate borders
  communicate_gauge_borders(conf);  
  communicate_gauge_edges(conf);
}

//create a 12 index point source
void prepare_source()
{
  int isloc=1;

  int lx[4];

  memset(original_source,0,sizeof(spincolor)*loc_vol);

  for(int idir=0;idir<4;idir++)
    {
      lx[idir]=source_pos[idir]-proc_coord[idir]*loc_size[idir];
      isloc=isloc && (lx[idir]>=0 && lx[idir]<loc_size[idir]);
    }

  int ivol=loclx_of_coord(lx);

  if(isloc)
    for(int ic=0;ic<3;ic++)
      for(int id=0;id<4;id++)
	original_source[ivol][ic][ic][id][id][0]=1;
}      

//perform the first inversion to produce the S0 for u and d
void calculate_S0()
{
  spincolor *temp_source=allocate_spincolor(loc_vol+loc_bord,"temporary_source");
  
  for(int ic_sour=0;ic_sour<3;ic_sour++)
    for(int id_sour=0;id_sour<4;id_sour++)
      { //take the source and put g5
	for(int ivol=0;ivol<loc_vol;ivol++)
	  {
	    get_spincolor_from_su3spinspin(temp_source[ivol],original_source[ivol],id_sour,ic_sour);
	    for(int id_sink=2;id_sink<4;id_sink++)
	      for(int ic_sink=0;ic_sink<3;ic_sink++)
		for(int ri=0;ri<2;ri++)
		  temp_source[ivol][id_sink][ic_sink][ri]=-temp_source[ivol][id_sink][ic_sink][ri];
	  }
		
	if(rank==0) printf("\n(S0) source index: id=%d, ic=%d\n",id_sour,ic_sour);
	
	//smerd the source
	dina_smearing(source,temp_source,conf,jacobi_kappa,jacobi_niter,source_pos[0]);
	
	//invert
	tinv-=take_time();
	inv_Q2_cg(solDD,source,NULL,conf,kappa,mass,nitermax,1,residue);	
	tinv+=take_time();
	reconstruct_doublet(sol_reco[0],sol_reco[1],solDD,conf,kappa,mass);
	
	for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
	  for(int ivol=0;ivol<loc_vol;ivol++)
	    {
	      //put the anti-periodic condition on the propagator
	      int dt=glb_coord_of_loclx[ivol][0]-source_pos[0];
	      double arg=M_PI*dt/glb_size[0];
	      complex phase={cos(arg),sin(arg)};
	      spincolor temp;
	      
	      unsafe_spincolor_prod_complex(temp,sol_reco[r][ivol],phase);
	      
	      put_spincolor_into_su3spinspin(S0[r][ivol],temp,id_sour,ic_sour);
	    }
      }

  if(rank==0) printf("inversions finished\n");
  
  //put the (1+ig5)/sqrt(2) factor
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ivol=0;ivol<loc_vol;ivol++)
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  rotate_spinspin_to_physical_basis(S0[r][ivol][ic1][ic2],!r,!r);
  
  if(rank==0) printf("rotations performed\n");
  
  free(temp_source);
}

//Calculate the proton contraction for a single point
void point_proton_contraction(spinspin contr,su3spinspin SU,su3spinspin SD,dirac_matr gamma1,dirac_matr gamma2,dirac_matr gamma3,dirac_matr gamma4)
{
  memset(contr,0,sizeof(spinspin));
  
  for(int ga1=0;ga1<4;ga1++)
    for(int ga2=0;ga2<4;ga2++)
      if(Proj[0][ga1][ga2][0]!=0||Proj[0][ga1][ga2][1]!=0)
	{
	  int delta1=gamma2.pos[ga1];
	  int delta2=gamma4.pos[ga2];
	  
	  for(int a1=0;a1<3;a1++)
	    for(int b1=0;b1<3;b1++)
	      for(int c1=0;c1<3;c1++)
		if(epsilon[a1][b1][c1])
		  for(int a2=0;a2<3;a2++)
		    for(int b2=0;b2<3;b2++)
		      for(int c2=0;c2<3;c2++)
			if(epsilon[a2][b2][c2])
			  {
			    for(int al1=0;al1<4;al1++)
			      {
				complex ter1={0,0};
				
				for(int al2=0;al2<4;al2++)
				  {
				    int be1=gamma1.pos[al1];
				    int be2=gamma3.pos[al2];
				    
				    complex ter2;
				    unsafe_complex_prod(ter2,SU[a1][a2][al1][al2],SU[c1][c2][delta1][delta2]);
				    complex_subt_the_prod(ter2,SU[a1][c2][al1][delta2],SU[c1][a2][delta1][al2]);
				  
				    safe_complex_prod(ter2,SD[b1][b2][be1][be2],ter2);
				    complex_summ_the_prod(ter1,gamma3.entr[al2],ter2);
				  }
				
				int se=epsilon[a1][b1][c1]*epsilon[a2][b2][c2];
				if(se==1) complex_summ_the_prod(contr[ga1][ga2],gamma1.entr[al1],ter1);
				else      complex_subt_the_prod(contr[ga1][ga2],gamma1.entr[al1],ter1);
			      }
			  }
	  
	  complex temp;
	  unsafe_complex_prod(temp,contr[ga1][ga2],gamma2.entr[ga1]);
	  unsafe_complex_prod(contr[ga1][ga2],temp,gamma4.entr[ga2]);
	}
}

//calculate all the 2pts contractions
void calculate_all_2pts(char *path)
{
  //output file
  FILE *output=open_text_file_for_output(path);

  tcontr-=take_time();
  
  char pm_tag[2][2]={"+","-"};
  
  complex *loc_contr[3],*glb_contr[3];
  for(int nns=0;nns<3;nns++)
    {
      loc_contr[nns]=(complex*)malloc(sizeof(complex)*glb_size[0]);
      glb_contr[nns]=(complex*)malloc(sizeof(complex)*glb_size[0]);
    }
  
  spinspin ter;
  complex point_contr[3];

  dirac_matr o3,o4=base_gamma[0];
  dirac_prod(&o3,&gC,&base_gamma[5]);

  for(int icontr=0;icontr<nproton_2pt_contr;icontr++)
    {

      dirac_matr o1,o2=base_gamma[list_2pt_op2[icontr]];
      dirac_prod(&o1,&gC,&base_gamma[list_2pt_op1[icontr]]);
      
      for(int rlike=0;rlike<2;rlike++)
	for(int rdislike=0;rdislike<2;rdislike++)
	  {
	    
	    if(rank==0) fprintf(output," # Two point for rlike=%d, rdislike=%d\n",rlike,rdislike);
	    
	    //perform the proton contraction putting operators on the sink or on the source
	    for(int SS=0;SS<2;SS++)
	      {
		//reset output
		for(int nns=0;nns<3;nns++) memset(loc_contr[nns],0,sizeof(complex)*glb_size[0]);
		
		//local loop
		for(int loc_site=0;loc_site<loc_vol;loc_site++)
		  {
		    int glb_t=glb_coord_of_loclx[loc_site][0];
		    
		    if(SS==0) point_proton_contraction(ter,S0[rlike][loc_site],S0[rdislike][loc_site],o1,o2,o3,o4);
		    else point_proton_contraction(ter,S0[rlike][loc_site],S0[rdislike][loc_site],o3,o4,o1,o2);
		    
		    for(int nns=0;nns<3;nns++)
		      {
			trace_prod_spinspins(point_contr[nns],ter,Proj[nns]);
			complex_summ(loc_contr[nns][glb_t],loc_contr[nns][glb_t],point_contr[nns]);
		      }
		  }
		
		//final reduction
		for(int nns=0;nns<3;nns++) MPI_Reduce(loc_contr[nns],glb_contr[nns],2*glb_size[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		
		if(rank==0)
		  {
		    if(SS==0) fprintf(output," # %s%s-P5S0\n",gtag[list_2pt_op1[icontr]],gtag[list_2pt_op2[icontr]]);
		    else      fprintf(output," # P5S0-%s%s\n",gtag[list_2pt_op1[icontr]],gtag[list_2pt_op2[icontr]]);
		    for(int nns=0;nns<3;nns++)
		      {
			if(nns<2) fprintf(output,"# Contraction with (1%sg4)/2\n",pm_tag[nns]);
			else      fprintf(output,"# Contraction with (1+g4)_00(spinorial)/2\n");
			for(int tt=0;tt<glb_size[0];tt++)
			  {
			    int t=(tt+source_pos[0])%glb_size[0];
			    fprintf(output," %+016.16g\t%+016.16g\n",glb_contr[nns][t][0],glb_contr[nns][t][1]);
			  }
			fprintf(output,"\n");
		      }
		  }
	      }
	  }
    }
  
  tcontr+=take_time();

  if(rank==0)
    {
      printf("contractions finished\n");
      fclose(output);
    }

  for(int nns=0;nns<3;nns++)
    {    
      free(loc_contr[nns]);
      free(glb_contr[nns]);
    }
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
  
  for(int iconf=0;iconf<nconf;iconf++)
    {
      read_conf_and_put_antiperiodic(conf,conf_path[iconf],source_pos[0]);
      
      calculate_S0();
      calculate_all_2pts(out_path[iconf]);
    }

  tot_time+=take_time();

  ///////////////////////////////////////////

  if(rank==0)
    {
      printf("Total time: %g s\n",tot_time);
      printf("-inversion time: %g%s avg: %d s\n",tinv/tot_time*100,"%",(int)(tinv/nconf/12));
      printf("-contraction time: %g%s\n",tcontr/tot_time*100,"%");
    }

  close_appretto();

  return 0;
}
