//we are calculate the NEUTRON EDM, so the ddu, and then the 110

#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

//configuration
int nconf;
char **conf_path,**outp_path;
quad_su3 *conf;
double mass;
double kappa;

//source
int source_pos[4];
spincolor *source;
su3spinspin *original_source;

//the propagators
su3spinspin *S0[2];

//inverter
int nitermax;
double residue;
spincolor *solDD,*sol_reco[2];

spinspin Proj[2]; //projectors over N and N*
spinspin C5; //C*gamma5

//                      e_00x   e_01x    e_02x     e_10x    e_11x   e_12x     e_20x   e_21x    e_22x
int epsilon[3][3][3]={{{0,0,0},{0,0,1},{0,-1,0}},{{0,0,-1},{0,0,0},{1,0,0}},{{0,1,0},{-1,0,0},{0,0,0}}};

//timings
double tinv=0,tcontr=0,tot_time=0;

//output file
FILE *output;

//sequential source for the "like" case, traced with Proj[0]
void prepare_like_sequential_source(su3spinspin *seq,su3spinspin *S)
{
  memset(seq,0,sizeof(su3spinspin)*loc_vol);
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int a1=0;a1<3;a1++)
      for(int z=0;z<3;z++)
	for(int c1=0;c1<3;c1++)
	  if(epsilon[a1][z][c1])
	    for(int a2=0;a2<3;a2++)
	      for(int x=0;x<3;x++)
		for(int c2=0;c2<3;c2++)
		  if(epsilon[a2][x][c2])
		    for(int al1=0;al1<4;al1++)
		      for(int al2=0;al2<4;al2++)
			for(int tau=0;tau<4;tau++)
			  for(int rho=0;rho<4;rho++)
			    if((C5[al2][tau][0]||C5[al2][tau][1])&&(C5[al1][rho][0]||C5[al1][rho][1]))
			      for(int ga1=0;ga1<4;ga1++)
				for(int ga2=0;ga2<4;ga2++)
				  {
				    int se=epsilon[a1][x][c1]*epsilon[a2][z][c2];
				    
				    complex ter;
				    
				    unsafe_complex_conj_conj_prod(ter,S[ivol][a1][a2][al1][al2],S[ivol][c1][c2][ga1][ga2]);
				    complex_subt_the_conj_conj_prod(ter,S[ivol][a1][c2][al1][ga2],S[ivol][c1][a2][ga1][al2]);
				    
				    safe_complex_conj1_prod(ter,C5[al1][rho],ter);
				    safe_complex_conj1_prod(ter,C5[al2][tau],ter);
				    
				    safe_complex_prod(ter,Proj[0][ga1][ga2],ter);
				    
				    for(int ri=0;ri<2;ri++)
				      if(se==1) seq[ivol][x][z][tau][rho][ri]+=ter[ri];
				      else      seq[ivol][x][z][rho][rho][ri]-=ter[ri];
				  }  
}

void initialize_EDM(char *input_path)
{
  //C5
  complex ima={0,1};
  dirac_matr gC,migC,gC5;
  memset(C5,0,sizeof(spinspin));
  dirac_prod(&migC,&(base_gamma[2]),&(base_gamma[4]));
  unsafe_dirac_compl_prod(&gC,&migC,ima);
  dirac_prod(&gC5,&gC,&(base_gamma[5]));
  
  for(int id1=0;id1<4;id1++)
    {
      int id2=gC5.pos[id1];
      for(int ri=0;ri<2;ri++)
	C5[id1][id2][ri]=gC5.entr[id1][ri];
    }
  
  //Proj[0] and Proj[1]
  for(int nns=0;nns<2;nns++) memset(Proj[nns],0,sizeof(spinspin));
  for(int id1=0;id1<4;id1++)
    {
      int id2=base_gamma[4].pos[id1];
      
      Proj[0][id1][id1][0]=Proj[1][id1][id1][0]=0.5;
      complex_prod_with_real(Proj[0][id1][id2],base_gamma[4].entr[id1],+0.5);
      complex_prod_with_real(Proj[1][id1][id2],base_gamma[4].entr[id1],-0.5);
    }

  open_input(input_path);

  // 1) Information about the gauge conf
  
  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  //Init the MPI grid 
  init_grid();
  //Allocate the gauge Conf
  conf=allocate_quad_su3(loc_vol+loc_bord,"conf");
  //Read the gauge conf
  read_str_int("NGaugeConf",&nconf);
  conf_path=(char**)malloc(sizeof(char*)*nconf);
  outp_path=(char**)malloc(sizeof(char*)*nconf);
  for(int iconf=0;iconf<nconf;iconf++)
    {
      conf_path[iconf]=(char*)malloc(1024);
      outp_path[iconf]=(char*)malloc(1024);
      read_str_str("GaugeConfPath",conf_path[iconf],1024);
      read_str(outp_path[iconf],1024);
    }
  //Put border condition and communicate
  //double theta[4]={1,0,0,0};

  //put_boundaries_conditions(conf,theta,1,0);
  //Kappa
  read_str_double("Kappa",&(kappa));
  
  // 2) Source position and mass
  read_str_int("SourcePosition",&(source_pos[0]));
  for(int i=1;i<4;i++) read_int(&(source_pos[i]));
  read_str_double("Mass",&mass);

  //Residue
  read_str_double("Residue",&residue);
  //Number of iterations
  read_str_int("NiterMax",&nitermax);

  close_input();

  /////////////////////////////// Allocate the various spinors //////////////////
  
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
  read_local_gauge_conf(conf,conf_path);
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    if(glb_coord_of_loclx[ivol][0]==(tsource+glb_size[0]-1)%glb_size[0])
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  for(int ri=0;ri<2;ri++)
	    conf[ivol][0][ic1][ic2][ri]*=-1;
  communicate_gauge_borders(conf);  
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
  for(int ic_sour=0;ic_sour<3;ic_sour++)
    for(int id_sour=0;id_sour<4;id_sour++)
      { //take the source and put g5
	for(int ivol=0;ivol<loc_vol;ivol++)
	  {
	    get_spincolor_from_su3spinspin(source[ivol],original_source[ivol],id_sour,ic_sour);
	    for(int id_sink=2;id_sink<4;id_sink++)
	      for(int ic_sink=0;ic_sink<3;ic_sink++)
		for(int ri=0;ri<2;ri++)
		  source[ivol][id_sink][ic_sink][ri]=-source[ivol][id_sink][ic_sink][ri];
	  }
	
	tinv-=take_time();
	inv_Q2_cg(solDD,source,NULL,conf,kappa,mass,nitermax,1,residue);	
	tinv+=take_time();
	reconstruct_doublet(sol_reco[0],sol_reco[1],solDD,conf,kappa,mass);
	
	for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
	  for(int ivol=0;ivol<loc_vol;ivol++)
	    put_spincolor_into_su3spinspin(S0[r][ivol],sol_reco[r][ivol],id_sour,ic_sour);
      }
  
  if(rank==0) printf("inversions finished\n");
  
  //put the (1+ig5)/sqrt(2) factor
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ivol=0;ivol<loc_vol;ivol++)
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  rotate_spinspin_to_physical_basis(S0[r][ivol][ic1][ic2],!r,!r);
  
  if(rank==0) printf("rotations performed\n");
}

//Calculate the proton contraction for a single point
void point_proton_contraction(spinspin contr,su3spinspin SU,su3spinspin SD)
{
  memset(contr,0,sizeof(spinspin));
  
  for(int a1=0;a1<3;a1++)
    for(int b1=0;b1<3;b1++)
      for(int c1=0;c1<3;c1++)
	if(epsilon[a1][b1][c1])
	  for(int a2=0;a2<3;a2++)
	    for(int b2=0;b2<3;b2++)
	      for(int c2=0;c2<3;c2++)
		if(epsilon[a2][b2][c2])
		  for(int al1=0;al1<4;al1++)
		    for(int al2=0;al2<4;al2++)
		      for(int be1=0;be1<4;be1++)
			for(int be2=0;be2<4;be2++)
			  if((C5[al1][be1][0]||C5[al1][be1][1])&&(C5[al2][be2][0]||C5[al2][be2][1]))
			    for(int ga1=0;ga1<4;ga1++)
			      for(int ga2=0;ga2<4;ga2++)
				{
				  int se=epsilon[a1][b1][c1]*epsilon[a2][b2][c2];
				  
				  complex ter;
				  unsafe_complex_prod(ter,SU[a2][a1][al2][al1],SU[c2][c1][ga2][ga1]);
				  complex_subt_the_prod(ter,SU[a2][c1][al2][ga1],SU[c2][a1][ga2][al1]);
				  
				  safe_complex_prod(ter,SD[b2][b1][be2][be1],ter);
				  safe_complex_prod(ter,C5[al1][be1],ter);
				  safe_complex_prod(ter,C5[al2][be2],ter);
				  
				  for(int ri=0;ri<2;ri++)
				    if(se==1) contr[ga2][ga1][ri]+=ter[ri];
				    else      contr[ga2][ga1][ri]-=ter[ri];
				}
}

//calculate all the 2pts contractions
void calculate_all_2pts()
{
  tcontr-=take_time();
  
  char ftag[2][1]={"U","D"},pm_tag[2][1]={"+","-"};
  
  complex *loc_contr[2],*glb_contr[2];
  for(int nns=0;nns<2;nns++)
    {
      loc_contr[nns]=(complex*)malloc(sizeof(complex)*glb_size[0]);
      glb_contr[nns]=(complex*)malloc(sizeof(complex)*glb_size[0]);
    }
  
  spinspin ter;
  complex point_contr[2];
  
  for(int r1=0;r1<2;r1++)
    {
      int r2=!r1;
      
      for(int nns=0;nns<2;nns++) memset(loc_contr[nns],0,sizeof(complex)*glb_size[0]);
  
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	{
	  int glb_t=glb_coord_of_loclx[loc_site][0];
	  
	  point_proton_contraction(ter,S0[r1][loc_site],S0[r2][loc_site]);
	  
	  for(int nns=0;nns<2;nns++)
	    {
	      trace_prod_spinspins(point_contr[nns],ter,Proj[nns]);
	      complex_summ(loc_contr[nns][glb_t],loc_contr[nns][glb_t],point_contr[nns]);
	    }
	}
      
      MPI_Barrier(MPI_COMM_WORLD);
      
      for(int nns=0;nns<2;nns++) MPI_Reduce(loc_contr[nns],glb_contr[nns],2*glb_size[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      
      if(rank==0)
	{
	  fprintf(output," # Two point for %s%s%s\n",ftag[r1],ftag[r1],ftag[r2]);
	  for(int nns=0;nns<2;nns++)
	    {
	      fprintf(output,"# Contraction with (1%sg4)/2\n",pm_tag[nns]);
	      for(int tt=0;tt<glb_size[0];tt++)
		{
		  int t=(tt+source_pos[0])%glb_size[0];
		  fprintf(output,"%g %g\n",glb_contr[nns][t][0],glb_contr[nns][t][1]);
		}
	    }
	}
    }

  tcontr+=take_time();
  
  if(rank==0) printf("contractions ultimated\n");

  for(int nns=0;nns<2;nns++)
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
  
  initialize_EDM(arg[1]);
  
  ///////////////////////////////////////////
  
  prepare_source();
  
  for(int iconf=0;iconf<nconf;iconf++)
    {
      output=open_text_file_for_output(outp_path[iconf]);
      
      read_conf_and_put_antiperiodic(conf,conf_path[iconf],source_pos[0]);
      
      calculate_S0();
      calculate_all_2pts();
      
      if(rank==0) fclose(output);
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
