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

//the propagators
su3spinspin *S0[2];

//inverter
int nitermax;
double residue;
spincolor *solDD,*sol_reco[2];

spinspin Pp,Pm; //projectors over N and N*
spinspin C5; //C*gamma5

void point_proton_contraction(spinspin contr,su3spinspin SU,su3spinspin SD)
{//                        e_00   e_01     e_02      e_10     e_11    e_12      e_20    e_21     e_22
  int epsilon[3][3][3]={{{0,0,0},{0,0,1},{0,-1,0}},{{0,0,-1},{0,0,0},{1,0,0}},{{0,1,0},{-1,0,0},{0,0,0}}};
  
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
			  //if((C5[al1][be1][0]||C5[al1][be1][1])&&(C5[al2][be2][0]||C5[al2][be2][1]))
			    for(int ga1=0;ga1<4;ga1++)
			      for(int ga2=0;ga2<4;ga2++)
				{
				  int se=epsilon[a1][b1][c1]*epsilon[a2][b2][c2];
				  
				  complex ter;
				  unsafe_complex_prod(ter,SU[a2][a1][al2][al1],SU[c2][c1][ga2][ga1]);
				  complex_subt_the_prod(ter,SU[a2][c1][al2][ga1],SU[c2][a1][ga2][al1]);
				  
				  safe_complex_prod(ter,SD[be2][be1][b2][b1],ter);
				  safe_complex_prod(ter,C5[al1][be1],ter);
				  safe_complex_prod(ter,C5[al2][be2],ter);
				  
				  for(int ri=0;ri<2;ri++)
				    if(se==1) contr[ga1][ga2][ri]=ter[ri];
				    else      contr[ga1][ga2][ri]=-ter[ri];
				}
}

void summ_the_point_proton_contraction_wrong(complex contr_plus,complex contr_minus,su3spinspin SU,su3spinspin SD)
{
  const int e_sign[6]={1,-1,-1,1,1,1-1};
  const int e_ind[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
  const int cc_p[4]={1,0,3,2};
  const int cc_e[4]={-1,1,-1,1};

  complex contr0={0,0},contr4={0,0};

  for(int iperm1=0;iperm1<6;iperm1++)
    {
      int a1=e_ind[iperm1][0];
      int b1=e_ind[iperm1][1];
      int c1=e_ind[iperm1][2];

      int s1=e_sign[iperm1];

      for(int iperm2=0;iperm2<6;iperm2++)
	{
	  int a2=e_ind[iperm2][0];
	  int b2=e_ind[iperm2][1];
	  int c2=e_ind[iperm2][2];

	  int s2=e_sign[iperm2];

	  for(int alpha1=0;alpha1<4;alpha1++)
	    {
	      int beta1=cc_p[alpha1];
	      for(int alpha2=0;alpha2<4;alpha2++)
		{
		  int beta2=cc_p[alpha2];
		  complex F0={0,0},F4={0,0};
		  for(int gamma1=0;gamma1<4;gamma1++)
		    {
		      int gamma2_0=gamma1;
		      int gamma2_4=(gamma1+2)%4;

		      complex_summ_the_prod(F0,SU[a1][a2][alpha1][alpha2],SU[c1][c2][gamma1][gamma2_0]);
		      complex_subt_the_prod(F0,SU[a1][c2][alpha1][gamma2_0],SU[c1][a2][gamma1][alpha2]);

		      complex_summ_the_prod(F4,SU[a1][a2][alpha1][alpha2],SU[c1][c2][gamma1][gamma2_4]);
		      complex_subt_the_prod(F4,SU[a1][c2][alpha1][gamma2_4],SU[c1][a2][gamma1][alpha2]);
		    }
		  
		  int fact=s1*s2*cc_e[alpha1]*cc_e[alpha2];
		  
		  if(fact==1)
		    {
		      complex_summ_the_prod(contr0,SD[b1][b2][beta1][beta2],F0);
		      complex_summ_the_prod(contr4,SD[b1][b2][beta1][beta2],F4);
		    }
		  else
		    {
		      complex_subt_the_prod(contr0,SD[b1][b2][beta1][beta2],F0);
		      complex_subt_the_prod(contr4,SD[b1][b2][beta1][beta2],F4);
		    }
		}
	    }
	}
    }
  
  complex_summ(contr_plus,contr_plus,contr0);
  complex_summ(contr_plus,contr_plus,contr4);

  complex_summ(contr_minus,contr_minus,contr0);
  complex_subt(contr_minus,contr_minus,contr4);
}

void proton_contraction(complex *contr_plus,complex *contr_minus,su3spinspin *SU,su3spinspin *SD)
{
  complex *loc_plus=(complex*)malloc(sizeof(complex)*glb_size[0]);
  complex *loc_minus=(complex*)malloc(sizeof(complex)*glb_size[0]);

  memset(loc_plus,0,sizeof(complex)*glb_size[0]);
  memset(loc_minus,0,sizeof(complex)*glb_size[0]);

  spinspin ter;
  complex cp,cm;
  
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    {
      int glb_t=glb_coord_of_loclx[loc_site][0];
      
      point_proton_contraction(ter,SU[loc_site],SD[loc_site]);
      
      trace_prod_spinspins(cp,ter,Pp);
      trace_prod_spinspins(cm,ter,Pm);

      complex_summ(loc_plus[glb_t],loc_plus[glb_t],cp);
      complex_summ(loc_minus[glb_t],loc_minus[glb_t],cm);
    }

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0) printf("going to perform the reduction\n");  
  
  if(rank_tot>0)
    {
      MPI_Reduce(loc_plus,contr_plus,2*glb_size[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(loc_minus,contr_minus,2*glb_size[0],MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    }
  else
    {
      memcpy(contr_plus,loc_plus,glb_size[0]*sizeof(complex));
      memcpy(contr_minus,loc_minus,glb_size[0]*sizeof(complex));
    }
  if(rank==0) printf("reduction completed\n");

  free(loc_plus);
  free(loc_minus);
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
  
  //Pp and Pm
  memset(Pp,0,sizeof(spinspin));
  memset(Pm,0,sizeof(spinspin));
  for(int id1=0;id1<4;id1++)
    {
      int id2=base_gamma[4].pos[id1];
      
      Pp[id1][id1][0]=0.5;
      complex_prod_with_real(Pp[id1][id2],base_gamma[4].entr[id1],+0.5);
      
      Pm[id1][id1][0]=0.5;
      complex_prod_with_real(Pm[id1][id2],base_gamma[4].entr[id1],-0.5);
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
  
  source=allocate_spincolor(loc_vol+loc_bord,"source");
  solDD=allocate_spincolor(loc_vol+loc_bord,"solDD");

  sol_reco[0]=allocate_spincolor(loc_vol,"solution_reco[0]");
  sol_reco[1]=allocate_spincolor(loc_vol,"solution_reco[1]");

  S0[0]=allocate_su3spinspin(loc_vol,"S0[0]");
  S0[1]=allocate_su3spinspin(loc_vol,"S0[1]");
}

void create_point_source(spincolor *source,int ic,int id,int *gx)
{
  int isloc=1;

  int lx[4];

  memset(source,0,sizeof(spincolor)*loc_vol);

  for(int idir=0;idir<4;idir++)
    {
      lx[idir]=gx[idir]-proc_coord[idir]*loc_size[idir];
      isloc=isloc && (lx[idir]>=0 && lx[idir]<loc_size[idir]);
    }

  int ivol=loclx_of_coord(lx);

  if(isloc)
    {
      if(id<2) source[ivol][id][ic][0]=1;
      else source[ivol][id][ic][0]=-1;
    }
}      

int main(int narg,char **arg)
{
  double tinv=0,tcontr=0,tot_time=0;
 
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
  
  complex *contr_plus=(complex*)malloc(sizeof(complex)*glb_size[0]);
  complex *contr_minus=(complex*)malloc(sizeof(complex)*glb_size[0]);
      
  ///////////////////////////////////////////

  for(int iconf=0;iconf<nconf;iconf++)
    {
      read_local_gauge_conf(conf,conf_path[iconf]);

      for(int ivol=0;ivol<loc_vol;ivol++)
	if(glb_coord_of_loclx[ivol][0]==(source_pos[0]+glb_size[0]-1)%glb_size[0])
	  for(int ic1=0;ic1<3;ic1++)
	    for(int ic2=0;ic2<3;ic2++)
	      for(int ri=0;ri<2;ri++)
		conf[ivol][0][ic1][ic2][ri]*=-1;
      communicate_gauge_borders(conf);  
      
      for(int ic=0;ic<3;ic++)
	{
	  for(int id=0;id<4;id++)
	    {
	      create_point_source(source,ic,id,source_pos);
	      tinv-=take_time();
	      inv_Q2_cg(solDD,source,NULL,conf,kappa,mass,nitermax,1,residue);	
	      tinv+=take_time();
	      reconstruct_doublet(sol_reco[0],sol_reco[1],solDD,conf,kappa,mass);
	      
	      for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
		for(int ivol=0;ivol<loc_vol;ivol++)
		  put_spincolor_into_su3spinspin(S0[r][ivol],sol_reco[r][ivol],id,ic);
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
      
      tcontr-=take_time();
      proton_contraction(contr_plus,contr_minus,S0[0],S0[1]);
      tcontr+=take_time();
      
      if(rank==0) printf("contractions ultimated\n");
      
      FILE *fout=open_text_file_for_output(outp_path[iconf]);
      for(int tt=0;tt<glb_size[0];tt++)
	if(rank==0)
	  {
	    int t=(tt+source_pos[0])%glb_size[0];
	    fprintf(fout,"%g %g\t%g %g\n",contr_plus[t][0],contr_plus[t][1],contr_minus[t][0],contr_minus[t][1]);
	  }
      if(rank==0) fclose(fout);
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
