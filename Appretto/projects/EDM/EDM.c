#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

//configuration
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

void summ_the_point_proton_contraction(complex contr_plus,complex contr_minus,su3spinspin SU,su3spinspin SD)
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
		  
		  //printf("SD %d %d %d %d %g %g\n",b1,b2,beta1,beta2,SD[b1][b2][beta1][beta2][0],SD[b1][b2][beta1][beta2][1]);
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
  complex loc_plus[glb_size[0]],loc_minus[glb_size[0]];
  
  memset(loc_plus,0,sizeof(complex)*glb_size[0]);
  memset(loc_minus,0,sizeof(complex)*glb_size[0]);

  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    {
      int glb_t=glb_coord_of_loclx[loc_site][0];
      
      summ_the_point_proton_contraction(loc_plus[glb_t],loc_minus[glb_t],SU[loc_site],SD[loc_site]);
    }
  
  for(int t=0;t<glb_size[0];t++)
    for(int ri=0;ri<2;ri++)
      {
	loc_plus[t][ri]*=0.5;
	loc_minus[t][ri]*=0.5;
      }
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
}


void initialize_EDM(char *input_path)
{
  open_input(input_path);

  // 1) Information about the gauge conf
  
  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  //Init the MPI grid 
  init_grid();
  //Allocate the gauge Conf
  conf=allocate_quad_su3(loc_vol+loc_bord,"conf");
  //Read the gauge conf
  char conf_path[1024];
  read_str_str("GaugeConfPath",conf_path,1024);
  read_local_gauge_conf(conf,conf_path);
  //Put border condition and communicate
  double theta[4]={1,0,0,0};
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);  
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
  int ivol=loclx_of_coord(gx);
  int isloc=1;

  int lx[4];

  memset(source,0,sizeof(spincolor)*loc_vol);

  for(int idir=0;idir<4;idir++)
    {
      lx[idir]=gx[idir]-proc_coord[idir]*loc_size[idir];
      isloc=isloc && (lx[idir]>=0 && lx[idir]<loc_size[idir]);
    }

  if(isloc)
    {
      if(id<2) source[ivol][id][ic][0]=1;
      else source[ivol][id][ic][0]=-1;
    }
}      

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  initialize_EDM(arg[1]);
  
  ///////////////////////////////////////////

  for(int ic=0;ic<3;ic++)
    {
      for(int id=0;id<4;id++)
	{
	  create_point_source(source,ic,id,source_pos);
	  inv_Q2_cg(solDD,source,NULL,conf,kappa,mass,nitermax,1,residue);
	  reconstruct_doublet(sol_reco[0],sol_reco[1],solDD,conf,kappa,mass);

	  for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
	    for(int ivol=0;ivol<loc_vol;ivol++)
	      put_spincolor_into_su3spinspin(S0[r][ivol],sol_reco[r][ivol],id,ic);
	}
    }

  //put the (1+ig5)/sqrt(2) factor
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ivol=0;ivol<loc_vol;ivol++)
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  rotate_spinspin_to_physical_basis(S0[r][ivol][ic1][ic2],!r,!r);
  
  complex contr_plus[glb_size[0]],contr_minus[glb_size[0]];

  proton_contraction(contr_plus,contr_minus,S0[0],S0[1]);

  FILE *fout=open_text_file_for_output("proton");
  for(int t=0;t<glb_size[0];t++)
    if(rank==0)
      fprintf(fout,"%g %g\t%g %g\n",contr_plus[t][0],contr_plus[t][1],contr_minus[t][0],contr_minus[t][1]);
  if(rank==0) fclose(fout);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
