#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

void density_profile(double *rho,spincolor *sp)
{
  int L=glb_size[1];

  int n[L];

  memset(rho,0,sizeof(double)*L);
  memset(n,0,sizeof(int)*L);

  for(int l=0;l<loc_vol;l++)
    if(glb_coord_of_loclx[l][0]==0)
      {
	double d=0;
	for(int mu=1;mu<4;mu++)
	  {
	    int x=glb_coord_of_loclx[l][mu];
	    if(x>=L/2) x-=L;
	    d+=x*x;
	  }
	d=sqrt(d);
	n[(int) d]++;
	
	for(int id=0;id<4;id++)
	  for(int ic=0;ic<3;ic++)
	    for(int ri=0;ri<2;ri++)
	      rho[(int)d]+=sp[l][id][ic][ri]*sp[l][id][ic][ri];
      }
  
  for(int d=0;d<L;d++) if(n[d]) rho[d]/=n[d];
}

int main(int narg,char **arg)
{
  char filename[1024];

  //basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  read_str_str("Filename",filename,1024);

  close_input();
  
  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////
  
  //allocate and read conf
  quad_su3 *orig_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"orig_conf");
  quad_su3 *smea_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"smea_conf");
  read_local_gauge_conf(orig_conf,filename);

  ape_smearing(smea_conf,orig_conf,1,0.1);
  printf("gauge conf smeared\n");

  //allocate and generate the source
  spincolor *origi_sp=allocate_spincolor(loc_vol+loc_bord,"orig_spincolor");
  memset(origi_sp,0,sizeof(spincolor)*loc_vol);
  if(rank==0) origi_sp[0][0][0][0]=1;
  
  //allocate and smeard
  spincolor *smear_sp=allocate_spincolor(loc_vol+loc_bord,"smear_spincolor");
  dina_smearing(smear_sp,origi_sp,smea_conf,4,50,0);
  
  int L=glb_size[1];
  double n_or[L],n_sm[L];
  density_profile(n_or,origi_sp);
  density_profile(n_sm,smear_sp);

  FILE *fout=fopen("/tmp/prof","w");
  for(int d=0;d<L;d++)
    {
      n_or[d]=sqrt(n_or[d]);
      n_sm[d]=sqrt(n_sm[d]);
      fprintf(fout,"%d %g %g\n",d,n_or[d],n_sm[d]);
    }
  fclose(fout);
    
  ///////////////////////////////////////////
  
  close_appretto();

  return 0;
}
