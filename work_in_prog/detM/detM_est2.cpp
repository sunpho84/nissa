#include "nissa.h"

const double rad2=1.414213562373095048801688724209;

void stoc_source(spincolor *spinore)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
        {
	  spinore[ivol][id][ic][0]=pm_one(ivol)/rad2;
	  spinore[ivol][id][ic][1]=pm_one(ivol)/rad2;
	}
}

void scal_prod(complex c,spincolor *a,spincolor *b)
{
  complex loc_c={0,0};
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	complex_summ_the_conj1_prod(loc_c,a[ivol][id][ic],b[ivol][id][ic]);
  MPI_Allreduce(loc_c,c,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

int main(int narg,char **arg)
{
  double kappa,mass;
  int seed=0;
  char filename[1024];

  //basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L);

  read_str_int("Seed",&seed);
  read_str_double("Kappa",&kappa);
  read_str_double("Mass",&mass);
  read_str_str("Filename",filename,1024);
  
  close_input();

  //Initialize the random generator
  init_random(seed);

  ///////////////////////////////////////////

  quad_su3 *conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol));
  read_ildg_gauge_conf(conf,filename);
  
  spincolor *eta=allocate_spincolor(loc_vol+bord_vol,"eta");
  spincolor *source=allocate_spincolor(loc_vol+bord_vol,"source");
  spincolor *phi=allocate_spincolor(loc_vol+bord_vol,"phi");
  
  int ns=100;
  complex tr[ns];
  
  for(int is=0;is<ns;is++)
    {
      stoc_source(eta);
      for(int ivol=0;ivol<loc_vol;ivol++)
	for(int id=0;id<2;id++)
	  for(int ic=0;ic<3;ic++)
	    for(int ri=0;ri<2;ri++)
	      {
		source[ivol][id][ic][ri]=eta[ivol][id][ic][ri];
		source[ivol][id+2][ic][ri]=-eta[ivol][id+2][ic][ri];
	      }
      
      inv_Q2_cg(phi,source,NULL,conf,kappa,mass,1000,5,1.e-10);
      apply_Q(source,phi,conf,kappa,mass);
      
      scal_prod(tr[is],eta,source);
      printf("%d %lg %lg\n",is,tr[is][0],tr[is][1]);
    }
  
  
  return 0;
}
