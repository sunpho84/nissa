#include "appretto.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();
  
  glb_size[0]=8;
  glb_size[1]=4;
  
  //Init the MPI grid 
  init_grid();
  
  spincolor *in=appretto_malloc("in",loc_vol,spincolor);
  spincolor *out=appretto_malloc("out",loc_vol,spincolor);

  start_loc_rnd_gen(121);
  //generate_undiluted_source(in,RND_Z4,-1);
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	{
	  in[ivol][id][ic][0]=glb_coord_of_loclx[ivol][1]*100000+glb_coord_of_loclx[ivol][2]*10000+glb_coord_of_loclx[ivol][3]*1000+glb_coord_of_loclx[ivol][0]*100+id*10+ic;
	  in[ivol][id][ic][1]=0;
	}
  
  fft4d((complex*)out,(complex*)in,12,1);
  spincolor *ori=appretto_malloc("ori",loc_vol,spincolor);
  fft4d((complex*)ori,(complex*)out,12,-1);

  double loc_diff=0,loc_norm=0;
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
      {
	for(int mu=0;mu<4;mu++) printf("%d ",glb_coord_of_loclx[ivol][mu]);
	printf("| %lg %lg ",in[ivol][id][ic][0],in[ivol][id][ic][1]);
	printf("| %lg %lg ",out[ivol][id][ic][0],out[ivol][id][ic][1]);
	printf("| %lg %lg ",ori[ivol][id][ic][0],ori[ivol][id][ic][1]);
	printf("| %lg %lg\n",ori[ivol][id][ic][0]-in[ivol][id][ic][0],ori[ivol][id][ic][1]-in[ivol][id][ic][1]);
	loc_diff+=pow(in[ivol][id][ic][0]-ori[ivol][id][ic][0],2)+pow(in[ivol][id][ic][1]-ori[ivol][id][ic][1],2);
	loc_norm+=pow(in[ivol][id][ic][0],2)+pow(in[ivol][id][ic][1],2);
      }
  double glb_diff,glb_norm;
  MPI_Allreduce(&loc_diff,&glb_diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&loc_norm,&glb_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  master_printf("%lg %lg\n",glb_diff,glb_norm);
  
  appretto_free(ori);
  appretto_free(out);
  appretto_free(in);
  
  close_appretto();
  
  return 0;
}
