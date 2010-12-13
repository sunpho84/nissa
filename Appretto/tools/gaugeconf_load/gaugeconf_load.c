#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

void test(quad_su3 *conf,as2t tmunu)
{
  double gplaq; 

  gplaq=global_plaquette(conf);
  if(rank==0) printf("plaq: %.10g\n",gplaq);

  su3_as2t *clov=(su3_as2t*)malloc(sizeof(su3_as2t)*loc_vol);
  
  clover_term(clov,conf);
  
  double partial_summ=0;
  complex temp;
  for(int X=0;X<loc_vol;X++)
    {
      as2t_saturate(temp,tmunu,clov[X][0][0]);
      partial_summ+=temp[0];
    }
  free(clov);
  double total_leaves_summ;
  MPI_Reduce(&partial_summ,&total_leaves_summ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if(rank==0) printf("%.12g\n",total_leaves_summ);
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

  quad_su3 *conft=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+loc_bord+loc_edge));
  read_local_gauge_conf(conft,filename);  
  communicate_gauge_borders(conft);
  communicate_gauge_edges(conft);
  test(conft,smunu[0][0]);

  if(rank_tot==1)
    {
      dirac_matr temp_gamma[4];
      temp_gamma[0]=base_gamma[4];
      temp_gamma[1]=base_gamma[1];
      temp_gamma[2]=base_gamma[3];
      temp_gamma[3]=base_gamma[2];
      for(int d=0;d<4;d++)
	{
	  temp_gamma[2].entr[d][0]=-temp_gamma[2].entr[d][0];
	  temp_gamma[2].entr[d][1]=-temp_gamma[2].entr[d][1];
	}

      as2t tmunu;
      for(int d=0;d<6;d++) tmunu[d][0]=tmunu[d][1]=0;
      int d=0;
      for(int mu=0;mu<4;mu++)
	for(int nu=mu+1;nu<4;nu++)
	  {
	    dirac_matr te;
	    dirac_prod(&te,&(temp_gamma[mu]),&(temp_gamma[nu]));
	    if(te.pos[0]==0)
	      {
		tmunu[d][0]+=te.entr[0][0];
		tmunu[d][1]+=te.entr[0][1];
	      }
	    d++;
	  }
      quad_su3 *conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+loc_bord+loc_edge));
      ac_rotate_gaugeconf(conf,conft,1);
      communicate_gauge_borders(conf);
      communicate_gauge_edges(conf);
      test(conf,tmunu);
      free(conf);
    }

  free(conft);
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
