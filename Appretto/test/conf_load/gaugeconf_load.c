#include "appretto.h"

void test(quad_su3 *conf,as2t tmunu)
{
  master_printf("plaq: %.10g\n",global_plaquette(conf));
  
  as2t_su3 *clov=appretto_malloc("clov",loc_vol,as2t_su3);
  
  Pmunu_term(clov,conf);
  
  double partial_summ=0;
  complex temp;
  as2t as;
  for(int X=0;X<loc_vol;X++)
    {
      for(int ias=0;ias<6;ias++)
	{
	  as[ias][0]=clov[X][ias][0][0][0];
	  as[ias][1]=clov[X][ias][0][0][1];
	}
      as2t_saturate(temp,tmunu,as);
      partial_summ+=temp[0];
    }
  appretto_free(clov);
  double total_leaves_summ;
  MPI_Reduce(&partial_summ,&total_leaves_summ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  master_printf("%.12g\n",total_leaves_summ);
}

int main(int narg,char **arg)
{
  char filename[1024];

  //basic mpi initialization
  init_appretto();

  if(narg<2) crash("Use: %s input_file\n",arg[0]);

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  read_str_str("Filename",filename,1024);

  close_input();

  //Init the MPI grid 
  init_grid();

  ///////////////////////////////////////////

  quad_su3 *conft=appretto_malloc("conft",loc_vol+loc_bord+loc_edge,quad_su3);
  quad_su3 *conf=appretto_malloc("conf",loc_vol+loc_bord+loc_edge,quad_su3);
  {
    dirac_matr temp_gamma[4];
    temp_gamma[0]=base_gamma[4];
    temp_gamma[1]=base_gamma[1];
    temp_gamma[2]=base_gamma[2];
    temp_gamma[3]=base_gamma[3];
    
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
    read_gauge_conf(conft,filename);  
    communicate_gauge_borders(conft);
    communicate_gauge_edges(conft);
    test(conft,tmunu);
  }

  if(rank_tot==1)
    {
      printf("Rotating\n");
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
      ac_rotate_gauge_conf(conf,conft,1);
      communicate_gauge_borders(conf);
      communicate_gauge_edges(conf);
      test(conf,tmunu);
    }

  appretto_free(conf);
  appretto_free(conft);
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
