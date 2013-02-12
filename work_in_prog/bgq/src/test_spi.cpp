#include "nissa.h"

void test_spi_comm()
{
  spi_comm_t a;
  set_lx_spi_comm(a,sizeof(double));

  //fill in the surface the rank
  for(int ivol=0;ivol<bord_vol;ivol++)
    {
      int in=0;
      for(int mu=0;mu<4;mu++)
	{
	  int up=loclx_neighup[ivol+loc_vol][mu];
	  int dw=loclx_neighdw[ivol+loc_vol][mu];
	  if(up>=0&&up<loc_vol) in=up;
	  if(dw>=0&&dw<loc_vol) in=dw;
	}
      ((double*)a.send_buf)[ivol]=glblx_of_loclx[in];
      ((double*)a.recv_buf)[ivol]=-10;
    }

  spi_start_comm(a);
  spi_comm_wait(a);
  
  char path[100];
  sprintf(path,"conf_%d",rank);
  FILE *fout=fopen(path,"w");
  
  int ibord=0;
  for(int mu=0;mu<4;mu++)
    {
      int rup=rank_neighup[mu];
      int rdw=rank_neighdw[mu];
      
      for(int ivol=0;ivol<loc_vol;ivol++)
	{
	  int lup=loclx_neighup[ivol][mu];
	  int ldw=loclx_neighdw[ivol][mu];
	  
	  if(lup>=loc_vol)
	    {
	      int bup=lup-loc_vol;
	      fprintf(fout,"ivol %d, bord %d, up: %d, exp %d\n",ivol,bup,(int)(((double*)(a.recv_buf))[bup]),glblx_of_bordlx[bup]);
	    }
	  if(ldw>=loc_vol)
	    {
	      int bdw=ldw-loc_vol;
	      fprintf(fout,"ivol %d, bord %d, dw: %d, exp %d\n",ivol,bdw,(int)(((double*)(a.recv_buf))[bdw]),glblx_of_bordlx[bdw]);
	    }
	}
    }
  
  
  /*
  for(int ibord=0;ibord<bord_vol;ibord++)
    fprintf(fout,"ibord %d: %lg\n",ibord,((double*)(a.recv_buf))[ibord]);
  */
  fclose(fout);
  
  unset_spi_comm(a);
  
  master_printf("test passed\n");
}


int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  init_grid(16,16);
  
  init_spi();
  
  test_spi_comm();
  
  close_nissa();
  
  return 0;
}
