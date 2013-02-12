#include "nissa.h"

void test_spi_comm()
{
  //allocate the vector
  spincolor *vec=nissa_malloc("vec",loc_vol+bord_vol,spincolor);
  vector_reset(vec);
  nissa_loc_vol_loop(ivol)
    vec[ivol][0][0][0]=glblx_of_loclx[ivol];
  
  //allocate the communicator
  spi_communicate_lx_borders(vec,spi_lx_spincolor_comm,sizeof(spincolor));
  
  for(int ivol=loc_vol;ivol<loc_vol+bord_vol;ivol++)
    {
      int obt=vec[ivol][0][0][0];
      int exp=glblx_of_bordlx[ivol-loc_vol];
      if(obt!=exp) crash("for ivol %d expected %d obtained %d",ivol,exp,obt);
    }
  
  nissa_free(vec);
  
  master_printf("test passed\n");
}


int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  init_grid(16,16);
  
  test_spi_comm();
  
  close_nissa();
  
  return 0;
}
