#include "nissa.h"

void test_lx_spi_comm()
{
  //allocate the vector
  spincolor *vec=nissa_malloc("vec",loc_vol+bord_vol,spincolor);
  vector_reset(vec);
  nissa_loc_vol_loop(ivol)
    vec[ivol][0][0][0]=glblx_of_loclx[ivol];
  
  communicate_lx_spincolor_borders(vec);
  
  for(int ivol=loc_vol;ivol<loc_vol+bord_vol;ivol++)
    {
      int obt=vec[ivol][0][0][0];
      int exp=glblx_of_bordlx[ivol-loc_vol];
      if(obt!=exp) crash("for ivol %d expected %d obtained %d",ivol,exp,obt);
    }
  
  nissa_free(vec);
  
  master_printf("lx_test passed\n");
}

void test_eo_spi_comm()
{
  int eo=ODD;
  
  //allocate the vector
  color *vec=nissa_malloc("vec",loc_volh+bord_volh,color);
  vector_reset(vec);
  nissa_loc_volh_loop(loceo)
    {
      int loclx=loclx_of_loceo[ODD][loceo];
      vec[loceo][0][0]=glblx_of_loclx[loclx];
    }
  
  //allocate the communicator
  communicate_od_color_borders(vec);
  
  for(int ieo=loc_volh;ieo<loc_volh+bord_volh;ieo++)
    {
      int bordeo=ieo-loc_volh;
      int ilx=loclx_of_loceo[ODD][ieo];
      int obt=vec[ieo][0][0];
      int exp=glblx_of_bordlx[ilx-loc_vol];
      if(obt!=exp) crash("for ieo %d (bordeo %d) expected %d obtained %d\n",ieo,bordeo,exp,obt);
    }
  
  nissa_free(vec);
  
  master_printf("evn test passed\n");
}

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  init_grid(16,16);
  
  test_lx_spi_comm();
  test_eo_spi_comm();
    
  close_nissa();
  
  return 0;
}
