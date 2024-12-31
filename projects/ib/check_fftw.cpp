#include <nissa.hpp>

using namespace nissa;

void in_main(int narg,char **arg)
{
    
  glb_size[0]=glb_size[1]=glb_size[2]=2;
  glb_size[3]=6;
  init_grid(-1,-1);
  
  start_loc_rnd_gen(1324);
  
  complex *vin=nissa_malloc("vin",loc_vol,complex);
  complex *vout=nissa_malloc("vout",loc_vol,complex);
  for(int ivol=0;ivol<loc_vol;ivol++)
    rnd_get_Z4(vin[ivol],loc_rnd_gen+ivol);
  
  int ncpp=1;
  int sign=1;
  int normalize=0;
  fft4d(vout,vin,all_dirs,ncpp,sign,normalize);
  fft4d(vout,vout,all_dirs,ncpp,-sign,!normalize);
  
  double_vector_subtassign((double*)vout,(double*)vin,loc_vol*2);
  double nvon=double_vector_glb_norm2(vout,2);
  double nvin=double_vector_glb_norm2(vin,2);
  
  MASTER_PRINTF("norms: %lg %lg\n",nvon,nvin);
  // for(int ivol=0;ivol<loc_vol;ivol++)
  //   {
  //     MASTER_PRINTF("el[%d][RE]: %lg %lg\n",ivol,vin[ivol][RE],vout[ivol][RE]);
  //     MASTER_PRINTF("el[%d][IM]: %lg %lg\n",ivol,vin[ivol][IM],vout[ivol][IM]);
  //   }
  nissa_free(vin);
  nissa_free(vout);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}

