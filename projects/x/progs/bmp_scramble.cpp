#include <math.h>

#include "nissa.hpp"
using namespace nissa;
using namespace std;

#include "../src/routines/bmp.hpp"

color *image;
bmpfile bmp;

//initialize the program
void init_prog(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<3) CRASH("use %s filein fileout",arg[0]);
  
  //read the file
  read_bmpfile(bmp,arg[1]);
  
  glb_size[0]=bmp.info_header.width;
  glb_size[1]=bmp.info_header.height;
  glb_size[2]=glb_size[3]=1;
  
  //init the grid
  init_grid(0,0);
  
  //alloc image
  image=nissa_malloc("image",loc_vol,color);
  vector_reset(image);
  
  //copy the image
  NISSA_LOC_VOL_LOOP(ivol)
    for(int ic=0;ic<3;ic++)
      image[ivol][ic][0]=(bmp.data)[ivol*3+ic];
}

//close the program
void close_calc()
{
  nissa_free(bmp.data);
  nissa_free(image);
  
  close_nissa();
}
/*
//save the image
void save(char *path)
{
  //copy the image
  NISSA_LOC_VOL_LOOP(ivol)
    for(int ic=0;ic<3;ic++)
      bmp.data[ivol*3+ic]=(unsigned int) 
}
*/
int main(int narg,char **arg)
{
  init_prog(narg,arg);
  
  close_calc();
  
  return 0;
}
