#include <mpi.h>

#include "../../../../src/nissa.h"

#include "../types/types.h"

//reorder a read corr16
void reorder_read_corr16(corr16 *c)
{
  int *order=nissa_malloc("order",loc_vol,int);

  int x[4];
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
        for(x[3]=0;x[3]<loc_size[3];x[3]++)
          {
            int isour=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
            int idest=loclx_of_coord(x);
            order[isour]=idest;
          }

  reorder_vector((char*)c,order,loc_vol,sizeof(corr16));
  nissa_free(order);
}

void read_corr16(corr16 *v,char *path)
{
  read_real_vector((double*)v,path,"scidac-binary-data",sizeof(corr16)/sizeof(double));
  reorder_read_corr16(v);
}

void write_corr16(char *path,corr16 *v,int prec=64)
{
  //Open the file
  ILDG_File file=ILDG_File_open_for_write(path);
  
  //order things as expected
  corr16 *temp=nissa_malloc("temp",loc_vol,corr16);

  int x[4],isour,idest;
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
        for(x[3]=0;x[3]<loc_size[3];x[3]++)
          {
            idest=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
            isour=loclx_of_coord(x);

            memcpy(temp[idest],v[isour],sizeof(corr16));
          }
  
  //Write the binary data
  write_double_vector(file,(double*)temp,sizeof(corr16)/sizeof(double),prec,"scidac-binary-data");

  nissa_free(temp);
  
  ILDG_File_close(file);
}

