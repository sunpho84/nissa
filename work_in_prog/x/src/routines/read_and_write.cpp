#include <mpi.h>

#include "../../../../src/nissa.h"

#include "../types/types.h"

void write_corr16(char *path,corr16 *v,int prec)
{
  //Open the file
  ILDG_File file=ILDG_File_open_for_write(path);
  
  //Write the binary data
  write_double_vector(file,(double*)v,sizeof(corr16)/sizeof(double),prec,"scidac-binary-data");

  ILDG_File_close(file);
}

void read_corr16(corr16 *v,char *path)
{read_real_vector((double*)v,path,"scidac-binary-data",sizeof(corr16)/sizeof(double));}
