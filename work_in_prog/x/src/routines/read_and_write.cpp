#include <lemon.h>
#include <mpi.h>

#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/IO/writer.h"
#include "../../../../src/IO/reader.h"
#include "../../../../src/base/global_variables.h"
#include "../../../../src/geometry/geometry_lx.h"
#include "../../../../src/base/debug.h"
#include "../../../../src/base/routines.h"
#include "../../../../src/base/vectors.h"

#include "../types/types.h"

void write_corr16(char *path,corr16 *v,int prec)
{
  //Open the file
  MPI_File *writer_file=nissa_malloc("Writer_file",1,MPI_File);
  int ok=MPI_File_open(cart_comm,path,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,writer_file);
  if(ok!=MPI_SUCCESS) crash("Couldn't open for writing the file: '%s'\n",path);
  
  MPI_File_set_size(*writer_file,0);
  LemonWriter *writer=lemonCreateWriter(writer_file,cart_comm);
  
  //Write the binary data
  write_double_vector(writer,(char*)v,"scidac-binary-data",sizeof(corr16)/sizeof(double),prec);

  //Close the file
  lemonDestroyWriter(writer);
  MPI_File_close(writer_file);
  nissa_free(writer_file);
}

void read_corr16(corr16 *v,char *path)
{read_real_vector((double*)v,path,"scidac-binary-data",sizeof(corr16)/sizeof(double));}
