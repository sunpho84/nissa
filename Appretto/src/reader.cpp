#pragma once

#include <lemon.h>
#include <lime.h>
#include "revert_endianess.cpp"

using namespace std;

//Read a whole spincolor
void read_spincolor(char *path,spincolor *spinor)
{
  char *header_type=NULL;

  //Open the file
  MPI_File *reader_file=new MPI_File;
  MPI_File_open(cart_comm,path,MPI_MODE_RDONLY,MPI_INFO_NULL,reader_file);
  LemonReader *reader=lemonCreateReader(reader_file,cart_comm);

  while(lemonReaderNextRecord(reader)!=LIME_EOF)
    {
      header_type=lemonReaderType(reader);
      if(rank==0) cout<<header_type<<endl;
      if (strcmp("scidac-binary-data",header_type)==0) if(rank==0) cout<<"Trovato"<<endl;
    }

  lemonDestroyReader(reader);
  MPI_File_close(reader_file);
}
